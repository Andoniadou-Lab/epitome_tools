import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scipy.sparse as sp
import scanpy as sc
import gc
import anndata
import polars as pl
from pathlib import Path
import scipy.io
import zipfile
import decoupler as dc
import shutil
from zipfile import ZipInfo

from concurrent.futures import ProcessPoolExecutor
anndata.settings.allow_write_nullable_strings = True


#####---------------------------------------------------------------------
##Creating large merged object from all pituitary datasets
#####---------------------------------------------------------------------

def remove_redundant_cols(adata):
  def drop_mask(cols):
      # True = drop: name contains '0' (but isn't exactly '0') or contains '-1'
      return np.array(
          [('0' in str(c) and str(c) != '0') or ('-1' in str(c)) for c in cols],
          dtype=bool,
      )

  adata.obs = adata.obs.loc[:, ~drop_mask(adata.obs.columns)]
  adata.var = adata.var.loc[:, ~drop_mask(adata.var.columns)]

  # rename obs from 1-len
  adata.obs.index = (np.arange(len(adata.obs)) + 1).astype(str)

  # turn to string so saving won't throw error
  adata.obs.columns = adata.obs.columns.astype(str)
  adata.var.columns = adata.var.columns.astype(str)

  for col in adata.obs.columns:
    adata.obs[col] = adata.obs[col].astype(str)

  adata.obsm = {}
  adata.layers = {}
  return adata


def remove_low_quality_clusters(adata, scvi_labels_to_keep, SRA_ID, df, ignore_core=False):

  author = df[df["SRA_ID"]==SRA_ID]["Author"].values[0]
  GEO = df[df["SRA_ID"]==SRA_ID]["GEO"].values[0]

  if ignore_core: #if not core sample, therefore not in scvi labels to keep, still keep the data by skipping filtering
    print(f"Ignoring core filtering for {SRA_ID}")
    return adata
  
  cell_ids = adata.obs['0'].astype(str) if '0' in adata.obs.columns else (adata.obs[0].astype(str) if 0 in adata.obs.columns else adata.obs_names.astype(str))

  if author == "Rebboah et al. (2025)":
    adata.obs["barcode"] = GEO + ":" + adata.obs["bc"].astype(str) + "_" + adata.obs["subpool"].astype(str) + "-1"
  elif author == "Sochodolsky et al. (2026)":
    adata.obs["barcode"] = SRA_ID + ":" + cell_ids + "-1"
  else:
    adata.obs["barcode"] = GEO + ":" + cell_ids + "-1"

  shape1 = adata.shape
  adata = adata[adata.obs["barcode"].isin(scvi_labels_to_keep[scvi_labels_to_keep["SRA_ID"]==SRA_ID]["barcode"].values)].copy()
  shape2 = adata.shape
  print(f"Removed {shape1[0]-shape2[0]} cells")

  return adata


def generate_merged_object(nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                           remaining_path="/analysis/adata_new_assigned_0828.h5ad",
                           back_up_remaining_path="/analysis/adata_assigned.h5ad",
                           atlas_source_path='/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
                           output_path="/content/drive/MyDrive/pituitary_atlas/Unstructured/wt_adata_merged_post_scvi_0902.h5ad",
                           filtering_based_on_barcodes = False,
                           scvi_labels_to_keep = None,
                           modality = ['sc','sn','multi_rna'],
                           normal = [1],
                           n_cells_threshold = 300,
                           ID_col='SRA_ID',
                           include_new_samples = False, #samples that are Core == 2
                           cell_type_col='new_cell_type',
                           target_cols_to_keep = [
                               'SRA_ID',
                               'GEO',
                               'Modality',
                               'Author',
                               'Name',
                               'Conditions',
                               'Normal',
                               'Background',
                               'Comp_sex',
                               'Estrous_stage',
                               'Age_numeric',
                               'Release'
                           ]
                           ):
    
    gc.collect()
    df = pd.read_excel(atlas_source_path)
    pituitary_atlas = df.copy()
    df = df[df['Modality'].isin(modality)]
    df = df[df["Normal"].isin(normal)]
    df = df[df["n_cells"]>n_cells_threshold]
    df.reset_index(drop=True, inplace=True)
    adata_merged_flag = False
    for i in range(len(df)-1,-1,-1):
      try:
          
          ID = df[ID_col][i]
          core_value = df[df[ID_col]==ID]["Core"].values[0]
          #convert to int
          core_value = int(core_value)
          print(ID)
          try:
            adata = sc.read(f"{nac_path}{ID}{remaining_path}")
          except:
            print("Falling back to adata_assigned")
            adata = sc.read(f"{nac_path}{ID}{back_up_remaining_path}")

          #ensure X is sparse!
          if not sp.issparse(adata.X):
              adata.X = sp.csr_matrix(adata.X)
          print("Loaded adata shape:", adata.shape)

          # 1. Add all metadata to obs
          adata.obs[ID_col] = ID

          for col in target_cols_to_keep:
              if col in df.columns:
                  adata.obs[col] = df[df[ID_col]==ID][col].values[0]

          # 2. Add orig_index
          adata.obs["orig_index"] = adata.obs.index.astype(str)

          # 2b. Keep the barcode so cells can be matched to the UMAP reference later
          if "0" in adata.obs.columns:
              adata.obs["cell_barcode"] = adata.obs["0"].astype(str)
          elif {"bc", "subpool"}.issubset(adata.obs.columns):  # Parse samples
              adata.obs["cell_barcode"] = adata.obs["bc"].astype(str) + "_" + adata.obs["subpool"].astype(str)
          else:
              adata.obs["cell_barcode"] = adata.obs_names.astype(str)

          core_action = False
          if include_new_samples and core_value != 0 and core_value != 1:
            core_action = True #true, therefore ignore filter and keep all cells
             
          if filtering_based_on_barcodes:
            adata = remove_low_quality_clusters(adata, scvi_labels_to_keep, ID, df,ignore_core=core_action)
          
          if not include_new_samples and core_value != 0 and core_value != 1:
            print(f"Not including non-core samples, removing {ID}")
            continue


          # add new_cell_type to target_cols_to_keep if it exists in adata.obs
          if cell_type_col in adata.obs.columns:
              if cell_type_col not in target_cols_to_keep:
                  target_cols_to_keep.append(cell_type_col)

          # 4. Subset safely using .copy() to decouple from original index structure
          adata.obs = adata.obs[target_cols_to_keep + ["cell_barcode"]].copy()

          # 5. Rename WITHOUT inplace=True
          adata.obs = adata.obs.rename(columns={cell_type_col: "cell_type"})
          print("Filtered obs to only keep relevant columns and renamed cell_type successfully")


          var_set = set(adata.var_names.astype(str))

          if "1" in var_set and "2" in var_set and "3" in var_set:
              print("The problematic sample is:", ID)
              break

          if adata_merged_flag == False:
              adata_merged = adata.copy()
              adata_merged_flag = True
          else:
              adata_merged = anndata.concat([adata_merged, adata], join="outer", uns_merge=None)
              adata_merged.obs_names_make_unique()
          print(adata_merged.shape)
          del adata

          #sleep 5 seconds
          gc.collect()
          time.sleep(2)
          

      #except print error
      except Exception as e:
        import traceback
        print(f"Error processing {ID}: {type(e).__name__} - {e}")
        traceback.print_exc()  # Prints the exact line number where it failed
        continue

    print(adata_merged.shape)
    remove_redundant_cols(adata_merged).write(output_path, compression='gzip')


#####---------------------------------------------------------------------
##Generate a matrix of cell counts per cell type across datasets.
#####---------------------------------------------------------------------

def generate_cell_count_matrix(df,
                            nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                           remaining_path="/analysis/adata_new_assigned_0828.h5ad",
                           back_up_remaining_path="/analysis/adata_assigned.h5ad",
                           filtering_based_on_barcodes = False,
                           scvi_labels_to_keep = None,
                           min_cells=50,
                           cell_type_col='new_cell_type',
                           cell_type_col_backup='assignments',
                           ID_col='SRA_ID',
                           include_new_samples = False

                             ):


    """
    Generate a matrix of cell counts per cell type across datasets.

    """

    

    all_cell_types = set()
    dataset_counts = {}
    total_counts = {}

    # First pass: collect all cell types and counts
    print("First pass: collecting cell types and counts...")
    for i in range(len(df)):
        ID = df[ID_col][i]
        print(ID)
        core_value = df[df[ID_col]==ID]["Core"].values[0]
        core_value = int(core_value)
        core_action = False
        if include_new_samples and core_value != 0 and core_value != 1:
          core_action = True

        try:
            print(f"\nProcessing {ID}...")
            try:
              adata = sc.read(f"{nac_path}{ID}{remaining_path}")
              if filtering_based_on_barcodes:
                adata = remove_low_quality_clusters(adata, scvi_labels_to_keep, ID, df, ignore_core=core_action)
            except:
              print("Falling back to adata_assigned")
              adata = sc.read(f"{nac_path}{ID}{back_up_remaining_path}")
              if filtering_based_on_barcodes:
                adata = remove_low_quality_clusters(adata, scvi_labels_to_keep, ID, df, ignore_core=core_action)

            try:
              cell_type_counts = adata.obs[cell_type_col].astype(str).value_counts()
            except:
              cell_type_counts = adata.obs[cell_type_col_backup].astype(str).value_counts()
            valid_types = cell_type_counts[cell_type_counts >= min_cells]

            # Store counts and update cell types set
            dataset_counts[ID] = valid_types
            all_cell_types.update(valid_types.index)
            total_counts[ID] = len(adata)

            print(f"Found {len(valid_types)} cell types with ≥{min_cells} cells")
            print(f"Total cells: {total_counts[ID]}")

        except Exception as e:
            print(f"Error processing {ID}: {str(e)}")
            continue

    # Sort cell types for consistent column order
    all_cell_types = sorted(list(all_cell_types))
    print(f"\nTotal unique cell types across all datasets: {len(all_cell_types)}")

    # Create matrix with datasets as rows and cell types as columns
    matrix_data = []
    for ID in dataset_counts.keys():
        row_data = {'dataset': ID}
        counts = dataset_counts[ID]

        # Fill in counts for each cell type
        for cell_type in all_cell_types:
            row_data[cell_type] = int(counts.get(cell_type, 0))

        matrix_data.append(row_data)

    # Convert to polars DataFrame
    cell_count_matrix = pl.DataFrame(matrix_data)

    return cell_count_matrix, total_counts


def save_count_matrix(matrix, total_counts, output_dir="/content/results"):
    """
    Save the cell count matrix and related information.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Save main count matrix
    matrix.write_csv(f"{output_dir}/cell_type_counts.csv")

    # Save total counts
    total_counts_df = pl.DataFrame({
        'dataset': list(total_counts.keys()),
        'total_cells': list(total_counts.values())
    })
    total_counts_df.write_csv(f"{output_dir}/total_cell_counts.csv")



def dataset_ct_matrix(atlas_source_path = '/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
                      output_dir="output_directory",
                      Modality = ['sc','sn','multi_rna'],
                      n_cells_threshold = 300,

                      nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                      remaining_path="/analysis/adata_new_assigned_0828.h5ad",
                      back_up_remaining_path="/analysis/adata_assigned.h5ad",
                      filtering_based_on_barcodes = False,
                      scvi_labels_to_keep = None,
                      min_cells=50,
                      cell_type_col='new_cell_type',
                      cell_type_col_backup='assignments',
                      ID_col='SRA_ID',
                      include_new_samples = False
                      ):
  #to display users at the start how many cell of each cell type we have
  import pandas as pd
  df = pd.read_excel(atlas_source_path)
  df = df[df['Modality'].isin(Modality)]
  df = df[df["n_cells"]>n_cells_threshold]
  df.reset_index(drop=True, inplace=True)
  #reset
  df = df.reset_index(drop=True)
  count_matrix, totals = generate_cell_count_matrix(
                              df,
                              nac_path=nac_path,
                              remaining_path=remaining_path,
                              back_up_remaining_path=back_up_remaining_path,
                              filtering_based_on_barcodes = filtering_based_on_barcodes,
                              scvi_labels_to_keep = scvi_labels_to_keep,
                              min_cells=min_cells,
                              cell_type_col=cell_type_col,
                              cell_type_col_backup=cell_type_col_backup,
                              ID_col=ID_col,
                              include_new_samples = include_new_samples
                              )
  save_count_matrix(count_matrix, totals, output_dir)

#---------------------------------------------------------------------
# Generate values for dotplot construction
#---------------------------------------------------------------------


def get_expression_metrics(adata,
                           min_abundance=50,
                           cell_type_col='new_cell_type',
                           cell_type_col_backup='assignments',
                            min_expr=0,
                            name=''):


    """Calculate expression metrics with genes as columns"""
    # Get cell types and filter by minimum cell count (50)
    try:
      cell_type_counts = adata.obs[cell_type_col].astype(str).value_counts()
    except:
      cell_type_counts = adata.obs[cell_type_col_backup].astype(str).value_counts()
    valid_cell_types = cell_type_counts[cell_type_counts >= min_abundance].index.values
    print(f"Cell types with >={min_abundance} cells: {len(valid_cell_types)} out of {len(cell_type_counts)}")
    genes = adata.var_names
    print(f"Number of genes: {len(genes)}")

    # Initialize result arrays
    n_cell_types = len(valid_cell_types)
    n_genes = len(genes)
    means = np.zeros((n_cell_types, n_genes))
    fracs = np.zeros((n_cell_types, n_genes))

    # Convert sparse matrix to dense once
    X = adata.X.toarray() if sp.issparse(adata.X) else np.array(adata.X)
    print(f"Matrix shape: {X.shape}")

    # Calculate metrics for each valid cell type
    for i, ct in enumerate(valid_cell_types):
        try:
          mask = adata.obs[cell_type_col].astype(str) == ct
          ct_expr = X[mask]
        except:
          mask = adata.obs[cell_type_col_backup].astype(str) == ct
          ct_expr = X[mask]
        print(f"  Cell type {ct}: {sum(mask)} cells")

        # Vectorized calculations
        means[i] = np.mean(ct_expr, axis=0)
        fracs[i] = np.mean(ct_expr > min_expr, axis=0)

    # Create DataFrames with SRA_ID prefixed cell types
    prefixed_cell_types = [f"{name}_{ct}" if name else ct for ct in valid_cell_types]

    # Create dictionary for DataFrame construction
    means_dict = {'cell_type': prefixed_cell_types}
    fracs_dict = {'cell_type': prefixed_cell_types}

    for j, gene in enumerate(genes):
        means_dict[gene] = means[:, j]  # Use all rows for this gene column
        fracs_dict[gene] = fracs[:, j]  # Use all rows for this gene column

    means_df = pl.DataFrame(means_dict)
    fracs_df = pl.DataFrame(fracs_dict)

    return means_df, fracs_df




def process_all_datasets(df,
                         nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                           remaining_path="/analysis/adata_new_assigned_0828.h5ad",
                           back_up_remaining_path="/analysis/adata_assigned.h5ad",
                           filtering_based_on_barcodes = False,
                           scvi_labels_to_keep = None,
                           ID_col='SRA_ID',
                           min_abundance=50,
                           cell_type_col='new_cell_type',
                           cell_type_col_backup='assignments',
                            min_expr=0,
                            include_new_samples = False
                           ):
    
    """Process datasets with genes as columns"""
    all_matrix1s = []
    all_matrix2s = []
    all_genes = set()

    # First pass: collect all genes
    print("First pass: collecting all genes...")
    for i in range(len(df)):
        ID = df[ID_col][i]
        print(f"\nProcessing {ID}...")
        try:
            try:
              adata = sc.read(f"{nac_path}{ID}{remaining_path}")
            except:
              print("Falling back to adata_assigned")
              adata = sc.read(f"{nac_path}{ID}{back_up_remaining_path}")
            all_genes.update(adata.var_names)
        except Exception as e:
            print(f"Error reading genes from {ID}: {str(e)}")
            continue

    all_genes = sorted(list(all_genes))
    print(f"Total unique genes across all datasets: {len(all_genes)}")

    # Second pass: process datasets and align genes
    print("\nSecond pass: processing datasets...")
    for i in range(len(df)):
        ID = df[ID_col][i]
        core_value = df[df[ID_col]==ID]["Core"].values[0]
        core_value = int(core_value)
        core_action = False
        if include_new_samples and core_value != 0 and core_value != 1:
          core_action = True
        try:
            print(f"\nProcessing {ID}...")
            try:
              adata = sc.read(f"{nac_path}{ID}{remaining_path}")
              if filtering_based_on_barcodes:
                adata = remove_low_quality_clusters(adata, scvi_labels_to_keep, ID, df, ignore_core=core_action)
              
            except:
              print("Falling back to adata_assigned")
              adata = sc.read(f"{nac_path}{ID}{back_up_remaining_path}")
              if filtering_based_on_barcodes:
                adata = remove_low_quality_clusters(adata, scvi_labels_to_keep, ID, df, ignore_core=core_action)


            if ("1" in adata.var_names) or ("10" in adata.var_names) or (10 in adata.var_names) or (1 in adata.var_names):
              #throw error and break
              raise Exception("1 or 10 in var_names")
              break

            #total norm and log1p
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

            # Get metrics with genes as columns
            matrix_1, matrix_2 = get_expression_metrics(adata, 
                                                        min_abundance=min_abundance,
                                                        cell_type_col= cell_type_col,
                                                        cell_type_col_backup=cell_type_col_backup,
                                                        min_expr=min_expr,
                                                        name=ID)
            print(f"Initial matrix shapes - M1: {matrix_1.shape}, M2: {matrix_2.shape}")

            # Add missing genes with zeros
            missing_genes = set(all_genes) - set(matrix_1.columns[1:])  # Skip cell_type column
            print(f"Adding {len(missing_genes)} missing genes")
            if missing_genes:
                new_cols = {gene: pl.lit(0.0) for gene in missing_genes}
                matrix_1 = matrix_1.with_columns(**new_cols)
                matrix_2 = matrix_2.with_columns(**new_cols)

            # Ensure columns are in the same order
            col_order = ['cell_type'] + sorted(all_genes)
            matrix_1 = matrix_1.select(col_order)
            matrix_2 = matrix_2.select(col_order)

            all_matrix1s.append(matrix_1)
            all_matrix2s.append(matrix_2)

            print(f"Completed {ID}, shape: {matrix_1.shape}")

        except Exception as e:
            print(f"Error processing {ID}: {str(e)}")
            continue

    # Combine matrices by concatenating rows
    if all_matrix1s and all_matrix2s:
        final_matrix1 = pl.concat(all_matrix1s, how="vertical")
        final_matrix2 = pl.concat(all_matrix2s, how="vertical")
        return final_matrix1, final_matrix2
    else:
        raise ValueError("No datasets were successfully processed")



def save_results(final_matrix1, final_matrix2, output_dir="results"):
    """Save results with genes as columns"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    for idx, (matrix, name) in enumerate([(final_matrix1, "matrix1"), (final_matrix2, "matrix2")]):
        print(f"Saving {name}")

        # Save cell types list
        cell_types = matrix.select("cell_type").to_series().to_list()
        with open(f"{output_dir}/{name}_rows.tsv", 'w') as f:
            f.write('\n'.join(cell_types))

        # Save gene names (columns)
        gene_cols = [col for col in matrix.columns if col != "cell_type"]
        with open(f"{output_dir}/{name}_genes.tsv", 'w') as f:
            f.write('\n'.join(gene_cols))

        # Save matrix data
        data_matrix = matrix.drop("cell_type").to_numpy()
        sparse_matrix = sp.csr_matrix(data_matrix)
        scipy.io.mmwrite(f"{output_dir}/{name}.mtx", sparse_matrix)


def dotplot_matrix(atlas_source_path = '/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
                  output_dir="results",
                  Modality = ['sc','sn','multi_rna'],
                  n_cells_threshold = 300,
                  nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                           remaining_path="/analysis/adata_new_assigned_0828.h5ad",
                           back_up_remaining_path="/analysis/adata_assigned.h5ad",
                           filtering_based_on_barcodes = False,
                           scvi_labels_to_keep = None,
                           ID_col='SRA_ID',
                           min_abundance=50,
                           cell_type_col='new_cell_type',
                           cell_type_col_backup='assignments',
                            min_expr=0,
                            include_new_samples = False
                  ):

  df = pd.read_excel(atlas_source_path)
  pituitary_atlas = df.copy()
  df = df[df['Modality'].isin(Modality)]
  df = df[df["n_cells"]>n_cells_threshold]
  df.reset_index(drop=True, inplace=True)

  final_matrix1, final_matrix2 = process_all_datasets(df,
                              nac_path=nac_path,
                           remaining_path=remaining_path,
                           back_up_remaining_path=back_up_remaining_path,
                           filtering_based_on_barcodes = filtering_based_on_barcodes,
                           scvi_labels_to_keep = scvi_labels_to_keep,
                           ID_col=ID_col,
                           min_abundance=min_abundance,
                           cell_type_col=cell_type_col,
                           cell_type_col_backup=cell_type_col_backup,
                            min_expr=min_expr,
                            include_new_samples = include_new_samples

                              )
  save_results(final_matrix1, final_matrix2, output_dir=output_dir)



#---------------------------------------------------------------------
# Generate values for cell type proportion analysis
#---------------------------------------------------------------------

def get_cell_type_abundances(adata,
                             cell_type_col='new_cell_type',
                              cell_type_col_backup='assignments'
                              ):
    """Calculate cell type abundances as percentages"""
    try:
      cell_types = adata.obs[cell_type_col].astype(str).value_counts()
    except:
      cell_types = adata.obs[cell_type_col_backup].astype(str).value_counts()
    total_cells = len(adata)

    # Calculate percentages
    percentages = (cell_types / total_cells) * 100
    return percentages.index.values, percentages.values

def process_datasets_abundances(df,
                              nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                            remaining_path="/analysis/adata_new_assigned_0828.h5ad",
                            back_up_remaining_path="/analysis/adata_assigned.h5ad",
                            scvi_labels_to_keep = None,
                            ID_col='SRA_ID',
                            cell_type_col='new_cell_type',
                            cell_type_col_backup='assignments',
                            filtering_based_on_barcodes = False,
                            include_new_samples = False
                            ):
    """Process datasets and create abundance matrix"""
    all_cell_types = set()

    # First pass: collect all cell types
    print("First pass: collecting all cell types...")
    for i in range(len(df)):
        ID = df[ID_col][i]
        try:
            try:
              adata = sc.read(f"{nac_path}{ID}{remaining_path}")
            except:
              print("Falling back to adata_assigned")
              adata = sc.read(f"{nac_path}{ID}{back_up_remaining_path}")
            cell_types, _ = get_cell_type_abundances(adata,
                                                      cell_type_col=cell_type_col,
                                                      cell_type_col_backup=cell_type_col_backup)
            all_cell_types.update(cell_types)
        except Exception as e:
            print(f"Error reading {ID}: {str(e)}")
            continue

    all_cell_types = sorted(list(all_cell_types))
    n_cell_types = len(all_cell_types)
    print(f"Total unique cell types: {n_cell_types}")

    # Create cell type to index mapping
    cell_type_to_idx = {ct: i for i, ct in enumerate(all_cell_types)}

    # Second pass: fill abundance matrix
    abundance_matrix = np.zeros((n_cell_types, len(df)))
    valid_samples = []

    print("\nSecond pass: calculating abundances...")
    for i in range(len(df)):
        ID = df[ID_col][i]
        core_value = df[df[ID_col]==ID]["Core"].values[0]
        core_value = int(core_value)
        core_action = False
        if include_new_samples and core_value != 0 and core_value != 1:
          core_action = True
        try:
            print(f"Processing {ID}...")
            try:
              adata = sc.read(f"{nac_path}{ID}{remaining_path}")
              if filtering_based_on_barcodes:
                adata = remove_low_quality_clusters(adata, scvi_labels_to_keep, ID, df, ignore_core=core_action)

            except:
              print("Falling back to adata_assigned")
              adata = sc.read(f"{nac_path}{ID}{back_up_remaining_path}")
              if filtering_based_on_barcodes:
                adata = remove_low_quality_clusters(adata, scvi_labels_to_keep, ID, df, ignore_core=core_action)

            cell_types, percentages = get_cell_type_abundances(adata,
                                                              cell_type_col=cell_type_col,
                                                              cell_type_col_backup=cell_type_col_backup)
            # Fill abundance matrix
            for ct, pct in zip(cell_types, percentages):
                abundance_matrix[cell_type_to_idx[ct], i] = pct

            valid_samples.append(ID)
            print(f"Completed {ID}")

        except Exception as e:
            print(f"Error processing {ID}: {str(e)}")
            continue

    return abundance_matrix, all_cell_types, valid_samples

def save_abundance_results(abundance_matrix, cell_types, samples, output_dir="abundance_results"):
    """Save abundance results in MTX format"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Save cell types (rows)
    with open(f"{output_dir}/abundance_rows.tsv", 'w') as f:
        f.write('\n'.join(cell_types))

    # Save sample IDs (columns)
    with open(f"{output_dir}/abundance_cols.tsv", 'w') as f:
        f.write('\n'.join(samples))

    # Save abundance matrix
    sparse_matrix = sp.csr_matrix(abundance_matrix)
    scipy.io.mmwrite(f"{output_dir}/abundance.mtx", sparse_matrix)


def cell_type_proportions(atlas_source_path = '/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
                         Modality = ['sc','sn','multi_rna'],
                         n_cells_threshold = 300,
                         nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                            remaining_path="/analysis/adata_new_assigned_0828.h5ad",
                            back_up_remaining_path="/analysis/adata_assigned.h5ad",
                            scvi_labels_to_keep = None,
                            ID_col='SRA_ID',
                            cell_type_col='new_cell_type',
                            cell_type_col_backup='assignments',
                            filtering_based_on_barcodes = False,
                          output_dir="abundance_results",
                          include_new_samples = False
                         ):

  df = pd.read_excel(atlas_source_path)
  df = df[df['Modality'].isin(Modality)]
  df = df[df["n_cells"]>n_cells_threshold]
  df.reset_index(drop=True, inplace=True)
  #reset
  df = df.reset_index(drop=True)
  # Usage
  matrix, cell_types, samples = process_datasets_abundances(df,
                              nac_path=nac_path,
                            remaining_path=remaining_path,
                            back_up_remaining_path=back_up_remaining_path,
                            scvi_labels_to_keep = scvi_labels_to_keep,
                            ID_col=ID_col,
                            cell_type_col=cell_type_col,
                            cell_type_col_backup=cell_type_col_backup,
                            filtering_based_on_barcodes = filtering_based_on_barcodes,
                            include_new_samples = include_new_samples

                              )
  save_abundance_results(matrix, cell_types, samples, output_dir=output_dir)

#---------------------------------------------------------------------
# Processing adatas for epitome h5 files
#---------------------------------------------------------------------

def simplify_and_process_adatas(df,
                                output_dir="/content/drive/MyDrive/pituitary_atlas/Unstructured/epitome_h5_files_rna",
                                nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                                remaining_path="/analysis/adata_new_assigned_0828.h5ad",
                                back_up_remaining_path="/analysis/adata_assigned.h5ad",
                                ID_col='SRA_ID',
                                calc_umaps=True,
                                cell_type_col='new_cell_type',
                                cell_type_col_backup='assignments'
                                ):
    
    """Simplify and process adatas, then save as h5ad files."""

    for i in range(len(df)):
        ID = df[ID_col][i]

        if Path(f"{output_dir}/{ID}_processed.h5ad").exists():
            print(f"Skipping {ID} as it has already been processed")
            continue
        try:
            print(f"\nProcessing {ID}...")
            try:
              adata = sc.read(f"{nac_path}/{ID}{remaining_path}")
            except:
              print("Falling back to adata_assigned")
              adata = sc.read(f"{nac_path}/{ID}{back_up_remaining_path}")

            # Simplify var and obs
            adata.var = adata.var.iloc[:, 0:0]  # Keep only index
            try:
              adata.obs = adata.obs[[cell_type_col]]
            except:
              try:
                adata.obs = adata.obs[[cell_type_col_backup]]
                #rename to "new_cell_type"
                adata.obs.rename(columns={cell_type_col_backup: cell_type_col}, inplace=True)

              except:
                adata.obs = adata.obs[["cell_type"]]
                #rename to "new_cell_type"
                adata.obs.rename(columns={"cell_type": cell_type_col}, inplace=True)

            #sc.pp.filter_genes(adata, min_cells=20)
            # Processing pipeline

            adata.obs.rename(columns={cell_type_col: "cell_type"}, inplace=True)

            adata.layers["counts"] = adata.X.copy()

            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            adata.layers["log1p"] = adata.X.copy()
            sc.pp.highly_variable_genes(adata, n_top_genes=2000)
            if calc_umaps:
              sc.pp.scale(adata, max_value=10)
              sc.tl.pca(adata, svd_solver='arpack')
              sc.pp.neighbors(adata, n_neighbors=25, n_pcs=42)
              sc.tl.leiden(adata, flavor='igraph')
              sc.tl.umap(adata)
            else:
              print("Skipping UMAP calculation")

            adata.X = adata.layers["log1p"].copy()

            sc.pl.umap(adata, color=['cell_type'])

            #make sure X is a sparse matrix
            if not sp.issparse(adata.X):
                adata.X = sp.csr_matrix(adata.X)

            #repeat for all layers
            for layer in adata.layers:
                if not sp.issparse(adata.layers[layer]):
                    adata.layers[layer] = sp.csr_matrix(adata.layers[layer])
            # Save processed adata
            output_path = f"{output_dir}/{ID}_processed.h5ad"
            adata.write_h5ad(output_path, compression='gzip')
            print(f"Saved to {output_path}")

        except Exception as e:
            print(f"Error processing {ID}: {str(e)}")
            continue


#wrap function to call save_all_h5s
def save_all_h5s(atlas_source_path = '/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
                  Modality = ['sc','sn','multi_rna'],
                  n_cells_threshold = 300,
                  output_dir="/content/drive/MyDrive/pituitary_atlas/Unstructured/epitome_h5_files_rna",
                                nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                                remaining_path="/analysis/adata_new_assigned_0828.h5ad",
                                back_up_remaining_path="/analysis/adata_assigned.h5ad",
                  ID_col='SRA_ID',
                  calc_umaps=True,
                  cell_type_col='new_cell_type',
                  cell_type_col_backup='assignments',

                  ):
   
  df = pd.read_excel(atlas_source_path)
  pituitary_atlas = df.copy()
  df = df[df['Modality'].isin(Modality)]
  df = df[df["n_cells"]>n_cells_threshold]
  df.reset_index(drop=True, inplace=True)
  #reset
  df = df.reset_index(drop=True)
  df

  
  # Usage
  simplify_and_process_adatas(df,
                              output_dir=output_dir,
                              nac_path=nac_path,
                              remaining_path=remaining_path,
                              back_up_remaining_path=back_up_remaining_path,
                              ID_col=ID_col,
                              calc_umaps=calc_umaps,
                              cell_type_col=cell_type_col,
                              cell_type_col_backup=cell_type_col_backup

                              )


#function to bundle up all h5 files into zip files of 50 each
def bundle_up_all_h5s(source_folder='/content/drive/MyDrive/pituitary_atlas/Unstructured/epitome_h5_files_rna/',
                      base_output_zip='/content/drive/MyDrive/pituitary_atlas/Unstructured/epitome_h5_files_rna/epitome_grouped_h5_files_0828'):
  """Bundles h5 files into zip files of specified chunk size."""
   
  # Get the list of all files in the folder
  all_files = [f for f in os.listdir(source_folder) if os.path.isfile(os.path.join(source_folder, f))]
  len(all_files)

  # Group files into chunks of 50
  chunk_size = 15
  file_chunks = [all_files[i:i + chunk_size] for i in range(0, len(all_files), chunk_size)]

  # Create a zip file for each chunk of files
  for idx, chunk in enumerate(file_chunks):
      zip_filename = f'{base_output_zip}_group_{idx + 1}.zip'
      zip_filepath = os.path.join(source_folder, zip_filename)

      with zipfile.ZipFile(zip_filepath, 'w', zipfile.ZIP_DEFLATED) as zipf:
          for file in chunk:
              file_path = os.path.join(source_folder, file)
              zipf.write(file_path, os.path.basename(file_path))  # Write file to the zip archive

      print(f'Created zip file: {zip_filename}')



def zip_summary_pdfs(atlas_source_path = '/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
                     Modality = ['sc','sn','multi_rna'],
                     n_cells_threshold = 300,
                     output_dir="/content/drive/MyDrive/pituitary_atlas/Unstructured/epitome_pdf_files",
                     output_zip_name="summary_pdfs.zip",
                     pdf_path_prefix="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                     pdf_path_suffix="/analysis/summary.pdf"
                     ):
    """Zips summary PDFs for all SRA IDs in the DataFrame, using SRA_ID in filename."""
    df = pd.read_excel(atlas_source_path)
    pituitary_atlas = df.copy()
    df = df[df['Modality'].isin(Modality)]
    df = df[df["n_cells"]>n_cells_threshold]
    df.reset_index(drop=True, inplace=True)

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_zip_path = os.path.join(output_dir, output_zip_name)

    with zipfile.ZipFile(output_zip_path, 'w') as zipf:
        for sra_id in df["SRA_ID"]:
            pdf_path = f"{pdf_path_prefix}{sra_id}{pdf_path_suffix}"
            if os.path.exists(pdf_path):
                # Create new filename with SRA_ID
                new_filename = f"summary_{sra_id}.pdf"
                with open(pdf_path, "rb") as f:
                  zip_info = ZipInfo(new_filename)
                  zipf.writestr(zip_info, f.read())
                print(f"Added {pdf_path} as {new_filename} to zip")
            else:
                print(f"PDF not found: {pdf_path}")

    print(f"Zip file created at: {output_zip_path}")



import os
import re
import gc
import zipfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import scanpy as sc


# --- 1. Worker (must stay at module level for ProcessPoolExecutor) ---
def save_gene_batch_parquet(gene_list, X_batch, obs_names, out_dir="genes_parquet"):
    """
    Worker: receives ONLY the columns for its genes (already sliced),
    so the full matrix is never re-sent to each process.
    """
    if hasattr(X_batch, "toarray"):
        X_batch = X_batch.toarray()

    batch_df = pd.DataFrame(X_batch, index=obs_names, columns=gene_list)

    file_paths = []
    for gene in gene_list:
        safe_name = re.sub(r"[^\w.\-]", "_", gene)  # avoid '/' etc. in filenames
        file_path = os.path.join(out_dir, f"{safe_name}.parquet")
        batch_df[[gene]].to_parquet(file_path, compression="zstd")
        file_paths.append(file_path)

    del batch_df, X_batch
    gc.collect()
    return file_paths


def _make_cell_ids(d):
    """cell_id = SRA_ID + '_' + barcode (from obs['cell_barcode'], obs['0'], obs[0], or obs_names)."""
    if "cell_barcode" in d.obs.columns:
        ids = d.obs["cell_barcode"].astype(str)
    elif "0" in d.obs.columns:
        ids = d.obs["0"].astype(str)
    elif 0 in d.obs.columns:
        ids = d.obs[0].astype(str)
    else:
        ids = pd.Series(d.obs_names.astype(str), index=d.obs_names)
    d.obs["cell_id"] = d.obs["SRA_ID"].astype(str) + "_" + ids


# --- 2. Main function ---
def generate_atlas_umap(
    atlas_path="/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx",
    dataset_with_coords="/content/drive/MyDrive/pituitary_atlas/Unstructured/scanvi_all_cells_0828.h5ad",
    wt_dataset="/content/drive/MyDrive/pituitary_atlas/Unstructured/wt_adata_merged_post_scvi_0902.h5ad",
    mut_dataset=None,  # optional: None or "" -> WT only
    output_total_counts_csv="/content/drive/MyDrive/pituitary_atlas/Unstructured/total_counts.csv",
    zip_output_path="/content/drive/MyDrive/pituitary_atlas/Unstructured/adata_export_large_umap_1015.zip",
    batch_size=100,
    max_workers=2,  # keep low in Colab to avoid OOM
    atlas_meta_cols=("Comp_sex", "10X version", "Modality", "Normal"),  # taken from the atlas sheet if present
):
    # ---- Load ----
    adata = sc.read_h5ad(dataset_with_coords)

    paths = {k: p for k, p in {"wt": wt_dataset, "mut": mut_dataset}.items() if p}
    if not paths:
        raise ValueError("Provide at least one expression dataset (wt_dataset and/or mut_dataset).")
    for k, p in paths.items():
        if not os.path.exists(p):
            raise FileNotFoundError(f"{k} dataset not found: {p}")
    expr_adatas = {k: sc.read_h5ad(p) for k, p in paths.items()}
    for d in expr_adatas.values():
        d.var_names_make_unique()  # duplicate gene names break concat

    # ---- 1. Cell IDs ----
    for d in [adata, *expr_adatas.values()]:
        _make_cell_ids(d)

    # Parse samples: barcode column is "nan" -> use bc + subpool instead
    if "0" in adata.obs.columns and {"bc", "subpool"}.issubset(adata.obs.columns):
        nans = adata.obs["0"].astype(str) == "nan"
        adata.obs.loc[nans, "cell_id"] = (
            adata.obs.loc[nans, "SRA_ID"].astype(str) + "_"
            + adata.obs.loc[nans, "bc"].astype(str) + "_"
            + adata.obs.loc[nans, "subpool"].astype(str)
        )

    # ---- 2. Merge (or single dataset) ----
    if len(expr_adatas) == 1:
        (label, merged_adata), = expr_adatas.items()
        merged_adata.obs["batch"] = label
    else:
        merged_adata = sc.concat(list(expr_adatas.values()), label="batch", keys=list(expr_adatas.keys()), index_unique="-")
    merged_adata.obs_names_make_unique()
    merged_adata.var_names_make_unique()
    del expr_adatas
    gc.collect()

    # ---- 3. Map UMAP coords ----
    if not {"UMAP1", "UMAP2"}.issubset(adata.obs.columns):
        if "X_umap" not in adata.obsm:
            raise KeyError("Reference has neither obs['UMAP1'/'UMAP2'] nor obsm['X_umap'].")
        umap = np.asarray(adata.obsm["X_umap"])
        adata.obs["UMAP1"] = umap[:, 0]
        adata.obs["UMAP2"] = umap[:, 1]
        print("Using UMAP coordinates from obsm['X_umap'].")
    dup = adata.obs["cell_id"].duplicated(keep=False)
    if dup.any():
        print(f"WARNING: {dup.sum()} reference cells share a cell_id "
              f"(e.g. {adata.obs.loc[dup, 'cell_id'].iloc[0]!r}); they are ambiguous and will not be matched.")
    coords = adata.obs.loc[~dup].set_index("cell_id")[["UMAP1", "UMAP2"]]
    del adata
    gc.collect()
    merged_adata.obs["UMAP1"] = merged_adata.obs["cell_id"].map(coords["UMAP1"])
    merged_adata.obs["UMAP2"] = merged_adata.obs["cell_id"].map(coords["UMAP2"])
    n_before = merged_adata.n_obs
    merged_adata = merged_adata[~merged_adata.obs["UMAP1"].isna()].copy()
    print(f"Matched {merged_adata.n_obs}/{n_before} cells to UMAP coordinates.")
    if merged_adata.n_obs == 0:
        raise ValueError("No cells matched the UMAP file - check cell_id construction.")

    # ---- 4. Total counts ----
    merged_adata.obs["ncounts"] = np.asarray(merged_adata.X.sum(axis=1)).ravel()
    merged_adata.obs["ncounts"].to_csv(output_total_counts_csv)

    # ---- 5. Atlas metadata (keeps original cell names) ----
    atlas = pd.read_excel(atlas_path).drop_duplicates(subset=["SRA_ID"])
    meta_cols = [c for c in atlas_meta_cols if c in atlas.columns]
    missing = [c for c in atlas_meta_cols if c not in atlas.columns]
    if missing:
        print(f"Atlas sheet has no {missing}; keeping any existing obs values for those.")
    atlas = atlas[meta_cols + ["SRA_ID"]]
    atlas["SRA_ID"] = atlas["SRA_ID"].astype(str)

    # Only drop obs columns that the atlas will replace (plus legacy 'Sex')
    to_drop = [c for c in meta_cols + ["Sex"] if c in merged_adata.obs.columns]
    obs = merged_adata.obs.drop(columns=to_drop)
    obs["SRA_ID"] = obs["SRA_ID"].astype(str)
    new_obs = obs.merge(atlas, on="SRA_ID", how="left")
    new_obs.index = obs.index  # merge resets the index; restore cell names
    merged_adata.obs = new_obs

    # ---- 6. Prepare for export ----
    if hasattr(merged_adata.X, "format") and merged_adata.X.format != "csc":
        merged_adata.X = merged_adata.X.tocsc()  # fast column slicing

    obs_path, umap_path, gene_dir = "obs.parquet", "umap.parquet", "genes_parquet"
    merged_adata.obs.to_parquet(obs_path, compression="gzip")
    merged_adata.obs[["UMAP1", "UMAP2"]].to_parquet(umap_path, compression="gzip")

    # ---- 7. Parallel gene export ----
    os.makedirs(gene_dir, exist_ok=True)
    X = merged_adata.X
    obs_names = merged_adata.obs_names
    gene_names = merged_adata.var_names.to_list()
    n_genes = len(gene_names)

    print(f"Starting parallel export of {n_genes} genes...")
    gene_files = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                save_gene_batch_parquet,
                gene_names[i:i + batch_size],
                X[:, i:i + batch_size],  # only this batch's columns are sent
                obs_names,
                gene_dir,
            )
            for i in range(0, n_genes, batch_size)
        ]
        for f in futures:
            gene_files.extend(f.result())

    # ---- 8. Zip + cleanup ----
    with zipfile.ZipFile(zip_output_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(obs_path)
        zf.write(umap_path)
        for fp in gene_files:
            zf.write(fp)

    os.remove(obs_path)
    os.remove(umap_path)
    for fp in gene_files:
        os.remove(fp)
    os.rmdir(gene_dir)

    print(f"Success! Saved to {zip_output_path}")
    return merged_adata




##########---------------------------------------------------------------------
#ATAC
##########---------------------------------------------------------------------

def processing_atac_samples(atac_all_cells,
                            atlas_source_path = '/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
                            target_folder='/content/drive/MyDrive/pituitary_atlas/Unstructured/epitome_h5_files_atac_0911_final/',
                            target_folder_pdfs="/content/drive/MyDrive/pituitary_atlas/Unstructured/epitome_pdf_files_atac_0911_final/",
                            output_zip_name_path = "/content/drive/MyDrive/pituitary_atlas/Unstructured/epitome_h5_files_atac_0911_final.zip",
                            output_zip_name_pdfs_path = "/content/drive/MyDrive/pituitary_atlas/Unstructured/epitome_pdf_files_atac_0911_final.zip"

):

  pituitary_atlas = pd.read_excel(atlas_source_path)
  #make this folder if doesnt exist
  if not os.path.exists(target_folder):
      os.makedirs(target_folder)

  if not os.path.exists(target_folder_pdfs):
      os.makedirs(target_folder_pdfs)

  for geo in df["GEO"].values:

    print(geo)
    #try:
    adata = atac_all_cells[atac_all_cells.obs["sample"]==geo].copy()
    adata.obs["GEO"] = adata.obs["sample"]

    n_cells = adata.shape[0]
    print(adata.shape)
    #find this GEO value in pituitary_atlas, and add this to the col n_cells
    pituitary_atlas.loc[pituitary_atlas["GEO"] == geo, "n_cells"] = n_cells
    pituitary_atlas.loc[pituitary_atlas["GEO"] == geo, "most_recent_workflow"] = "v_0.01"

    adata.obs = adata.obs.merge(pituitary_atlas, on="GEO").copy()

    print(adata.shape)
    #write adata
    # Convert all columns in obs to strings if needed
    for col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype(str)

    #change index to str
    adata.obsm["_scvi_extra_categorical_covs"].index = ["barcode_" + str(i) for i in adata.obs.index.astype(int)]
    adata.obs.index = ["barcode_" + str(i) for i in adata.obs.index.astype(int)]

    adata = adata.copy()
    #write adata
    adata.write(target_folder+geo+".h5ad")

    # save atac PDFs
    try:
      pdf_path = "/content/drive/MyDrive/pituitary_atlas/Unstructured/from_pc_0910/"+geo+"_summary.pdf"
      shutil.copy(pdf_path, target_folder_pdfs+geo+".pdf")
    except:
      print("No pdf found for "+geo)

    #except:
      print("Error reading "+geo)
      pituitary_atlas.loc[pituitary_atlas["GEO"] == geo, "n_cells"] = 0
      pituitary_atlas.loc[pituitary_atlas["GEO"] == geo, "most_recent_workflow"] = "none"
      continue

  #resave back to atlas
  pituitary_atlas.to_excel('/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx', index=False)


  
  # Zip the folder
  shutil.make_archive(output_zip_name_path.replace('.zip', ''), 'zip', target_folder)
  shutil.make_archive(output_zip_name_pdfs_path.replace('.zip', ''), 'zip', target_folder_pdfs)
  print(f"Zipped h5 files to {output_zip_name_path}")
  print(f"Zipped pdf files to {output_zip_name_pdfs_path}")


##############
#### Pseudobulking RNA
##############
##############
#### Pseudobulking RNA
##############

def _safe_int(x, default=None):
  """int() that tolerates NaN / blanks / strings like '1.0'."""
  try:
    v = pd.to_numeric(x, errors='coerce')
    return default if pd.isna(v) else int(v)
  except Exception:
    return default


def _read_first_nonempty_h5ad(nac_path, ID, suffixes):
  """Try each suffix in order; return the first readable AnnData with >0 cells, else None."""
  for p in suffixes:
    if not p:
      continue
    path = f"{nac_path}{ID}{p}"
    if not os.path.isfile(path):
      continue
    try:
      ad = sc.read_h5ad(path)
    except Exception as e:
      print(f"Could not read {path}: {type(e).__name__} - {e}")
      continue
    if ad.n_obs > 0:
      return ad
    print(f"{path} has 0 cells, trying next")
  return None


def generate_pseudobulks(atlas_source_path = '/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
                         modality = ['sc','sn','multi_rna'],
                         min_cells_threshold = 300,
                         scvi_labels_to_keep = None,
                         filtering_based_on_barcodes = True,
                         ID_col = 'SRA_ID',
                         nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
                         remaining_path="/analysis/adata_new_assigned_0828.h5ad",
                         remaining_path_backup="/analysis/adata_assigned.h5ad",
                         sample_col='sample',
                         groups_col='new_cell_type',
                         groups_col_backups=('assignments', 'cell_type'),
                         min_cells=50,
                         pdata_output_path_suffix="/analysis/pdata_assigned_0828.h5ad",
                         include_new_samples = False,
                         plot_hist = True
                         ):

  if filtering_based_on_barcodes and scvi_labels_to_keep is None:
    raise ValueError("filtering_based_on_barcodes=True needs scvi_labels_to_keep (a DataFrame with 'SRA_ID' and 'barcode').")

  df = pd.read_excel(atlas_source_path)
  df = df[df['Modality'].isin(modality)].reset_index(drop=True)

  missing = []

  for i in range(len(df)-1,-1,-1):
    row = df.iloc[i]
    ID = row[ID_col]
    if pd.isna(ID):
      continue
    ID = str(ID)
    SRA_ID = str(row['SRA_ID'])
    GEO = str(row['GEO'])
    author = row['Author']
    normal = _safe_int(row['Normal'])
    n_cells = pd.to_numeric(row['n_cells'], errors='coerce')
    core_value = _safe_int(row['Core'], default=-1)
    is_core = core_value in (0, 1)
    core_action = include_new_samples and not is_core  # True -> keep all cells, skip barcode filter

    print(ID)

    #if normal == 2 then skip
    if normal == 2:
      continue
    #if n_cells < threshold skip (unknown n_cells is attempted)
    if pd.notna(n_cells) and n_cells < min_cells_threshold:
      continue

    #skip if pdata already exists
    out_path = f"{nac_path}{ID}{pdata_output_path_suffix}"
    if os.path.isfile(out_path):
      print(f"Skipping {ID} as pdata already exists")
      continue

    #skip non-core samples unless asked to include them
    if not include_new_samples and not is_core:
      print(f"Skipping non-core sample {ID}")
      continue

    try:
      #read adata assigned (primary path, then backup; skip missing, unreadable or empty files)
      adata_assigned = _read_first_nonempty_h5ad(nac_path, ID, (remaining_path, remaining_path_backup))
      if adata_assigned is None:
        print(f" Missing/empty entry for {ID}")
        missing.append(ID)
        continue

      if plot_hist:
        cell_counts = np.asarray(adata_assigned.X.sum(axis=1)).ravel()
        plt.hist(cell_counts, color='red')
        plt.title(ID)
        plt.show()

      if filtering_based_on_barcodes and not core_action:
        obs = adata_assigned.obs
        cell_ids = obs['0'].astype(str) if '0' in obs.columns else (obs[0].astype(str) if 0 in obs.columns else pd.Series(adata_assigned.obs_names.astype(str), index=obs.index))

        if author == "Rebboah et al. (2025)":
          obs["barcode"] = GEO + ":" + obs["bc"].astype(str) + "_" + obs["subpool"].astype(str) + "-1"
        elif author == "Sochodolsky et al. (2026)":
          obs["barcode"] = SRA_ID + ":" + cell_ids + "-1"
        else:
          obs["barcode"] = GEO + ":" + cell_ids + "-1"

        keep = scvi_labels_to_keep.loc[scvi_labels_to_keep["SRA_ID"].astype(str) == SRA_ID, "barcode"].values
        n_before = adata_assigned.n_obs
        adata_assigned = adata_assigned[adata_assigned.obs["barcode"].isin(keep)].copy()
        print(f"Removed {n_before - adata_assigned.n_obs} cells")

        if adata_assigned.n_obs == 0:
          print(f"No cells left after barcode filtering for {ID}")
          missing.append(ID)
          continue

      # 1. Pick the cell type column (fall back if groups_col is missing)
      ct_col = next((c for c in (groups_col, *groups_col_backups) if c in adata_assigned.obs.columns), None)
      if ct_col is None:
        print(f"No cell type column ({groups_col} or {groups_col_backups}) for {ID}")
        missing.append(ID)
        continue
      if ct_col != groups_col:
        print(f"Using '{ct_col}' as cell type column for {ID}")
      cell_types = adata_assigned.obs[ct_col].astype(str)

      # 2. Restrict obs ONLY to the grouping columns to avoid metadata conflicts inside decoupler
      adata_assigned.obs = pd.DataFrame({sample_col: ID, groups_col: cell_types.values},
                                        index=adata_assigned.obs_names)

      # 3. Make observation/var index unique to prevent broadcasting crashes
      adata_assigned.obs_names_make_unique()
      adata_assigned.var_names_make_unique()

      # 4. Pseudobulk
      pdata_assigned = dc.pp.pseudobulk(
          adata_assigned,
          sample_col=sample_col,
          groups_col=groups_col,
          mode='sum'
      )

      # 5. Drop groups with too few cells (column name differs between decoupler versions)
      n_cells_col = "psbulk_cells" if "psbulk_cells" in pdata_assigned.obs.columns else "psbulk_n_cells"
      pdata_assigned = pdata_assigned[pd.to_numeric(pdata_assigned.obs[n_cells_col], errors='coerce') >= min_cells].copy()
      # store under both names so downstream code works with either decoupler version
      pdata_assigned.obs["psbulk_n_cells"] = pdata_assigned.obs[n_cells_col]

      if pdata_assigned.n_obs == 0:
        print(f"No cell types with >= {min_cells} cells for {ID}")
        missing.append(ID)
        continue

      # 6. Make obs/var writeable and save
      pdata_assigned.obs.columns = pdata_assigned.obs.columns.astype(str)
      pdata_assigned.var.columns = pdata_assigned.var.columns.astype(str)
      for col in pdata_assigned.obs.columns:
        pdata_assigned.obs[col] = pdata_assigned.obs[col].astype(str)

      os.makedirs(os.path.dirname(out_path), exist_ok=True)
      pdata_assigned.write(out_path)
      print(f"ALL DONE for {ID}")

    except Exception as e:
      import traceback
      print(f"Error processing {ID}: {type(e).__name__} - {e}")
      traceback.print_exc()
      missing.append(ID)
      continue

  print("Missing/skipped samples (no file, empty, no cell type column, too few cells, or error):")
  print(missing)
  return missing


def read_merge_save_pbulk(
      atlas_source_path = '/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
      output_path="/content/drive/MyDrive/pituitary_atlas/Unstructured/pdatas0828.h5ad",
      Modality = ['sc','sn','multi_rna'],
      total_n_cells_threshold = 200,
      nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
      remaining_path="/analysis/pdata_assigned_0828.h5ad",
      cell_type_col='new_cell_type',
      backup_cell_type_col='assignments',
      ID_col='SRA_ID',
      sample_col='sample'
      ):

  # whole-sample pseudobulks have no cell type: pass cell_type_col=None
  # (cell_type_col equal to ID_col or sample_col is treated the same way)
  if cell_type_col in (None, "", ID_col, sample_col):
    cell_type_col = None
  if backup_cell_type_col in (None, "", ID_col, sample_col):
    backup_cell_type_col = None
  ct_cols = [cell_type_col] if cell_type_col else []

  df = pd.read_excel(atlas_source_path)
  df = df[df['Modality'].isin(Modality)].reset_index(drop=True)
  df[ID_col] = df[ID_col].astype(str)

  pdatas = []
  failed = []
  for i in range(len(df)):
    row = df.iloc[i]
    #if normal == 2 then skip
    if _safe_int(row['Normal']) == 2:
      continue
    n_cells = pd.to_numeric(row['n_cells'], errors='coerce')
    if pd.notna(n_cells) and n_cells < total_n_cells_threshold:
      continue
    ID = row[ID_col]
    print(ID)

    path = f"{nac_path}{ID}{remaining_path}"
    if not os.path.isfile(path):
      print(f"No pdata for {ID}")
      failed.append(ID)
      continue
    try:
      p = sc.read_h5ad(path)
    except Exception as e:
      print(f"Error reading {ID}: {type(e).__name__} - {e}")
      failed.append(ID)
      continue
    if p.n_obs == 0:
      print(f"Empty pdata for {ID}")
      failed.append(ID)
      continue

    obs = p.obs
    # harmonise cell count column across decoupler versions
    if "psbulk_n_cells" not in obs.columns:
      if "psbulk_cells" in obs.columns:
        obs["psbulk_n_cells"] = obs["psbulk_cells"]
      else:
        obs["psbulk_n_cells"] = np.nan
    # harmonise cell type column
    if cell_type_col and cell_type_col not in obs.columns:
      if backup_cell_type_col and backup_cell_type_col in obs.columns:
        obs[cell_type_col] = obs[backup_cell_type_col]
      else:
        print(f"No {cell_type_col}/{backup_cell_type_col} column in pdata for {ID}")
        failed.append(ID)
        continue
    # harmonise sample column
    if sample_col not in obs.columns:
      obs[sample_col] = ID
    obs[ID_col] = ID

    p.obs = obs[list(dict.fromkeys([sample_col, *ct_cols, "psbulk_n_cells", ID_col]))].copy()
    p.var_names_make_unique()
    pdatas.append(p)

  print(f"Loaded {len(pdatas)} pdatas")
  if not pdatas:
    raise ValueError("No pseudobulk files could be loaded - check nac_path / remaining_path.")

  pdatas_merged = anndata.concat(pdatas, join='outer', fill_value=0)
  pdatas_merged.obs_names_make_unique()

  # change nans to 0 in X (dense or sparse)
  X = pdatas_merged.X
  if sp.issparse(X):
    X.data[np.isnan(X.data)] = 0
  else:
    pdatas_merged.X = np.nan_to_num(np.asarray(X, dtype=float), nan=0.0)

  pdatas_merged.obs["psbulk_n_cells"] = pd.to_numeric(pdatas_merged.obs["psbulk_n_cells"], errors='coerce')
  print(f"Total cells pseudobulked: {pdatas_merged.obs['psbulk_n_cells'].sum():.0f}")

  # fill missing / "nan" cell types and samples
  obs = pdatas_merged.obs.copy()
  for col in (sample_col, *ct_cols):
    obs[col] = obs[col].astype(object)
    bad = obs[col].isna() | obs[col].astype(str).isin(["nan", "None", ""])
    if col == sample_col:
      obs.loc[bad, col] = obs.loc[bad, ID_col]
    elif bad.any():
      print(f"WARNING: {bad.sum()} pseudobulks have no {cell_type_col}")
    obs[col] = obs[col].astype(str)

  # add metadata from the atlas sheet (keep cell names; one row per sample)
  meta = df.drop_duplicates(subset=[ID_col])
  meta = meta.drop(columns=[c for c in (sample_col, *ct_cols, "psbulk_n_cells") if c in meta.columns and c != ID_col])
  new_obs = obs[[sample_col, *ct_cols, "psbulk_n_cells"]].merge(
      meta, left_on=sample_col, right_on=ID_col, how='left')
  new_obs.index = obs.index
  pdatas_merged.obs = new_obs

  cols_to_remove = ['published', 'species',
       'Sample_wrong_sex', 'Sample_missing_sex', 'Sample_wrong_barcode_kit',
       'Sample_switched_up', 'Sample_kept_for_analysis',
       'Study_Seurat_version', 'Study_CR_version', 'Study_minimum_n',
       'Study_minimum_group', 'pseudoaligned', 'passed_qc', 'filtering_junk',
       'median_cellassign_prob', 'most_recent_workflow', 'n_cells',
       'processed_atac', 'passed_qc_tcc', 'tcc_pseudobulked', 'Core',
       'Notes 1', 'Notes 2']
  pdatas_merged.obs = pdatas_merged.obs.drop(columns=cols_to_remove, errors='ignore')

  # remove columns that contain "Unnamed" or _index
  cols = pdatas_merged.obs.columns.astype(str)
  pdatas_merged.obs = pdatas_merged.obs.loc[:, ~(cols.str.contains('Unnamed') | cols.str.contains('_index'))]

  pdatas_merged.obs.columns = pdatas_merged.obs.columns.astype(str)
  pdatas_merged.var.columns = pdatas_merged.var.columns.astype(str)
  pdatas_merged.var.index = pdatas_merged.var.index.astype(str)
  pdatas_merged.obs.index = pdatas_merged.obs.index.astype(str)

  # object columns (mixed types / NaN) -> str so h5ad writing doesn't fail
  for col in pdatas_merged.obs.columns:
    if pdatas_merged.obs[col].dtype == object:
      pdatas_merged.obs[col] = pdatas_merged.obs[col].astype(str)
  for col in pdatas_merged.var.columns:
    if pdatas_merged.var[col].dtype == object:
      pdatas_merged.var[col] = pdatas_merged.var[col].astype(str)

  os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
  pdatas_merged.write(output_path)
  print(f"Saved {pdatas_merged.shape} to {output_path}")
  if failed:
    print("Samples without a usable pdata file:")
    print(failed)
  return pdatas_merged