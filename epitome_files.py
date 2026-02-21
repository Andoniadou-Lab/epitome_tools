
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

from concurrent.futures import ProcessPoolExecutor
anndata.settings.allow_write_nullable_strings = True


#####---------------------------------------------------------------------
##Creating large merged object from all pituitary datasets
#####---------------------------------------------------------------------

def remove_redundant_cols(adata):
  # Use a lambda function to check if column name contains '0' but is not exactly '0'
  filter_condition = lambda x: '0' in str(x) and str(x) != '0'

  # Apply the filter to both adata.obs and adata.var
  adata.obs = adata.obs.loc[:, ~adata.obs.columns.map(filter_condition)]
  adata.var = adata.var.loc[:, ~adata.var.columns.map(filter_condition)]
  adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.contains('-1')]
  adata.var = adata.var.loc[:, ~adata.var.columns.str.contains('-1')]

  #rename obs from 1-len
  adata.obs.index = [f"{i+1}" for i in range(len(adata.obs))]

  #turn to string so saving won't throw error
  adata.obs.columns = adata.obs.columns.astype(str)
  adata.var.columns = adata.var.columns.astype(str)

  for col in adata.obs.columns:
    adata.obs[col] = adata.obs[col].astype(str)

  #remove obsm and layers to make object more lightweight

  adata.obsm = {}
  adata.layers = {}

  return adata


def remove_low_quality_clusters(adata, scvi_labels_to_keep, SRA_ID, df, ignore_core=False):

  author = df[df["SRA_ID"]==SRA_ID]["Author"].values[0]
  GEO = df[df["SRA_ID"]==SRA_ID]["GEO"].values[0]

  if ignore_core: #if not core sample, therefore not in scvi labels to keep, still keep the data by skipping filtering
    print(f"Ignoring core filtering for {SRA_ID}")
    return adata
  
  elif author == "Rebboah et al. (2025)":

    adata.obs["barcode"] = [GEO + ":" + adata.obs["bc"][i] + "_" + adata.obs["subpool"][i] + "-1"  for i in range(len(adata.obs))]
    shape1 = adata.shape
    adata = adata[adata.obs["barcode"].isin(scvi_labels_to_keep[scvi_labels_to_keep["SRA_ID"]==SRA_ID]["barcode"].values)]
    shape2 = adata.shape
    print(f"Removed {shape1[0]-shape2[0]} cells")

  elif author == "Sochodolsky et al. (2026)":
     
    adata.obs["barcode"] = [SRA_ID+":" + adata.obs['0'][i] + "-1" for i in range(len(adata.obs))]
    shape1 = adata.shape
    adata = adata[adata.obs["barcode"].isin(scvi_labels_to_keep[scvi_labels_to_keep["SRA_ID"]==SRA_ID]["barcode"].values)]
    shape2 = adata.shape
    print(f"Removed {shape1[0]-shape2[0]} cells")

  else:
    adata.obs["barcode"] = [GEO+":" + adata.obs['0'][i] + "-1" for i in range(len(adata.obs))]
    shape1 = adata.shape
    adata = adata[adata.obs["barcode"].isin(scvi_labels_to_keep[scvi_labels_to_keep["SRA_ID"]==SRA_ID]["barcode"].values)]
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
                           include_new_samples = False #samples that are Core == 2
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
          adata.obs[ID_col] = ID
          print("Added SRA_ID to obs")
          adata.obs["GEO"] = df[df[ID_col]==ID]["GEO"].values[0]
          print("Added GEO to obs")
          adata.obs["assay"] = df[df[ID_col]==ID]["Modality"].values[0]
          print("Added Modality to obs")

          adata.obs["orig_index"] = adata.obs.index.astype(str)
          core_action = False
          if include_new_samples and core_value == 2:
            core_action = True #true, therefore ignore filter and keep all cells
             
          if filtering_based_on_barcodes:
            adata = remove_low_quality_clusters(adata, scvi_labels_to_keep, ID, df,ignore_core=core_action)
          
          if not include_new_samples and core_value == 2:
            print(f"Not including non-core samples, removing {ID}")
            continue
            

          if ("1" in adata.var_names and "2" in adata.var_names and "3" in adata.var_names) or (1 in adata.var_names and 2 in adata.var_names and 3 in adata.var_names.astype(str)):
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
          print(f"Error processing {ID}: {str(e)}")
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
        core_action = False
        if include_new_samples and core_value == 2:
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
        core_action = False
        if include_new_samples and core_value == 2:
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
        core_action = False
        if include_new_samples and core_value == 2:
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
              adata.obs = adata.obs[["new_cell_type"]]
            except:
              try:
                adata.obs = adata.obs[["assignments"]]
                #rename to "new_cell_type"
                adata.obs.rename(columns={"assignments": "new_cell_type"}, inplace=True)

              except:
                adata.obs = adata.obs[["cell_type"]]
                #rename to "new_cell_type"
                adata.obs.rename(columns={"cell_type": "new_cell_type"}, inplace=True)

            #sc.pp.filter_genes(adata, min_cells=20)
            # Processing pipeline

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

            sc.pl.umap(adata, color=['new_cell_type'])

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
                  calc_umaps=True

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
                              calc_umaps=calc_umaps

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
                zipf.write(pdf_path, arcname=new_filename)
                print(f"Added {pdf_path} as {new_filename} to zip")
            else:
                print(f"PDF not found: {pdf_path}")

    print(f"Zip file created at: {output_zip_path}")


import pandas as pd
import scanpy as sc
import os
import gc
import zipfile
from concurrent.futures import ProcessPoolExecutor

def save_gene_batch_parquet(gene_list, sparse_X, obs_names, var_names):
    """
    Worker function: Receives only the essential data pieces.
    This avoids NameErrors and keeps RAM usage predictable.
    """
    # Find positions of requested genes
    gene_indices = [var_names.get_loc(g) for g in gene_list]
    
    # Extract columns and convert to dense only for this small batch
    X_batch = sparse_X[:, gene_indices].toarray()

    batch_df = pd.DataFrame(
        data=X_batch,
        index=obs_names,
        columns=gene_list
    )

    file_paths = []
    for gene in gene_list:
        file_path = f"genes_parquet/{gene}.parquet"
        # Using zstd for speed and efficiency
        batch_df[[gene]].to_parquet(file_path, compression="zstd")
        file_paths.append(file_path)

    del batch_df, X_batch
    gc.collect()
    return file_paths

# --- 2. The Main Function ---
def generate_atlas_umap(
    atlas_path='/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
    dataset_with_coords="/content/drive/MyDrive/pituitary_atlas/Unstructured/scanvi_all_cells_0828.h5ad",
    wt_dataset="/content/drive/MyDrive/pituitary_atlas/Unstructured/wt_adata_merged_post_scvi_0902.h5ad",
    mut_dataset="/content/drive/MyDrive/pituitary_atlas/Unstructured/mut_adata_merged_post_scvi_0902.h5ad",
    output_total_counts_csv="/content/drive/MyDrive/pituitary_atlas/Unstructured/total_counts.csv",
    zip_output_path="/content/drive/MyDrive/pituitary_atlas/Unstructured/adata_export_large_umap_1015.zip",
):
    # Load and process data
    adata = sc.read_h5ad(dataset_with_coords)
    adata_wt = sc.read_h5ad(wt_dataset)
    adata_mut = sc.read_h5ad(mut_dataset)

    # 1. Create IDs
    for d in [adata_wt, adata_mut, adata]:
        d.obs["cell_id"] = d.obs["SRA_ID"].astype(str) + "_" + d.obs["0"].astype(str)
    
    # Handle NaNs for parse samples
    nans = adata.obs["0"] == "nan"
    adata.obs.loc[nans, "cell_id"] = adata.obs.loc[nans, "SRA_ID"].astype(str) + "_" + \
                                     adata.obs.loc[nans, "bc"].astype(str) + "_" + \
                                     adata.obs.loc[nans, "subpool"].astype(str)

    # 2. Merge data using the new concat method
    merged_adata = sc.concat([adata_wt, adata_mut], label="batch")
    merged_adata.obs_names_make_unique() # STOPS THE WARNING AND RAM SPIKES
    merged_adata = merged_adata.copy()

    # 3. Map UMAPs
    mapping_umap1 = dict(zip(adata.obs["cell_id"], adata.obs["UMAP1"]))
    mapping_umap2 = dict(zip(adata.obs["cell_id"], adata.obs["UMAP2"]))
    merged_adata.obs["UMAP1"] = merged_adata.obs["cell_id"].map(mapping_umap1)
    merged_adata.obs["UMAP2"] = merged_adata.obs["cell_id"].map(mapping_umap2)
    merged_adata = merged_adata[~merged_adata.obs["UMAP1"].isna()].copy()

    # 4. Total Counts
    merged_adata.obs["ncounts"] = merged_adata.X.sum(axis=1).A1
    merged_adata.obs["ncounts"].to_csv(output_total_counts_csv)

    # 5. Metadata merge (your atlas logic)
    pituitary_atlas = pd.read_excel(atlas_path).drop_duplicates(subset=["SRA_ID"])
    pituitary_atlas = pituitary_atlas[["Comp_sex","10X version","Modality","Normal","SRA_ID"]]
    
    for col in ["Comp_sex","10X version","Modality","Normal","Sex"]:
        if col in merged_adata.obs.columns:
            merged_adata.obs.drop(col, axis=1, inplace=True)
            
    merged_adata.obs = merged_adata.obs.merge(pituitary_atlas, on="SRA_ID", how="left")

    # 6. Prepare for Export
    if hasattr(merged_adata.X, 'format') and merged_adata.X.format != 'csc':
        merged_adata.X = merged_adata.X.tocsc()

    obs_path = "obs.parquet"
    merged_adata.obs.to_parquet(obs_path, compression="gzip")
    merged_adata.obs[["UMAP1", "UMAP2"]].to_parquet("umap.parquet", compression="gzip")

    # 7. PARALLEL GENE EXPORT
    os.makedirs("genes_parquet", exist_ok=True)
    
    # We pass the sparse matrix and names ONLY. 
    # This uses a fraction of the RAM compared to passing the whole adata object.
    sparse_X = merged_adata.X
    obs_names = merged_adata.obs_names
    var_names = merged_adata.var_names
    
    gene_names = var_names.to_list()
    BATCH_SIZE = 100
    gene_batches = [gene_names[i:i + BATCH_SIZE] for i in range(0, len(gene_names), BATCH_SIZE)]

    # Limit workers in Colab to prevent OOM
    MAX_WORKERS = 2 

    print(f"Starting parallel export of {len(gene_names)} genes...")
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # We use a list comprehension to pass the specific objects to each task
        futures = [executor.submit(save_gene_batch_parquet, batch, sparse_X, obs_names, var_names) 
                   for batch in gene_batches]
        all_file_paths_lists = [f.result() for f in futures]

    gene_files = [path for sublist in all_file_paths_lists for path in sublist]

    # 8. Zip and Cleanup
    with zipfile.ZipFile(zip_output_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(obs_path)
        zf.write("umap.parquet")
        for file_path in gene_files:
            zf.write(file_path)

    # Cleanup temp files
    os.remove(obs_path)
    os.remove("umap.parquet")
    for f_path in gene_files:
        os.remove(f_path)
    os.rmdir("genes_parquet")

    print(f"Success! Saved to {zip_output_path}")



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
                          min_cells=50,
                          pdata_output_path_suffix="/analysis/pdata_assigned_0828.h5ad",
                          include_new_samples = False

                         ):

  df = pd.read_excel(atlas_source_path)
  #modality
  df = df[df['Modality'].isin(modality)]
  #reset
  df.reset_index(drop=True, inplace=True)

  missing = []

  for i in range(len(df)-1,-1,-1):
    SRA_ID = df.iloc[i]['SRA_ID']
    GEO = df.iloc[i]['GEO']
    ID = df.iloc[i][ID_col]
    core_value = df[df[ID_col]==ID]["Core"].values[0]
    core_action = False
    if include_new_samples and core_value == 2:
      core_action = True

    print(ID)

    author = df.iloc[i]['Author']
    normal = df.iloc[i]['Normal']

    #if normal == 2 then skip
    if normal == 2:
      continue
    #if ["n_cells"]<300 skip
    n_cells = df.iloc[i]['n_cells']
    if n_cells < min_cells_threshold:
      continue
    
    #if f"{nac_path}{ID}{pdata_output_path_suffix}" exists
    if os.path.exists(f"{nac_path}{ID}{pdata_output_path_suffix}"):
      print(f"Skipping {ID} as pdata already exists")
      continue

     #read adata assigned
    try:
      print(f"Processing {ID}...")
      try:
        adata_assigned = sc.read(f"{nac_path}{ID}{remaining_path}")
      except:
        print("Falling back to adata_assigned")
        adata_assigned = sc.read(f"{nac_path}{ID}{remaining_path_backup}")
    except:
       print(f" Missing entry for {ID}")
       missing.append(ID)
       continue
    
    adata_assigned.obs["sample"] = ID
    cell_counts = adata_assigned.X.sum(axis=1)

    import matplotlib.pyplot as plt
    plt.hist(cell_counts,color='red')
    plt.show()

    if filtering_based_on_barcodes and core_action==False:

      if author == "Rebboah et al. (2025)":

        adata_assigned.obs["barcode"] = [GEO + ":" + adata_assigned.obs["bc"][i] + "_" + adata_assigned.obs["subpool"][i] + "-1"  for i in range(len(adata_assigned.obs))]
        shape1 = adata_assigned.shape
        adata_assigned = adata_assigned[adata_assigned.obs["barcode"].isin(scvi_labels_to_keep[scvi_labels_to_keep["SRA_ID"]==SRA_ID]["barcode"].values)]
        shape2 = adata_assigned.shape
        print(f"Removed {shape1[0]-shape2[0]} cells")


      else:
        adata_assigned.obs["barcode"] = [GEO+":" + adata_assigned.obs['0'][i] + "-1" for i in range(len(adata_assigned.obs))]
        shape1 = adata_assigned.shape
        adata_assigned = adata_assigned[adata_assigned.obs["barcode"].isin(scvi_labels_to_keep[scvi_labels_to_keep["SRA_ID"]==SRA_ID]["barcode"].values)]
        shape2 = adata_assigned.shape
        print(f"Removed {shape1[0]-shape2[0]} cells")

    #if sample_col not in obs add it 
    if sample_col not in adata_assigned.obs.columns:
      adata_assigned.obs[sample_col] = ID

    
    pdata_assigned = dc.get_pseudobulk(
        adata_assigned,
        sample_col=sample_col,
        groups_col=groups_col,
        mode='sum',
        min_cells=50,
        min_counts=0
    )

    pdata_assigned.obs.columns = pdata_assigned.obs.columns.astype(str)
    pdata_assigned.var.columns = pdata_assigned.var.columns.astype(str)

    for col in pdata_assigned.obs.columns:
      pdata_assigned.obs[col] = pdata_assigned.obs[col].astype(str)
    pdata_assigned.write(f"{nac_path}{ID}{pdata_output_path_suffix}")
    print(f"ALL DONE for {ID}")
  
    
  print("Missing assignments files for the following SRA_IDs:")
  print(missing)





def read_merge_save_pbulk(
      atlas_source_path = '/content/drive/MyDrive/pituitary_atlas/source_table/pituitary_atlas.xlsx',
      output_path="/content/drive/MyDrive/pituitary_atlas/Unstructured/pdatas0828.h5ad",
      Modality = ['sc','sn','multi_rna'],
      total_n_cells_threshold = 200,
      nac_path="/content/drive/MyDrive/pituitary_atlas/processed/nac/",
      remaining_path="/analysis/pdata_assigned_0828.h5ad",
      cell_type_col='new_cell_type',
      backup_cell_type_col='assignments'
      ):
  
  df = pd.read_excel(atlas_source_path)
  #in df keep only sc sn or multi_rna
  df = df[df['Modality'].isin(Modality)]
  #reset
  df.reset_index(drop=True, inplace=True)

  pdatas= []
  for i in range(len(df)):
    normal = df.iloc[i]['Normal']
    #if normal == 2 then skip
    if normal == 2:
      continue
    #if ["n_cells"]<300 skip
    n_cells = df.iloc[i]['n_cells']
    if n_cells < total_n_cells_threshold:
      continue
    SRA_ID = df.iloc[i]['SRA_ID']
    print(SRA_ID)

    try:
      pdata_assigned = sc.read(f"{nac_path}{SRA_ID}{remaining_path}")
      pdatas.append(pdata_assigned)
    except:
      print(f"Error reading {SRA_ID}")
      continue

  print(len(pdatas))

  pdatas_merged = anndata.concat(pdatas, join='outer')
  #change nans to 0 in pdatas
  pdatas_merged.X[np.isnan(pdatas_merged.X)] = 0

  #turn ["psbulk_n_cells"] to ints
  pdatas_merged.obs["psbulk_n_cells"] = pdatas_merged.obs["psbulk_n_cells"].astype(float)
  pdatas_merged.obs["psbulk_n_cells"].sum()

  # 1. Fill 'new_cell_type' where it equals the string "nan"
  #if there is backup_cell_type_col use that to fill in
  if backup_cell_type_col not in pdatas_merged.obs.columns:
    pdatas_merged.obs[backup_cell_type_col] = pdatas_merged.obs[cell_type_col]
  mask_cell_type = pdatas_merged.obs[cell_type_col].isna()
  pdatas_merged.obs.loc[mask_cell_type, cell_type_col] = \
      pdatas_merged.obs.loc[mask_cell_type, backup_cell_type_col]

  # 2. Fill 'sample' where it equals the string "nan"
  mask_sample = pdatas_merged.obs["sample"].isna()
  pdatas_merged.obs.loc[mask_sample, "sample"] = \
      pdatas_merged.obs.loc[mask_sample, "SRA_ID"]

  #add metadata from df
  #only keep sample and new_cell_type
  pdatas_merged.obs =pdatas_merged.obs[['sample',cell_type_col,"psbulk_n_cells"]]
  #rename SRA_ID to sample
  pdatas_merged.obs = pd.merge(pdatas_merged.obs, df, left_on='sample', right_on='SRA_ID')

  #from pdatas.obs remove columns that contain "Unnamed" or _index
  pdatas_merged.obs = pdatas_merged.obs.loc[:, ~pdatas_merged.obs.columns.str.contains('Unnamed')]
  pdatas_merged.obs = pdatas_merged.obs.loc[:, ~pdatas_merged.obs.columns.str.contains('_index')]

  pdatas_merged.obs.columns = pdatas_merged.obs.columns.astype(str)
  pdatas_merged.var.columns = pdatas_merged.var.columns.astype(str)
  pdatas_merged.var.index = pdatas_merged.var.index.astype(str)
  pdatas_merged.obs.index= pdatas_merged.obs.index.astype(str)

  # Convert all columns in obs to strings if needed
  for col in pdatas_merged.obs.columns:
      if pdatas_merged.obs[col].apply(type).nunique() > 1:
          pdatas_merged.obs[col] = pdatas_merged.obs[col].astype(str)

  # Now, try saving the AnnData object
  pdatas_merged.write(f"{output_path}")














