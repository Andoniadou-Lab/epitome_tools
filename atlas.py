import subprocess
import os
import pandas as pd
import subprocess
import json
from scipy.io import mmread, mmwrite
import numpy as np
from scipy.sparse import csc_matrix
from scipy.stats import median_abs_deviation as mad
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tempfile
import scipy.sparse as sp
import scanpy as sc
import time

def sra_downloader(raw_path,df):
    """
    Downloads SRA datasets from NCBI using wget and ffq.
    Args:
        raw_path (str): Path to the directory where datasets will be downloaded.
        df (pd.DataFrame): DataFrame containing SRA IDs and experiment types.
    """

    datasets_already_downloaded = os.listdir(raw_path)

    #iterating through the table and getting datasets
    for i in range(len(df)-1,-1,-1):
        SRA_ID = df.iloc[i]['SRA_ID']
        experiment_type = df.iloc[i]['Modality']

        if SRA_ID in datasets_already_downloaded:
            print("SRA_ID already downloaded")
            continue

        if (not os.path.exists(raw_path+ SRA_ID)) and ((experiment_type == "sc") or (experiment_type == "sn")):
            os.mkdir(raw_path + SRA_ID)

            print("Folder created")

        target_folder = raw_path + SRA_ID
        download_command = f"wget -P {target_folder} /path/to/target/folder $(ffq --ncbi {SRA_ID} | jq -r '.[] | .url' | tr '\n' ' ')"
        #run command
        print(download_command)
        subprocess.run(download_command, shell=True)


def fastq_downloader(raw_path, df):
    """
    Downloads fastq files from SRA datasets.
    Args:
        raw_path (str): Path to the directory where datasets are stored.
        df (pd.DataFrame): DataFrame containing SRA IDs and experiment types.
    """

    for i in range(len(df)-1,-1,-1):
        SRA_ID = df.iloc[i]['SRA_ID']
        experiment_type = df.iloc[i]['Modality']
        outdir = raw_path + SRA_ID
        files = os.listdir(outdir)
        for file in files:
            print(file)
            if ("lite" in file) or ("sra" in file):
                #check if there are fastq files with starting with file
                if any(x.startswith(file) and x.endswith('.fastq.gz') for x in files):
                    print("fastq files found")
                    #or it has bam string in it
                elif "bam" in file:
                    print("bam file found")
                else:
                    file = raw_path + SRA_ID + '/' + file
                    fastqdump_command= f"fastq-dump --outdir {outdir} --gzip --split-files {file}"
                    subprocess.run(fastqdump_command, shell=True)


def run_kb_count_nac(nac_path,raw_path, df, species="mouse"):
    """ Runs kb count for nascent RNA datasets.
    Args:
        nac_path (str): Path to the directory where nascent RNA datasets will be stored.
        raw_path (str): Path to the directory where raw datasets are stored.
        df (pd.DataFrame): DataFrame containing SRA IDs, experiment types, and technologies.
    """

    kbref_command_nac = f"kb ref --workflow=nac  -d {species} -i index_nac.idx -g t2g_nac.txt -g t2g_nac.txt -c1 cdna.txt -c2 nascent.txt -f1 cdna.fasta -f2 nascent.fasta --overwrite"
    kbref_command_std = f"kb ref -d {species} -i index_std.idx -g t2g_std.txt --overwrite"

    #if /index_nac.idx doesnt exist, run kb ref
    if not os.path.exists("index_nac.idx"):
        print("Running kb ref for nac")
        subprocess.run(kbref_command_nac, shell=True)

    if not os.path.exists("index_std.idx"):
        print("Running kb ref for std")
        subprocess.run(kbref_command_std, shell=True)

    datasets_already_downloaded = os.listdir(nac_path)
    
    #iterating through the table and getting datasets
    for i in range(len(df)):
        SRA_ID = df.iloc[i]['SRA_ID']
        experiment_type = df.iloc[i]['Modality']
        tech = df.iloc[i]['10X version']
        author = df.iloc[i]['Author']

        #print all of these
        print("""
        SRA_ID: {}
        experiment_type: {}
        tech: {}
        """.format(SRA_ID,experiment_type,tech))

        #check if SRA_ID is already downloaded
        if SRA_ID in datasets_already_downloaded:
            print("SRA_ID already aligned for nac")
            continue
        
        all_fastqs = find_fastqs(raw_path, SRA_ID)
        #if there arent at least 2 fastqs then continue
        if len(all_fastqs) < 2:
            print("Not enough fastqs")
            continue

        if not os.path.exists(nac_path + SRA_ID):
            os.mkdir(nac_path + SRA_ID)

        if experiment_type == "sn":
            print("Running kb count")
            output_dir = f"{nac_path}{SRA_ID}/"
            kbcount_nac_command1 = f"kb count -i index_nac.idx -g t2g_nac.txt -c1 cdna.txt -c2 nascent.txt  -x {tech} --workflow=nac --sum=total --overwrite --h5ad -o {output_dir} {' '.join(all_fastqs)}"
            print(kbcount_nac_command1)
            subprocess.run(kbcount_nac_command1, shell=True)

        elif (experiment_type == "sc"):
            print("Running kb count")
            output_dir = f"{nac_path}{SRA_ID}/"
            kbcount_nac_command1 = f"kb count -i index_nac.idx -g t2g_nac.txt -c1 cdna.txt -c2 nascent.txt -x {tech} --workflow=nac --sum=total --overwrite --h5ad -o {output_dir} {' '.join(all_fastqs)}"
            print(kbcount_nac_command1)
            subprocess.run(kbcount_nac_command1, shell=True)
            
            #change nac to tcc
            output_dir_tcc = nac_path.replace("nac", "tcc")
            #make sure it exists
            if not os.path.exists(output_dir_tcc):
                os.mkdir(output_dir_tcc)
            kbcount_tcc_command = f"kb count -i index_std.idx -g t2g_std.txt --overwrite -x {tech} --tcc --h5ad -o {output_dir_tcc} {' '.join(all_fastqs)}"
            print(kbcount_tcc_command)
            subprocess.run(kbcount_tcc_command, shell=True)

        elif (experiment_type == "multi_rna"):
            print("Running kb count")
            output_dir = f"{nac_path}{SRA_ID}/"
            #from all_fastqs remove ones that contain "atac"
            all_fastqs = [x for x in all_fastqs if "atac" not in x]
            print(all_fastqs)
            kbcount_nac_command1 = f"kb count -i index_nac.idx -g t2g_nac.txt -c1 cdna.txt -c2 nascent.txt -x {tech} -w /10xMultiome/gex_737K-arc-v1.txt --workflow=nac --sum=total --overwrite --h5ad -o {output_dir} {' '.join(all_fastqs)}"
            print(kbcount_nac_command1)
            subprocess.run(kbcount_nac_command1, shell=True)

def find_fastqs(raw_path, SRA_ID):
    """ Finds fastq files in the specified directory for a given SRA ID.
    Args:
        raw_path (str): Path to the directory where datasets are stored.
        SRA_ID (str): SRA ID for which fastq files are to be found.
    Returns:
        list: List of fastq file paths.
    """

    directory = raw_path + SRA_ID
    files = os.listdir(directory)
    all_fastqs = []
    for file in files:
        if file.endswith('.fastq.gz'):
            all_fastqs.append(directory + '/' + file)

    lite_name = True if any(x for x in all_fastqs if 'lite' in x) else False
    if lite_name:
        id_tags = [x.split('lite')[0] for x in all_fastqs]

    SRR_name = True if any(x for x in all_fastqs if 'SRR' in x) else False
    if SRR_name:
        id_tags = [x.split('.1_')[0] for x in all_fastqs]

        if id_tags[0].count('/') > 1:
            id_tags = [x.split('_1')[0] for x in all_fastqs]

    #if neither
    elif (not lite_name) and (not SRR_name):
        id_tags = [x.split('_R1_')[0] for x in all_fastqs if '_R1_' in x]

    print(id_tags)
    all_fastqs_final = []
    unique_tags= list(set(id_tags))
    for tag in unique_tags:
        fastqs_with_tag = [x for x in all_fastqs if tag in x]
        if len(fastqs_with_tag) > 1:
            #sort by size and take the two largest
            fastqs_with_tag.sort(key=lambda x: os.path.getsize(x), reverse=True)
            all_fastqs_final.append(fastqs_with_tag[1])
            all_fastqs_final.append(fastqs_with_tag[0])

    all_fastqs = all_fastqs_final
    return all_fastqs


def process_matrix_files(matrix_file, barcodes_file,matrix_file_new,barcodes_file_new,num=10):
    # Read matrix data
    # Load matrix data
    matrix_data = mmread(matrix_file).tocsc()  # Load matrix as compressed sparse column format

    # Read barcodes
    barcodes = pd.read_csv(barcodes_file, header=None)
    #squeeze
    barcodes = barcodes.squeeze()

    # Compute the sum of counts for each cell (column)
    cell_sums = matrix_data.sum(axis=1).A1

    # Filter cells where the total counts are greater than 10,000
    filtered_cell_indices = cell_sums > num

    filtered_matrix_data = matrix_data[filtered_cell_indices, :]
    filtered_barcodes = barcodes[filtered_cell_indices]

    # Save the filtered matrix
    mmwrite(matrix_file_new, csc_matrix(filtered_matrix_data))

    # Save the filtered barcodes
    filtered_barcodes.to_csv(barcodes_file_new, header=False, index=False)

    # Load Matrix Market file as DataFrame
    matrix_data_df = pd.read_csv(matrix_file_new, sep=' ', header=2)

    # Sort by first then second column
    first_col = matrix_data_df.columns[0]
    second_col = matrix_data_df.columns[1]
    matrix_data_sorted = matrix_data_df.sort_values(by=[first_col, second_col])

    # Write back to .mtx file with the same %header as before
    # Write header lines
    with open(matrix_file_new, 'r') as f:
        header_lines = [f.readline() for _ in range(2)]  # Read the first three lines

    # Write sorted data to the matrix file
    with open(matrix_file_new, 'w') as f:
        f.writelines(header_lines)
        matrix_data_sorted.to_csv(f, sep=' ', index=False, header=True, float_format='%.6g')

    print(f"Sorted Matrix Market file saved as: {matrix_file_new}")


def filter_top_cells_mtx(matrix_file, barcodes_file, output_matrix_file, output_barcodes_file, top_n=50000):
    """
    Filter MTX matrix and corresponding barcode file to keep only the top N cells by total counts.
    Filters rows of the matrix (cells) while keeping all columns (genes).

    Parameters:
    -----------
    matrix_file : str
        Path to input MTX file
    barcodes_file : str
        Path to input barcodes file
    output_matrix_file : str
        Path to save filtered MTX file
    output_barcodes_file : str
        Path to save filtered barcodes file
    top_n : int, optional
        Number of top cells to keep (default: 50000)

    Returns:
    --------
    tuple
        (filtered_matrix, filtered_barcodes)
    """
    # Load the matrix in CSR format for efficient row operations
    matrix = mmread(matrix_file).tocsr()

    # Read barcodes
    barcodes = pd.read_csv(barcodes_file, header=None).squeeze()

    # Calculate row sums (total counts per cell)
    row_sums = np.array(matrix.sum(axis=1)).flatten()

    # Get indices of top N cells by total counts
    top_cell_indices = np.argsort(row_sums)[-top_n:]
    top_cell_indices.sort()  # Sort indices to maintain original order

    # Filter matrix and barcodes
    filtered_matrix = matrix[top_cell_indices, :]
    filtered_barcodes = barcodes.iloc[top_cell_indices]

    # Save filtered matrix
    mmwrite(output_matrix_file, filtered_matrix)

    # Save filtered barcodes
    filtered_barcodes.to_csv(output_barcodes_file, header=False, index=False)

    # Sort the matrix file if needed
    matrix_data_df = pd.read_csv(output_matrix_file, sep=' ', header=2)

    # Sort by first then second column
    first_col = matrix_data_df.columns[0]
    second_col = matrix_data_df.columns[1]
    matrix_data_sorted = matrix_data_df.sort_values(by=[first_col, second_col])

    # Write header lines and sorted data
    with open(output_matrix_file, 'r') as f:
        header_lines = [f.readline() for _ in range(2)]

    with open(output_matrix_file, 'w') as f:
        f.writelines(header_lines)
        matrix_data_sorted.to_csv(f, sep=' ', index=False, header=True, float_format='%.6g')

    return filtered_matrix, filtered_barcodes

def filtering_cells(df,raw_path="/content/drive/MyDrive/pituitary_atlas",root_path = "/content/",atlas_source_file = "pituitary_atlas.xlsx",force=False):
    max_retries = 3
    for i in range(len(df)-1,-1,-1):
        SRA_ID = df.iloc[i]['SRA_ID']
        if os.path.exists(f"{raw_path}/processed/nac/{SRA_ID}/analysis"):
            print(SRA_ID)
            try:
                #see if {raw_path}/processed/nac/{SRA_ID}/counts_unfiltered/ exists or continue
                if not os.path.exists(f"{raw_path}/processed/nac/{SRA_ID}/counts_unfiltered/"):
                    print(f"{raw_path}/processed/nac/{SRA_ID}/counts_unfiltered/ does not exist, skipping.")
                    continue
                barcodes= f"{raw_path}/processed/nac/{SRA_ID}/counts_unfiltered/cells_x_genes.barcodes.txt"
                matrix_file = f"{raw_path}/processed/nac/{SRA_ID}/counts_unfiltered/cells_x_genes.total.mtx"
                genes = f"{raw_path}/processed/nac/{SRA_ID}/counts_unfiltered/cells_x_genes.genes.names.txt"
                barcodes_out = f"{raw_path}/processed/nac/{SRA_ID}/analysis/out_42_cells_x_genes.barcodes.txt"
                output_mat = f"{raw_path}/processed/nac/{SRA_ID}/analysis/out_42_cells_x_genes.total.mtx"

                #if output_mat exists, continue
                if not force:
                    if os.path.exists(output_mat):
                        #set filtering_junk to 1
                        atlas =  pd.read_excel(f'{raw_path}/source_table/{atlas_source_file}')
                        atlas.loc[atlas['SRA_ID'] == SRA_ID, 'filtering_junk'] = 1
                        #save back
                        atlas.to_excel(f'{raw_path}/source_table/{atlas_source_file}', index=False)
                        print(f"{output_mat} already exists, skipping filtering.")
                        continue
                i =0

                while (not os.path.exists(barcodes_out) or not os.path.exists(output_mat)) and i < max_retries:
                    filter_top_cells_mtx(matrix_file, barcodes,output_mat,barcodes_out,40000)
                    #sleep 30 sec to allow for files to sync with google drive
                    time.sleep((i+1)*30)
                    i += 1
                
                
                #mx filter
                barcodes_out2 = f"{raw_path}/processed/nac/{SRA_ID}/analysis/out_43_cells_x_genes.barcodes.txt"
                output_mat2 = f"{raw_path}/processed/nac/{SRA_ID}/analysis/out_43_cells_x_genes.total.mtx"

                mxfilter_command = f"mx filter -c 2 -bi {barcodes_out} -bo {barcodes_out2} -o {output_mat2} {output_mat}"
                #run command
                for attempt in range(1, max_retries + 1):
                    result = subprocess.run(mxfilter_command, shell=True)

                    if result.returncode == 0:
                        print(f"Command succeeded on attempt {attempt}")
                        break
                    else:
                        print(f"Attempt {attempt} failed with return code {result.returncode}")
                else:
                    print("All attempts failed.")
                
                #inspect
                inspect_command = f"mx inspect -a all -o {root_path}inspect_check1.json -gi {genes} -bi {barcodes_out2} {output_mat2}"
                time.sleep(30) #sleep 30 sec to allow for files to sync with google drive
                subprocess.run(inspect_command, shell=True)
                time.sleep(30) #sleep 30 sec to allow for files to sync with google drive
                
                inspect=json.load(open(f'{root_path}inspect_check1.json'))
                ncells= inspect["ncells"]

                barcodes_out3 = f"{raw_path}/processed/nac/{SRA_ID}/analysis/out_44_cells_x_genes.barcodes.txt"
                output_mat3 = f"{raw_path}/processed/nac/{SRA_ID}/analysis/out_44_cells_x_genes.total.mtx"

                if ncells < 3000:
                    #sleep 30 sec to allow for files to sync with google drive
                    time.sleep(30) #sleep 30 sec to allow for files to sync with google drive
                    process_matrix_files(matrix_file, barcodes,output_mat3,barcodes_out3,num=1500)
                    time.sleep(30) #sleep 30 sec to allow for files to sync with google drive
                else:
                    i = 0
                    while not os.path.exists(barcodes_out3) or not os.path.exists(output_mat3) or i < max_retries:
                        #filter top cells
                        time.sleep(30) #sleep 30 sec to allow for files to sync with google drive
                        filter_top_cells_mtx(matrix_file, barcodes,output_mat3,barcodes_out3,int(int(ncells)+min(0.3*int(ncells),2000)))
                        time.sleep(30) #sleep 30 sec to allow for files to sync with google drive
                        i += 1

                barcodes_out4 = f"{raw_path}/processed/nac/{SRA_ID}/analysis/out_45_cells_x_genes.barcodes.txt"

                output_mat4 = f"{raw_path}/processed/nac/{SRA_ID}/analysis/out_45_cells_x_genes.total.mtx"

                genes_out = f"{raw_path}/processed/nac/{SRA_ID}/analysis/out_45_cells_x_genes.genes.names.txt"

                #if these all exist, print that we have successfully produced the output
                if os.path.exists(output_mat4) and os.path.exists(barcodes_out4) and os.path.exists(genes_out):
                    print(f"Successfully produced output for {SRA_ID}")

                if not force:
                    #if output_mat4  exist continue
                    if os.path.exists(output_mat4):
                        continue

                mxclean_command = f"mx clean -gi {genes} -go {genes_out} -bi {barcodes_out3} -bo {barcodes_out4} -o {output_mat4} {output_mat3}"
                for attempt in range(1, max_retries + 1):
                    time.sleep(30) #sleep 30 sec to allow for files to sync with google drive
                    result = subprocess.run(mxclean_command, shell=True)
                    time.sleep(30) #sleep 30 sec to allow for files to sync with google drive
                    if result.returncode == 0:
                        print(f"Command succeeded on attempt {attempt}")
                        break
                    else:
                        print(f"Attempt {attempt} failed with return code {result.returncode}")
                else:
                    print("All attempts failed.")

                #inspect
                
                inspect_command = f"mx inspect -a all -o {root_path}inspect_check1.json -gi {genes_out} -bi {barcodes_out4} {output_mat4}"
                time.sleep(30)
                subprocess.run(inspect_command, shell=True)
                time.sleep(30)
                #read json
                
                inspect=json.load(open(f'{root_path}inspect_check1.json'))

                #print ncells
                print("We have this many cells: ", inspect["ncells"])
                print("We have this many genes: ", inspect["ngenes"])

                #set filtering_junk to 1
                atlas =  pd.read_excel(f'{raw_path}/source_table/{atlas_source_file}')
                atlas.loc[atlas['SRA_ID'] == SRA_ID, 'filtering_junk'] = 1
                atlas.to_excel(f'{raw_path}/source_table/{atlas_source_file}', index=False)

            except Exception as e:
                #print explicit error
                print(f"Error for {SRA_ID}: {e}")


def qc1(adata,species,min_genes=800, min_counts=1000):
    #you could also use a whitelist of barcodes from the filtered barcodes for each sample - not even sure what I meant by this
    if species == "mouse":
      mt_tag="mt-"
      ribo_tag=("Rps","Rpl")
    elif species == "human":
      mt_tag="MT-"
      ribo_tag=("RPS","RPL")
    sc.pp.filter_cells(adata, min_genes = min_genes)
    sc.pp.filter_cells(adata, min_counts = min_counts)
    adata.var["mt"] = adata.var_names.str.startswith(mt_tag)
    adata.var["ribo"] = adata.var_names.str.startswith(ribo_tag)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, percent_top=[20], log1p=True)

    return adata

def mad_outlier(adata, metric, nmads, upper_only = False,value = False):
    if value and upper_only:
      upper_val = np.median(adata.obs[metric]) + nmads * mad(adata.obs[metric])
      return [upper_val,0]
    elif value and not upper_only:
      upper_val = np.median(adata.obs[metric]) + nmads * mad(adata.obs[metric])
      lower_val = np.median(adata.obs[metric]) - nmads * mad(adata.obs[metric])
      return [upper_val,lower_val]
    else:
      M = adata.obs[metric]

      if not upper_only:
          return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))

    return (M > np.median(M) + nmads * mad(M))

def violinplot(ax, adata, col, vals):
    values = adata.obs[col].values
    sns.violinplot(y=values, ax=ax)
    for v in vals:
        ax.axhline(v, color='red')
    ax.set_title(col)
    ax.set_ylabel('')

def qc2(adata, species, experiment="sc"):
    if experiment == "sc":
        if species == "mouse":
            malat_tag = "Malat1"
        elif species == "human":
            malat_tag = "MALAT1"

        adata.var["malat"] = adata.var_names.str.startswith(malat_tag)

        sc.pp.calculate_qc_metrics(adata, qc_vars=["malat"], inplace=True, percent_top=[20], log1p=True)

        mt_val = mad_outlier(adata, 'pct_counts_mt', 5, upper_only=True, value=True)
        ribo_val = mad_outlier(adata, 'pct_counts_ribo', 5, upper_only=True, value=True)
        malat_val = mad_outlier(adata, 'pct_counts_malat', 5, upper_only=True, value=True)
        log1p_n_genes_val = mad_outlier(adata, 'log1p_n_genes_by_counts', 4, upper_only=False, value=True)
        log1p_total_val = mad_outlier(adata, 'log1p_total_counts', 4, upper_only=False, value=True)
        top_genes_val = mad_outlier(adata, 'pct_counts_in_top_20_genes', 5, upper_only=True, value=True)

        # Create a 2x3 grid of subplots
        fig, axs = plt.subplots(2, 3, figsize=(15, 10))
        axs = axs.flatten()

        violinplot(axs[0], adata, 'pct_counts_mt', mt_val + [25])
        violinplot(axs[1], adata, 'pct_counts_ribo', ribo_val + [30])
        violinplot(axs[2], adata, 'pct_counts_malat', malat_val)
        violinplot(axs[3], adata, 'log1p_n_genes_by_counts', log1p_n_genes_val)
        violinplot(axs[4], adata, 'log1p_total_counts', log1p_total_val)
        violinplot(axs[5], adata, 'pct_counts_in_top_20_genes', top_genes_val)

        plt.tight_layout()

        adata = adata[adata.obs.pct_counts_mt < 25]
        adata = adata[adata.obs.pct_counts_ribo < 30]

        bool_vector = mad_outlier(adata, 'log1p_total_counts', 4) + \
                     mad_outlier(adata, 'log1p_n_genes_by_counts', 4) + \
                     mad_outlier(adata, 'pct_counts_in_top_20_genes', 5, upper_only=True) + \
                     mad_outlier(adata, 'pct_counts_mt', 5, upper_only=True) + \
                     mad_outlier(adata, 'pct_counts_ribo', 5, upper_only=True) + \
                     mad_outlier(adata, 'pct_counts_malat', 5, upper_only=True)

    else:
        mt_val = mad_outlier(adata, 'pct_counts_mt', 5, upper_only=True, value=True)
        ribo_val = mad_outlier(adata, 'pct_counts_ribo', 5, upper_only=True, value=True)
        log1p_n_genes_val = mad_outlier(adata, 'log1p_n_genes_by_counts', 4, upper_only=False, value=True)
        log1p_total_val = mad_outlier(adata, 'log1p_total_counts', 4, upper_only=False, value=True)
        top_genes_val = mad_outlier(adata, 'pct_counts_in_top_20_genes', 5, upper_only=True, value=True)

        # Create a 2x3 grid of subplots (5 plots, one empty)
        fig, axs = plt.subplots(2, 3, figsize=(15, 10))
        axs = axs.flatten()

        violinplot(axs[0], adata, 'pct_counts_mt', mt_val + [5])
        violinplot(axs[1], adata, 'pct_counts_ribo', ribo_val + [10])
        violinplot(axs[2], adata, 'log1p_n_genes_by_counts', log1p_n_genes_val)
        violinplot(axs[3], adata, 'log1p_total_counts', log1p_total_val)
        violinplot(axs[4], adata, 'pct_counts_in_top_20_genes', top_genes_val)
        axs[5].axis('off')  # Turn off the empty subplot

        plt.tight_layout()

        adata = adata[adata.obs.pct_counts_ribo < 10]
        adata = adata[adata.obs.pct_counts_mt < 5]
        bool_vector = mad_outlier(adata, 'log1p_total_counts', 4) + \
                     mad_outlier(adata, 'log1p_n_genes_by_counts', 4) + \
                     mad_outlier(adata, 'pct_counts_in_top_20_genes', 5)

    adata = adata[~bool_vector]
    adata.uns['cells_removed'] = sum(bool_vector)

    return adata, fig

def kneeplot(adata_filtered_final, adata_filtered_init, adata_raw):
    """
    Create a publication-quality knee plot to visualize cell filtering in single-cell RNA sequencing data.

    This function generates a scatter plot showing the relationship between cell rank and
    total UMI counts per cell. The plot follows publication-quality standards suitable for
    journals like Nature, with rank on the x-axis and UMI counts on the y-axis.

    Parameters
    ----------
    adata_filtered_final : AnnData
        Final filtered AnnData object containing selected cells
    adata_filtered_init : AnnData
        Initially filtered AnnData object
    adata_raw : AnnData
        Raw (unfiltered) AnnData object

    Returns
    -------
    tuple
        A tuple containing:
        - matplotlib.figure.Figure: Publication-quality knee plot
        - pandas.DataFrame: DataFrame containing the plotting data and cell classifications

    Notes
    -----
    The plot shows log10-transformed values for both axes:
    - x-axis: log10(rank of cells by total counts)
    - y-axis: log10(total counts per cell)
    """
    filtered_final_obs_names = adata_filtered_final.obs_names
    filtered_init_obs_names = adata_filtered_init.obs_names

    # Ranking barcodes
    adata_raw.obs["sum"] = adata_raw.X.sum(1)

    sum_obs_names_df = pd.DataFrame(adata_raw.obs["sum"])
    sum_obs_names_df.index = adata_raw.obs_names
    # Sort
    sum_obs_names_df = sum_obs_names_df.sort_values(by="sum", ascending=False)  # Changed to descending
    sum_obs_names_df["rank"] = range(1, len(sum_obs_names_df) + 1)  # Start rank from 1
    sum_obs_names_df["log10_sum"] = np.log10(sum_obs_names_df["sum"])
    sum_obs_names_df["log10_rank"] = np.log10(sum_obs_names_df["rank"])
    sum_obs_names_df["filtered_init"] = 0
    sum_obs_names_df.loc[filtered_init_obs_names, "filtered_init"] = 1
    sum_obs_names_df["filtered_final"] = 0
    sum_obs_names_df.loc[filtered_final_obs_names, "filtered_final"] = 1

    # Set the style for publication-quality figure
    plt.style.use('default')
    fig = plt.figure(figsize=(10, 6), dpi=300)

    # Create GridSpec for the main plot
    gs = GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])

    # Plot with Nature-style formatting
    # Plot black dots for all cells
    ax.scatter(sum_obs_names_df["log10_rank"],
              sum_obs_names_df["log10_sum"],
              color="black",
              s=10,
              alpha=0.3,
              label="All barcodes",
              rasterized=True)  # Rasterize for smaller file size

    # Plot thick light blue dots for initially filtered cells
    init_mask = sum_obs_names_df["filtered_init"] == 1
    ax.scatter(sum_obs_names_df.loc[init_mask, "log10_rank"],
              sum_obs_names_df.loc[init_mask, "log10_sum"],
              color="#4878CF",  # Nature-style blue
              s=100,
              alpha=0.5,
              label="Initial filtering",
              rasterized=True)

    # Plot small red dots for final filtered cells
    final_mask = sum_obs_names_df["filtered_final"] == 1
    ax.scatter(sum_obs_names_df.loc[final_mask, "log10_rank"],
              sum_obs_names_df.loc[final_mask, "log10_sum"],
              color="#E24A33",  # Nature-style red
              s=10,
              alpha=0.85,
              label="Final filtering",
              rasterized=True)

    # Customize the plot for publication quality
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Set font sizes
    SMALL_SIZE = 10
    MEDIUM_SIZE = 12
    BIGGER_SIZE = 14

    # Add labels and legend
    ax.set_xlabel("log10(Barcode rank)", fontsize=MEDIUM_SIZE)
    ax.set_ylabel("log10(Total UMI counts)", fontsize=MEDIUM_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=SMALL_SIZE)

    # Add color legend at top right
    leg = ax.legend(frameon=False, fontsize=SMALL_SIZE, markerscale=3,
                   bbox_to_anchor=(1.02, 1), loc='upper right')

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Add gridlines
    ax.grid(True, linestyle='--', alpha=0.3, which='major')

    # Add counts in the corner
    total_cells = len(sum_obs_names_df)
    init_cells = sum_obs_names_df["filtered_init"].sum()
    final_cells = sum_obs_names_df["filtered_final"].sum()

    stats_text = (f'n(total) = {total_cells:,}\n'
                 f'n(initial) = {int(init_cells):,}\n'
                 f'n(final) = {int(final_cells):,}')

    # Add text box with cell counts at bottom left
    plt.text(0.05, 0.05, stats_text,
             transform=ax.transAxes,
             verticalalignment='bottom',
             horizontalalignment='left',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
             fontsize=SMALL_SIZE)

    return fig, sum_obs_names_df

def create_assignment_barchart(assignments):
    """
    Create a bar chart of cell type assignments and return the figure.

    Parameters:
    -----------
    assignments : list or pandas.Series
        List of cell type assignments

    Returns:
    --------
    matplotlib.figure.Figure
        Figure containing the bar chart
    """
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    pd.Series(assignments).value_counts().plot(kind='bar', stacked=True, ax=ax)

    plt.title('Initial Cellassign Cell Types')
    plt.xlabel('Cell Type')
    plt.ylabel('Number of Cells')
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()

    return fig

def make_summary_pdf(atlas_df, SRA_ID, fig_kneeplot, violinplots, fig_umap, fig_dotplot, raw_dotplot, fig_barchart,sparsity_original,sparsity,res_sparsity, path="/content/drive/MyDrive/pituitary_atlas"):

    """
    Creates a more compact layout with larger knee plot
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Save individual plots
        knee_path = os.path.join(tmpdir, 'knee.png')
        violin_path = os.path.join(tmpdir, 'violin.png')
        bar_path = os.path.join(tmpdir, 'bar.png')
        umap_path = os.path.join(tmpdir, 'umap.png')
        dot_path = os.path.join(tmpdir, 'dot.png')
        raw_dot_path = os.path.join(tmpdir, 'raw_dot.png')

        fig_kneeplot.savefig(knee_path, bbox_inches='tight', dpi=500)
        violinplots.savefig(violin_path, bbox_inches='tight', dpi=500)
        fig_barchart.savefig(bar_path, bbox_inches='tight', dpi=500)
        fig_umap.savefig(umap_path, bbox_inches='tight', dpi=500)
        fig_dotplot.savefig(dot_path, bbox_inches='tight', dpi=400)
        raw_dotplot.savefig(raw_dot_path, bbox_inches='tight', dpi=400)

        # Create summary figure - reduced height
        fig = plt.figure(figsize=(20, 24))
        # 4 rows, reduced height ratios and spacing
        gs = GridSpec(5, 5, height_ratios=[0.25, 0.35, 0.35, 0.4,0.4], figure=fig)
        gs.update(wspace=0.04, hspace=0.08)  # Reduced spacing further

        # Metadata
        meta_ax = fig.add_subplot(gs[0, :])
        meta_ax.axis('off')
        meta_row = atlas_df[atlas_df['SRA_ID'] == SRA_ID].iloc[0]
        metadata_text = (
            f"Sample Information:\n"
            f"==================\n"
            f"Author: {meta_row['Author']}\n"
            f"SRA_ID: {meta_row['SRA_ID']}\n"
            f"GEO: {meta_row['GEO']}\n"
            f"Age_numeric: {meta_row['Age_numeric']}\n"
            f"Name: {meta_row['Name']}\n"
            f"Comp_sex: {meta_row['Comp_sex']}\n"
            f"Number of cells: {meta_row['n_cells']}\n"
            f"Initial sparsity: {sparsity_original}\n"
            f"Resulting sparsity: {res_sparsity}\n"
            f"Median CellAssign probability: {meta_row['median_cellassign_prob']:.3f}"
        )
        meta_ax.text(0.02, 0.95, metadata_text,
                    fontsize=10,
                    verticalalignment='top',
                    fontfamily='monospace',
                    bbox=dict(facecolor='white',
                             edgecolor='black',
                             alpha=0.8))

        # Violin plot and Knee plot row
        violin_ax = fig.add_subplot(gs[1, :3])  # Violin takes 3/5 of width
        violin_ax.imshow(plt.imread(violin_path))
        violin_ax.axis('off')

        knee_ax = fig.add_subplot(gs[1, 3:])  # Knee takes 2/5 of width
        knee_ax.imshow(plt.imread(knee_path))
        knee_ax.axis('off')

        # UMAP and Bar plot row
        umap_ax = fig.add_subplot(gs[2, :3])  # UMAP takes 3/5 of width
        umap_ax.imshow(plt.imread(umap_path))
        umap_ax.axis('off')

        bar_ax = fig.add_subplot(gs[2, 3:])  # Bar takes 2/5 of width
        bar_ax.imshow(plt.imread(bar_path))
        bar_ax.axis('off')

        #raw dotplot

        raw_dot_ax = fig.add_subplot(gs[3, :])
        raw_dot_ax.imshow(plt.imread(raw_dot_path))
        raw_dot_ax.axis('off')

        # Dot plot
        dot_ax = fig.add_subplot(gs[4, :])
        dot_ax.imshow(plt.imread(dot_path))
        dot_ax.axis('off')

        # Save final figure
        plt.savefig(f"{path}/processed/nac/{SRA_ID}/analysis/summary.pdf",
                    bbox_inches='tight', dpi=400, pad_inches=0)  # Reduced padding
        plt.close('all')

def get_sparsity(adata):
  if sp.issparse(adata.X):
    zero_entries = (adata.X == 0).sum()
    non_zero_entries = (adata.X != 0).sum()
    sparsity = zero_entries / (non_zero_entries + zero_entries)
  else:
    zero_entries = np.sum(adata.X == 0)
    non_zero_entries = np.sum(adata.X != 0)
    sparsity = zero_entries / (non_zero_entries + zero_entries)
  return sparsity

def shrink_high_values_normalized_df(df, column_name="ambient_profile_mrna", percentile=50, shrinkage_factor=0.5):
    """
    Shrink values above the specified percentile toward that percentile value and normalize
    the column so it sums to 1, while keeping zero values and lower values unchanged.

    Args:
        df (pd.DataFrame): A Pandas DataFrame containing the data.
        column_name (str): The name of the column to process.
        percentile (float): The percentile threshold (default 80).
        shrinkage_factor (float): Factor determining how much to shrink (0-1, default 0.5).

    Returns:
        pd.DataFrame: A new DataFrame with the specified column adjusted and normalized.
    """
    # Create a copy of the DataFrame
    df_copy = df.copy()

    # Extract the column to process
    data = df_copy[column_name]

    # Calculate the percentile threshold using only non-zero values
    non_zero_data = data[data != 0]
    percentile_value = np.percentile(non_zero_data, percentile)

    # Function to shrink values above the threshold
    def shrink_high_value(x):
        if x > percentile_value:
            # Shrink towards the percentile value
            return x + (percentile_value - x) * shrinkage_factor
        return x

    # Apply shrinkage only to values above the threshold
    adjusted_data = data.copy()
    # Only process non-zero values
    adjusted_data[data != 0] = non_zero_data.apply(shrink_high_value)

    # Normalize to sum to 1
    total_sum = adjusted_data.sum()
    if total_sum > 0:
        adjusted_data = adjusted_data / total_sum

    # Update the DataFrame copy
    df_copy[column_name] = adjusted_data

    return df_copy