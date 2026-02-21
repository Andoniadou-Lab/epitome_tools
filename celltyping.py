import xgboost as xgb
import joblib
from pathlib import Path
import scipy.sparse as sp
import numpy as np
import scanpy as sc
import pandas as pd
from collections import Counter

from importlib import resources
from pathlib import Path

def get_cell_typing_path(
    version: str,
    species: str,
    modality: str,
) -> tuple[Path, Path]:
    """
    Resolve packaged doublet model paths inside epitome_tools.
    
    Returns:
        model_path: Path to the XGBoost JSON model
        label_encoder_path: Path to the label encoder (joblib)
    """
    if modality == "rna":
        task = "cell_typing_rna"
        model_file = "rna_model.json"
        label_encoder_file = "label_encoder_rna.pkl"
    elif modality == "atac":
        task = "cell_typing_atac"
        model_file = "atac_model.json"
        label_encoder_file = "label_encoder_atac.pkl"
    else:
        raise ValueError(f"Unsupported modality: {modality}")

    base_path = resources.files("epitome_tools") / "models" / version / species / task

    model_path = base_path / model_file
    label_encoder_path = base_path / label_encoder_file

    # Optional sanity check
    for path in (model_path, label_encoder_path):
        if not path.exists():
            raise FileNotFoundError(f"Packaged file not found: {path}")

    return model_path, label_encoder_path


def load_celltype_model(version="v_0.02",species="mouse", modality="rna"):
    """
    Load the XGBoost model from the specified path.
    """
    model_path, label_encoder_path = get_cell_typing_path(version, species, modality)
    model = xgb.XGBClassifier()
    model.load_model(model_path)

    # Access the booster and retrieve the feature names
    booster = model.get_booster()

    # Get the feature names (these should match the features used for training)
    feature_names = model.feature_names_in_

    label_encoder = joblib.load(label_encoder_path)

    return model, label_encoder, feature_names



def prepare_matrix_celltype(adata, feature_names,active_assay="sc", nan_or_zero='nan'):
    """
    Prepare the AnnData object for prediction by selecting the relevant features.
    """

    adata_to_handle = adata.copy()

    # Define all potential assay features based on feature_names_sorted1
    assay_features = [f for f in feature_names if f.startswith('assay_')]

    #number of cells
    n_obs = adata_to_handle.n_obs

    #initialise assay data with zeros
    assay_data1 = np.zeros((n_obs, len(assay_features)), dtype=np.float32)

    # Create a mapping from assay feature name to its column index
    assay_feature_indices1 = {name: i for i, name in enumerate(assay_features)}

    #active assay must be one of sc, sn, multi_rna, if neither, set it as sc
    if active_assay not in ['sc', 'sn', 'multi_rna']:
        print(f"Warning: Active assay '{active_assay}' not recognized. Defaulting to 'sc'.")
        active_assay = 'sc'
    #add prefix
    active_assay = f"assay_{active_assay}_onehot"

    if active_assay in assay_feature_indices1:
        assay_data1[:, assay_feature_indices1[active_assay]] = 1.0

    # Convert assay data to sparse if original data is sparse
    assay_data_matrix1 = sp.csr_matrix(assay_data1) if sp.issparse(adata_to_handle.X) else assay_data1

    # --- 2. Combine Gene Expression and Assay Features for Model 1 ---
    # We need to combine these temporarily to easily subset later

    # Ensure original data is CSR for efficient column slicing if sparse
    if sp.issparse(adata_to_handle.X) and not isinstance(adata_to_handle.X, sp.csr_matrix):
        adata_X = adata_to_handle.X.tocsr().copy()
        print("Converted adata_orig.X to CSR format.")
    else:
        adata_X = adata_to_handle.X.copy()

    # Combine the matrices horizontally
    combined_X1 = sp.hstack([adata_X, assay_data_matrix1], format='csr') if sp.issparse(adata_X) else np.hstack([adata_X, assay_data1])

    # Create combined feature names list
    combined_feature_names1 = adata_to_handle.var_names.tolist() + assay_features

    # Create a mapping from the combined feature names to their column index
    combined_feature_indices1 = {name: i for i, name in enumerate(combined_feature_names1)}
    print(f"Combined matrix shape for model 1: {combined_X1.shape}")

    # Initialize the final matrix with NaNs (XGBoost can handle NaNs)
    if nan_or_zero == 'nan':
        X_final1 = np.full((n_obs, len(feature_names)), np.nan, dtype=np.float32)
    elif nan_or_zero == 'zero':
        X_final1 = np.zeros((n_obs, len(feature_names)), dtype=np.float32)

    # Create a mapping for the target feature order
    target_feature_indices1 = {name: i for i, name in enumerate(feature_names)}

    # Find which features required by the model are present in our combined data
    available_features1 = [f for f in feature_names if f in combined_feature_indices1]
    missing_features1 = [f for f in feature_names if f not in combined_feature_indices1]

    if missing_features1:
        print(f"Warning: {len(missing_features1)} features required by model 1 are missing from the data: {missing_features1[:5]}...") # Print first 5

    print(f"Found {len(available_features1)} available features out of {len(feature_names)} required for model 1.")

    # Get the column indices in the *combined* data for the available features - indices of where it is in the combined matrix
    source_indices1 = [combined_feature_indices1[f] for f in available_features1]

    # Get the column indices in the *final* matrix for these available features - indices of where the model expects it
    target_indices1 = [target_feature_indices1[f] for f in available_features1]

    # Fill the final matrix with data from the available features
    # Ensure data is dense for assignment; handle potential memory issues for large datasets
    if sp.issparse(combined_X1):
        # Slice sparse matrix efficiently and convert to dense for assignment
        X_final1[:, target_indices1] = combined_X1[:, source_indices1].toarray()
        print("Filled final matrix for model 1 from sparse data.")
    else:
        X_final1[:, target_indices1] = combined_X1[:, source_indices1]
        print("Filled final matrix for model 1 from dense data.")
    
    
    return X_final1


def perform_celltype_prediction(matrix, model, label_encoder, return_probas=True):
    """
    Perform cell type prediction using the provided model and label encoder.
    """
    probas = model.predict_proba(matrix)
    predicted_labels = model.predict(matrix)
    predicted_cell_types = label_encoder.inverse_transform(predicted_labels)

    if return_probas:
        return predicted_cell_types, probas
    else:
        return predicted_cell_types
    

def smoothing_cell_types(adata):
    """
    Checks for PCA embedding, generates if missing, finds 10 closest neighbors,
    and reassigns cell identity based on neighbor majority.
    """
    if 'X_pca' not in adata.obsm:
        sc.tl.pca(adata)

    sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_pca')

    new_cell_types = adata.obs['predicted_cell_type'].copy()

    for i in range(adata.n_obs):
        # Access neighbor indices from the connectivities matrix in adata.obsp
        neighbor_indices = adata.obsp['connectivities'].indices[
            adata.obsp['connectivities'].indptr[i]:
            adata.obsp['connectivities'].indptr[i+1]
        ]
        neighbor_cell_types = adata.obs['predicted_cell_type'].iloc[neighbor_indices]
        type_counts = Counter(neighbor_cell_types)

        for cell_type, count in type_counts.items():
            if count > 5:
                new_cell_types.iloc[i] = cell_type
                break

    adata.obs['cell_type_final'] = new_cell_types
    return adata