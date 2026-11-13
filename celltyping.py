import xgboost as xgb
import joblib
from pathlib import Path
import scipy.sparse as sp
import numpy as np
import scanpy as sc
import pandas as pd
from importlib import resources

def get_cell_typing_path(
    version: str,
    species: str,
    modality: str,
) -> tuple[Path, Path]:
    """
    Resolve packaged cell typing model paths inside epitome_tools.
    
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

    # Get the feature names (these should match the features used for training)
    feature_names = model.feature_names_in_

    label_encoder = joblib.load(label_encoder_path)

    return model, label_encoder, feature_names



def prepare_matrix_celltype(adata, feature_names,active_assay="sc", nan_or_zero='zero', verbose=False):
    """
    Build the model-ready feature matrix (genes + assay one-hot) aligned to
    the order expected by the cell typing model, filling missing features
    with NaN or 0 depending on `nan_or_zero`.
    """

    # Define all potential assay features based on feature_names_sorted1
    assay_features = [f for f in feature_names if f.startswith('assay_')]

    #number of cells
    n_obs = adata.n_obs

    #active assay must be one of sc, sn, multi_rna, if neither, set it as sc
    if active_assay not in ['sc', 'sn', 'multi_rna']:
        print(f"Warning: Active assay '{active_assay}' not recognized. Defaulting to 'sc'.")
        active_assay = 'sc'
    #add prefix
    active_assay = f"assay_{active_assay}_onehot"

    # Initialize the final matrix with NaNs (XGBoost can handle NaNs)
    if nan_or_zero == 'nan':
        X_final1 = np.full((n_obs, len(feature_names)), np.nan, dtype=np.float32)
    elif nan_or_zero == 'zero':
        X_final1 = np.zeros((n_obs, len(feature_names)), dtype=np.float32)

    # Create a mapping for the target feature order
    target_feature_indices1 = {name: i for i, name in enumerate(feature_names)}

    # Assay one-hot columns are always filled: 1 for the active assay, 0 otherwise
    for f in assay_features:
        X_final1[:, target_feature_indices1[f]] = 1.0 if f == active_assay else 0.0

    # Mapping from gene name to its column index in adata
    gene_feature_indices1 = {name: i for i, name in enumerate(adata.var_names)}
    assay_feature_set1 = set(assay_features)

    # Find which features required by the model are present in our data
    available_features1 = [f for f in feature_names if f in gene_feature_indices1 and f not in assay_feature_set1]
    missing_features1 = [f for f in feature_names if f not in gene_feature_indices1 and f not in assay_feature_set1]

    if missing_features1:
        print(f"Warning: {len(missing_features1)} features required by model 1 are missing from the data: {missing_features1[:5]}...") # Print first 5

    if verbose:
        print(f"Found {len(available_features1) + len(assay_features)} available features out of {len(feature_names)} required for model 1.")

    # Get the column indices in adata for the available features
    source_indices1 = [gene_feature_indices1[f] for f in available_features1]

    # Get the column indices in the *final* matrix for these available features - indices of where the model expects it
    target_indices1 = [target_feature_indices1[f] for f in available_features1]

    # Only the model's columns are sliced out and densified, so adata.X is never copied in full
    adata_X = adata.X
    if sp.issparse(adata_X):
        X_final1[:, target_indices1] = adata_X[:, source_indices1].toarray()
        if verbose:
            print("Filled final matrix for model 1 from sparse data.")
    else:
        X_final1[:, target_indices1] = adata_X[:, source_indices1]
        if verbose:
            print("Filled final matrix for model 1 from dense data.")

    return X_final1


def perform_celltype_prediction(matrix, model, label_encoder, return_probas=True):
    """
    Perform cell type prediction using the provided model and label encoder.
    """
    probas = model.predict_proba(matrix)
    predicted_labels = probas.argmax(axis=1)
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

    connectivities = adata.obsp['connectivities']
    cell_type_codes = pd.Categorical(adata.obs['predicted_cell_type'])
    n_types = len(cell_type_codes.categories)

    # One entry per stored neighbor, kept in connectivities storage order
    neighbor_rows = np.repeat(np.arange(adata.n_obs), np.diff(connectivities.indptr))
    neighbor_codes = cell_type_codes.codes[connectivities.indices]
    has_type = neighbor_codes >= 0
    neighbor_rows, neighbor_codes = neighbor_rows[has_type], neighbor_codes[has_type]

    type_counts = np.zeros((adata.n_obs, n_types), dtype=np.int64)
    np.add.at(type_counts, (neighbor_rows, neighbor_codes), 1)

    # A cell takes the first type, in neighbor order, seen more than 5 times among its neighbors
    qualifies = type_counts[neighbor_rows, neighbor_codes] > 5
    qualifying_rows, first_qualifying = np.unique(neighbor_rows[qualifies], return_index=True)
    new_cell_types.iloc[qualifying_rows] = np.asarray(cell_type_codes.categories)[
        neighbor_codes[qualifies][first_qualifying]
    ]

    adata.obs['cell_type_final'] = new_cell_types
    return adata