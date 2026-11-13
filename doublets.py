import numpy as np
import scipy.sparse as sp
import joblib
import xgboost as xgb
from importlib import resources
from pathlib import Path

def get_doublet_path(
    version: str,
    species: str,
    modality: str,
) -> tuple[Path, Path, Path]:
    """
    Resolve packaged doublet model paths inside epitome_tools.
    
    Returns:
        model_path: Path to the XGBoost JSON model
        label_encoder_path: Path to the label encoder (joblib)
        threshold_path: Path to the threshold (joblib)
    """
    if modality == "rna":
        task = "doublet_rna"
        model_file = "rna_model_binary.json"
        label_encoder_file = "rna_label_encoder_binary.pkl"
        threshold_file = "final_threshold.pkl"
    elif modality == "atac":
        task = "doublet_atac"
        model_file = "atac_model_binary.json"
        label_encoder_file = "atac_label_encoder_binary.pkl"
        threshold_file = "atac_final_threshold.pkl"
    else:
        raise ValueError(f"Unsupported modality: {modality}")

    base_path = resources.files("epitome_tools") / "models" / version / species / task

    model_path = base_path / model_file
    label_encoder_path = base_path / label_encoder_file
    threshold_path = base_path / threshold_file

    # Optional sanity check
    for path in (model_path, label_encoder_path, threshold_path):
        if not path.exists():
            raise FileNotFoundError(f"Packaged file not found: {path}")

    return model_path, label_encoder_path, threshold_path


def load_doublet_model(version="v_0.02",species="mouse", modality="rna"):
    """
    Load the XGBoost model for doublet prediction from the specified path.
    """

    model_path, label_encoder_path, threshold_path = get_doublet_path(version, species, modality)
    model = xgb.XGBClassifier()
    model.load_model(model_path)

    # Get the feature names (these should match the features used for training)
    feature_names = model.feature_names_in_

    label_encoder = joblib.load(label_encoder_path)
    
    threshold = joblib.load(threshold_path)

    return model, label_encoder, threshold, feature_names


def prepare_matrix_doublet(adata, feature_names,  active_assay="sc", nan_or_zero='zero', verbose=False):
    """
    Build the model-ready feature matrix (genes + assay one-hot + pct_counts_kept)
    aligned to the order expected by the doublet model, filling missing features
    with NaN or 0 depending on `nan_or_zero`.
    """

    assay_features1 = [f for f in feature_names if f.startswith('pct_counts_kept')]
    assay_features2 = [f for f in feature_names if f.endswith('one_hot')]

    assay_features = assay_features1 + assay_features2

    #number of cells
    n_obs = adata.n_obs

    #active assay must be one of sc, sn, multi_rna, if neither, set it as sc
    if active_assay not in ['sc', 'sn', 'multi_rna']:
        print(f"Warning: Active assay '{active_assay}' not recognized. Defaulting to 'sc'.")
        active_assay = 'sc'
        #add as prefix
    active_assay = f"{active_assay}_one_hot"

    if nan_or_zero == 'nan':
        X_final = np.full((n_obs, len(feature_names)), np.nan, dtype=np.float32)
    elif nan_or_zero == 'zero':
        X_final = np.zeros((n_obs, len(feature_names)), dtype=np.float32)

    target_feature_indices = {name: i for i, name in enumerate(feature_names)}

    # Assay columns are always filled, defaulting to 0
    for f in assay_features:
        X_final[:, target_feature_indices[f]] = 1.0 if f == active_assay else 0.0

    #add adata.obs['pct_counts_kept'] if exists
    if 'pct_counts_kept' in adata.obs.columns and 'pct_counts_kept' in target_feature_indices:
        if verbose:
            print("Adding 'pct_counts_kept' to assay data.")
            print(adata.obs['pct_counts_kept'].values[:5])
        X_final[:, target_feature_indices['pct_counts_kept']] = adata.obs['pct_counts_kept'].values.astype(np.float32)

    # Mapping from gene name to its column index in adata
    gene_feature_indices = {name: i for i, name in enumerate(adata.var_names)}
    assay_feature_set = set(assay_features)

    available_features = [f for f in feature_names if f in gene_feature_indices and f not in assay_feature_set]
    missing_features = [f for f in feature_names if f not in gene_feature_indices and f not in assay_feature_set]

    if missing_features:
        print(f"Warning: {len(missing_features)} features required by model 2 are missing: {missing_features[:5]}...")

    if verbose:
        print(f"Found {len(available_features) + len(assay_features)} available features out of {len(feature_names)} required for model 2.")
    source_indices = [gene_feature_indices[f] for f in available_features]

    target_indices = [target_feature_indices[f] for f in available_features]

    # Only the model's columns are sliced out and densified, so adata.X is never copied in full
    adata_X = adata.X
    if sp.issparse(adata_X):
        X_final[:, target_indices] = adata_X[:, source_indices].toarray()
        if verbose:
            print("Filled final matrix for model 2 from sparse data.")
    else:
        X_final[:, target_indices] = adata_X[:, source_indices]
        if verbose:
            print("Filled final matrix for model 2 from dense data.")
    return X_final

def perform_doublet_prediction(matrix, model, label_encoder, threshold, verbose=False):
    """
    Predict doublets, apply the calibrated threshold, and derive the final
    human-readable label and per-cell uncertainty.

    Returns:
        predicted_doublet_labels: Raw model predictions
        is_doublet: Boolean array, True = real/singlet, False = doublet (thresholded)
        doublet_score: Confidence score
        final_doublet_label: 'real' or 'doublet', derived from `is_doublet`
        doublet_uncertainty: 1 - max(P(doublet), P(real)) per cell
    """
    n_obs = matrix.shape[0]
    probas = model.predict_proba(matrix) # Get probabilities
    predicted_labels = probas.argmax(axis=1)
    predicted_doublet_labels = label_encoder.inverse_transform(predicted_labels)
    if verbose:
        print("Doublet prediction with model complete.")

    if probas.shape[1] > 1:
        doublet_probabilities = probas[:, 1]
    elif probas.shape[1] == 1:
        doublet_probabilities = probas[:, 0]
    else:
        doublet_probabilities = np.zeros(n_obs)
        print("Warning: Model 2 probas has shape < 1")

    doublet_score = 1 - doublet_probabilities
    # Defaults to real/singlet (True) unless the score crosses the calibrated threshold
    is_doublet = doublet_probabilities <= (1 - threshold)

    final_doublet_label = np.where(is_doublet, "real", "doublet")
    doublet_uncertainty = np.minimum(doublet_score, 1 - doublet_score)

    return predicted_doublet_labels, is_doublet, doublet_score, final_doublet_label, doublet_uncertainty