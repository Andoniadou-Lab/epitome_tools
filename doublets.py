import numpy as np
import scipy.sparse as sp
import joblib
import xgboost as xgb
from pathlib import Path


def load_doublet_model(model_path,label_encoder_path,threshold_path):
    """
    Load the XGBoost model for doublet prediction from the specified path.
    """
    model = xgb.XGBClassifier()
    model.load_model(model_path)

    # Access the booster and retrieve the feature names
    booster = model.get_booster()

    # Get the feature names (these should match the features used for training)
    feature_names = model.feature_names_in_

    label_encoder = joblib.load(label_encoder_path)
    
    threshold = joblib.load(threshold_path)

    return model, label_encoder, threshold, feature_names


def prepare_matrix_doublet(adata, feature_names):

    adata_to_handle = adata.copy()

    # Define all potential assay features based on feature_names_sorted1
    assay_features = [f for f in feature_names if f.startswith('total_')]

    #number of cells
    n_obs = adata_to_handle.n_obs

    #initialise assay data with zeros
    assay_data = np.zeros((n_obs, len(assay_features)), dtype=np.float32)

    # Create a mapping from assay feature name to its column index
    assay_feature_indices = {name: i for i, name in enumerate(assay_features)}

    # Convert assay data to sparse if original data is sparse
    assay_data_matrix = sp.csr_matrix(assay_data) if sp.issparse(adata_to_handle.X) else assay_data

    # --- 2. Combine Gene Expression and Assay Features for Model 1 ---
    # We need to combine these temporarily to easily subset later

    # Ensure original data is CSR for efficient column slicing if sparse
    if sp.issparse(adata_to_handle.X) and not isinstance(adata_to_handle.X, sp.csr_matrix):
        adata_X = adata_to_handle.X.tocsr().copy()
        print("Converted adata_orig.X to CSR format.")
    else:
        adata_X = adata_to_handle.X.copy()

    # Combine the matrices horizontally
    combined_X1 = sp.hstack([adata_X, assay_data_matrix], format='csr') if sp.issparse(adata_X) else np.hstack([adata_X, assay_data])

    # Create combined feature names list
    combined_feature_names = adata_to_handle.var_names.tolist() + assay_features

    # Create a mapping from the combined feature names to their column index
    combined_feature_indices = {name: i for i, name in enumerate(combined_feature_names)}
    print(f"Combined matrix shape for model 1: {combined_X1.shape}")

    X_final = np.full((n_obs, len(feature_names)), np.nan, dtype=np.float32)
    target_feature_indices = {name: i for i, name in enumerate(feature_names)}

    # Reuse combined_feature_indices1, available_features1, source_indices1 from Model 1
    available_features = [f for f in feature_names if f in combined_feature_indices] 
    missing_features = [f for f in feature_names if f not in combined_feature_indices]

    if missing_features:
        print(f"Warning: {len(missing_features)} features required by model 2 are missing: {missing_features[:5]}...")

    print(f"Found {len(available_features)} available features out of {len(feature_names)} required for model 2.")
    source_indices = [combined_feature_indices[f] for f in available_features] # Use combined_feature_indices1

    target_indices = [target_feature_indices[f] for f in available_features]

    if sp.issparse(combined_X1):
        # Slice sparse matrix efficiently and convert to dense for assignment
        X_final[:, target_indices] = combined_X1[:, source_indices].toarray()
        print("Filled final matrix for model 1 from sparse data.")
    else:
        X_final[:, target_indices] = combined_X1[:, source_indices]
        print("Filled final matrix for model 1 from dense data.")
    return X_final

def perform_doublet_prediction(matrix, model, label_encoder, threshold):

    n_obs = matrix.shape[0]
    predicted_labels = model.predict(matrix)
    predicted_doublet_labels = label_encoder.inverse_transform(predicted_labels)
    probas = model.predict_proba(matrix) # Get probabilities
    print("Doublet prediction with model complete.")
    is_doublet = np.full(n_obs, False, dtype=bool) # Default: not a doublet
    doublet_score = np.zeros(n_obs)

    if probas.shape[1] > 1:
      doublet_probabilities = probas[:, 1]
    elif probas.shape[1] == 1:
       doublet_probabilities = probas[:, 0]
    else:
       doublet_probabilities = np.zeros(n_obs)
       print("Warning: Model 2 probas has shape < 1")

    for i in range(n_obs):
        if predicted_doublet_labels[i] == 'doublet': # Only check if model2 predicts doublet
            doublet_score[i] = doublet_probabilities[i]
            if doublet_probabilities[i] < threshold:
                is_doublet[i] = True


    return predicted_doublet_labels, is_doublet, doublet_score