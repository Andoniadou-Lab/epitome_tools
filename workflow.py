import numpy as np
import pandas as pd
import scipy.sparse as sp
import joblib
from pathlib import Path
import xgboost as xgb
from .celltyping import load_celltype_model, prepare_matrix_celltype, perform_celltype_prediction,smoothing_cell_types
from .doublets import load_doublet_model, prepare_matrix_doublet, perform_doublet_prediction

def check_sample_compatibility_features(adata, feature_names, return_present=True,
                                        return_missing=True):
    varnames = adata.var_names

    #check percentage of features that are present in the adata object
    # hash-based lookup; np.isin on string arrays is very slow with ATAC-sized var_names
    present = pd.Index(feature_names).isin(varnames)
    present_percentage = np.sum(present) / len(feature_names) * 100
    print(f"Percentage of features present in adata: {present_percentage:.2f}%")
    if return_present:
        present_features = feature_names[present]
        print(f"Number of features present in adata: {len(present_features)}")
    else:
        present_features = None
    if return_missing:
        missing_features = feature_names[~present]
        print(f"Number of features missing in adata: {len(missing_features)}")
    else:
        missing_features = None

    #pass if at least 70% of the features are present
    if present_percentage < 70:
        print("Warning: Less than 70% of the features are present in the adata object.")
        passing = False
    else:
        print("At least 70% of the features are present in the adata object.")
        passing = True

    
    if return_present and return_missing:
        return passing, present_features, missing_features
    
    elif return_present:
        return passing, present_features
    elif return_missing:
        return passing, missing_features
    else:
        return passing
    

def check_sample_compatibility_normalization(adata, force=False):
    # Check if the model is for RNA or ATAC
    #check number of cells, if more than 50000, throw an error
    if adata.shape[0] > 50000:
        print(f"Error: The number of cells in the dataset is greater than 50000. Are you sure this has been filtered correctly?")
    
    #take info from first 10 cells
    first_20_cells = adata[:20, :].X
    #check if they have integer values or if 1.0 occurs more than 10% of the time
    integer_type = np.issubdtype(first_20_cells.dtype, np.integer)
    if integer_type:
        print(f"Warning: The dataset appears to be in integer format. Are you sure this has been normalized correctly?")
    
    not_normed = False
    one_values = np.sum(first_20_cells == 1.0) / first_20_cells.size
    if one_values > 0.05:
        not_normed = True
        print(f"Warning: The dataset has more than 5% of the values equal to 1.0. Are you sure this has been normalized correctly?")

    one_values = np.sum(first_20_cells == 1) / first_20_cells.size
    if one_values > 0.05:
        not_normed = True
        print(f"Warning: The dataset has more than 5% of the values equal to 1.0. Are you sure this has been normalized correctly?")

    #if all cells sum to nearly 10k, say it hasnt been logged
    not_logged = False
    if np.all(np.sum(first_20_cells, axis=1) > 9500) and np.all(np.sum(first_20_cells, axis=1) < 10500):
        print(f"Warning: Have you logged the dataset? The cells sum to nearly 10k.")
        not_logged = True
    
    #if force is True, return True
    passing = False
    if force:
        print(f"Warning: Force is set to True. Passing the dataset compatibility check.")
        passing = True
    elif not integer_type and one_values < 0.1 and not not_logged:
        print(f"The dataset has passed the compatibility check.")
        passing = True

    return passing, not_normed, not_logged


def pct_counts_kept(adata):
  counts_per_10k_current = adata.X.sum(axis=1) / 10000

  return np.asarray(counts_per_10k_current).flatten()


def calc_pct_counts_kept(adata, features):
    feature_indices = {name: i for i, name in enumerate(adata.var_names)}
    source_indices = [feature_indices[f] for f in features if f in feature_indices]

    # Slice only the model features so adata.X is never copied in full
    adata_X = adata.X
    if sp.issparse(adata_X):
        #undo log1p
        kept_X = adata_X[:, source_indices].expm1()
    else:
        kept_X = np.expm1(np.ascontiguousarray(adata_X[:, source_indices]))

    adata.obs["pct_counts_kept"] = np.asarray(kept_X.sum(axis=1) / 10000).flatten()
    return adata




def cell_type_workflow(adata_to_use, active_assay="sc",modality="rna",in_place=True, nan_or_zero='zero',smoothing=True, version="v_0.02",species="mouse", verbose=False):
    """
    Main workflow for cell type prediction. Adds `predicted_cell_type`,
    `predicted_cell_type_proba`, `cell_type_uncertainty`, per-class `proba_*`
    columns, and (if `smoothing`) `cell_type_final` to `adata.obs`.
    """

    if species=="human":
        print("Human model is still in development and is not recommended for use, results may not be accurate.")

    if in_place:
        adata = adata_to_use
    else:
        adata = adata_to_use.copy()
    
    adata.var_names_make_unique()

    model, label_encoder, feature_names = load_celltype_model(version=version, species=species, modality=modality)

    #checks
    check_sample_compatibility_features(adata, feature_names, return_present=False, return_missing=False)
    check_sample_compatibility_normalization(adata, force=False)

    # Prepare the matrix for cell type prediction
    X_final = prepare_matrix_celltype(adata, feature_names, active_assay=active_assay, nan_or_zero=nan_or_zero, verbose=verbose)
    # Perform cell type prediction
    predicted_cell_types, probas = perform_celltype_prediction(X_final, model, label_encoder)
    # Add predictions to adata
    
    adata.obs['predicted_cell_type'] = predicted_cell_types
    adata.obs['predicted_cell_type_proba'] = probas.max(axis=1)  # Store max probability
    adata.obs['cell_type_uncertainty'] = 1 - probas.max(axis=1)

    for i, label in enumerate(label_encoder.classes_):
        adata.obs[f'proba_{label}'] = probas[:, i]


    if smoothing:
        adata = smoothing_cell_types(adata)
    return adata


def doublet_workflow(adata_to_use,active_assay="sc",modality="rna",in_place=True, nan_or_zero='zero', version="v_0.02",species="mouse", verbose=False):
    """
    Main workflow for doublet prediction. Adds `init_predicted_doublet_epitome`,
    `thresholded_doublet_epitome`, `doublet_score_epitome`, `final_doublet_label`,
    and `doublet_uncertainty` to `adata.obs`.
    """

    if species=="human":
        print("Human model is still in development and is not recommended for use, results may not be accurate.")


    if in_place:
        adata = adata_to_use
    else:
        adata = adata_to_use.copy()
    adata.var_names_make_unique()


    model, label_encoder, threshold, feature_names = load_doublet_model(version=version, species=species, modality=modality)


    #checks
    check_sample_compatibility_features(adata, feature_names, return_present=False, return_missing=False)
    check_sample_compatibility_normalization(adata, force=False)

    adata = calc_pct_counts_kept(adata, feature_names)
    # Prepare the matrix for doublet prediction
    X_final = prepare_matrix_doublet(adata, feature_names, active_assay=active_assay, nan_or_zero=nan_or_zero, verbose=verbose)
    # Perform doublet prediction
    predicted_doublet_labels, is_doublet, doublet_score, final_doublet_label, doublet_uncertainty = perform_doublet_prediction(X_final, model, label_encoder, threshold, verbose=verbose)
    # Add predictions to adata
    
    adata.obs['init_predicted_doublet_epitome'] = predicted_doublet_labels
    adata.obs['thresholded_doublet_epitome'] = is_doublet
    adata.obs['doublet_score_epitome'] = doublet_score
    adata.obs['final_doublet_label'] = final_doublet_label
    adata.obs['doublet_uncertainty'] = doublet_uncertainty
    return adata

def celltype_doublet_workflow(adata, active_assay="sc", modality="rna", in_place=True, nan_or_zero='zero', smoothing=True, version="v_0.02",species="mouse", verbose=False):
    """
    Main workflow for cell type and doublet prediction.
    """
    if species=="human":
        print("Human model is still in development and is not recommended for use, results may not be accurate.")
        
    adata = cell_type_workflow(adata, active_assay=active_assay, modality=modality, in_place=in_place, nan_or_zero=nan_or_zero, smoothing=smoothing, version=version, species=species, verbose=verbose)
    adata = doublet_workflow(adata, active_assay=active_assay, modality=modality, in_place=in_place, nan_or_zero=nan_or_zero, version=version, species=species, verbose=verbose)
    return adata