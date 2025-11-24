import pandas as pd
import numpy as np
from sklearn.preprocessing import QuantileTransformer


def median_impute_flux_by_reaction(flux_df):
    """
    Imputes missing (zero) fluxes in flux_df by the median of nonzero fluxes 
    across GO terms (columns) for each reaction (row).

    Parameters:
    - flux_df: pd.DataFrame, shape (n_reactions, n_GO_terms)

    Returns:
    - imputed_flux_df: pd.DataFrame of the same shape, with missing values imputed
    """
    flux_imputed = flux_df.copy()

    for idx in flux_imputed.index:  # Iterate over reactions (rows)
        row = flux_imputed.loc[idx]
        active = row != 0
        if active.sum() > 0:
            median_value = row[active].median()
            flux_imputed.loc[idx, ~active] = median_value
        else:
            # Optional: If the whole row is zero, keep as is or set NaN
            pass  # or: flux_imputed.loc[idx, :] = np.nan

    return flux_imputed




def normalize_flux_data(
    X,
    method="zscore",              # Options: "zscore", "plog", "quantile"
    nonzero_only=False,           # Whether to normalize only nonzero values
    quantile_target="normal",     # Options: "normal" or "uniform" for quantile normalization
    log_eps=1e-6                  # Small epsilon for log1p to avoid log(0)
):
    X_norm = X.values  # Ensure numpy array

    for j in range(X_norm.shape[1]):  # Normalize per column (reaction flux)
        col = X.values[:, j]

        if method == "zscore":
            if nonzero_only:
                active = col != 0
                mean = col[active].mean() if active.sum() > 0 else 0.0
                std = col[active].std() + 1e-6 if active.sum() > 0 else 1.0
                X_norm[active, j] = (col[active] - mean) / std
            else:
                mean = col.mean()
                std = col.std() + 1e-6
                X_norm[:, j] = (col - mean) / std


        elif method == "plog":
            col_log = np.log1p(col + log_eps)
            if nonzero_only:
                active = col != 0
                mean = col_log[active].mean() if active.sum() > 0 else 0.0
                std = col_log[active].std() + 1e-6 if active.sum() > 0 else 1.0
                X_norm[active, j] = (col_log[active] - mean) / std
                X_norm[~active, j] = 0.0  # KEEP original zero positions!
            else:
                mean = col_log.mean()
                std = col_log.std() + 1e-6
                X_norm[:, j] = (col_log - mean) / std

        elif method == "quantile":
            transformer = QuantileTransformer(
                output_distribution=quantile_target,
                n_quantiles=min(X.shape[0], 1000),  # Safe upper limit for small datasets
                random_state=42
            )
            col_reshaped = col.reshape(-1, 1)
            if nonzero_only:
                active = col != 0
                if active.sum() > 1:  # Need at least two points to apply quantile transform
                    col_active = col[active].reshape(-1, 1)
                    col_active_norm = transformer.fit_transform(col_active).flatten()
                    X_norm[active, j] = col_active_norm
                else:
                    X_norm[active, j] = 0.0  # If only one nonzero, set to zero
            else:
                X_norm[:, j] = transformer.fit_transform(col_reshaped).flatten()
        else:
            raise ValueError(f"Unsupported method: {method}")

    if isinstance(X, pd.DataFrame):
       X_norm = pd.DataFrame(X_norm, columns=X.columns, index=X.index)

    return X_norm
