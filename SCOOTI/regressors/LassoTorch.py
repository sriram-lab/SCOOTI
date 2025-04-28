import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, Subset
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from tqdm import tqdm

class LassoDataset(Dataset):
    def __init__(self, X: pd.DataFrame, y_col: np.ndarray):
        self.X = torch.tensor(X.values, dtype=torch.float32)
        self.y = torch.tensor(y_col.reshape(-1, 1), dtype=torch.float32)

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]


class LassoRegressor:
    def __init__(self, l1_lambda=1e-3, lr=1e-3, epochs=300, batch_size=1000, verbose=True):
        self.l1_lambda = l1_lambda
        self.lr = lr
        self.epochs = epochs
        self.batch_size = batch_size
        self.verbose = verbose
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.feature_names = None
        self.target_names = None
        self.coef_df = None

    class LinearModel(nn.Module):
        def __init__(self, input_dim):
            super().__init__()
            self.linear = nn.Linear(input_dim, 1, bias=False)

        def forward(self, x):
            return self.linear(x)

    def _train_one_fold(self, X, y_col, train_idx):
        dataset = LassoDataset(X, y_col)
        loader = DataLoader(Subset(dataset, train_idx), batch_size=self.batch_size, shuffle=True)
        model = self.LinearModel(X.shape[1]).to(self.device)
        optimizer = optim.Adam(model.parameters(), lr=self.lr)
        loss_fn = nn.MSELoss()

        for epoch in range(self.epochs):
            model.train()
            for X_batch, y_batch in loader:
                X_batch, y_batch = X_batch.to(self.device), y_batch.to(self.device)
                optimizer.zero_grad()
                preds = model(X_batch)
                mse_loss = loss_fn(preds, y_batch)
                l1_penalty = sum(p.abs().sum() for p in model.parameters())
                loss = mse_loss + self.l1_lambda * l1_penalty
                loss.backward()
                optimizer.step()

                # 🔒 Clamp weights to be non-negative
                with torch.no_grad():
                    model.linear.weight.data.clamp_(min=0.0)

        with torch.no_grad():
            weights = model.linear.weight.view(-1).cpu().numpy()
        return weights

    def fit_cv(self, X: pd.DataFrame, y: pd.DataFrame, k=5):
        """Fit K models on K folds for each output column and average coefficients"""
        self.feature_names = X.columns.tolist()
        self.target_names = y.columns.tolist()
        n_targets = y.shape[1]
        n_features = X.shape[1]

        all_coefs = []

        for j in tqdm(range(n_targets)):
            #if self.verbose and j % 100 == 0:
            #    print(f"Fitting target {j+1}/{n_targets}")
            y_col = y.iloc[:, j].values
            kf = KFold(n_splits=k, shuffle=True, random_state=42)
            fold_coefs = []

            for train_idx, _ in kf.split(X):
                weights = self._train_one_fold(X, y_col, train_idx)
                fold_coefs.append(weights)

            avg_weights = np.mean(np.stack(fold_coefs), axis=0)
            all_coefs.append(avg_weights)

        coef_matrix = np.stack(all_coefs, axis=1)  # shape: n_features x n_targets
        self.coef_df = pd.DataFrame(coef_matrix, index=self.feature_names, columns=self.target_names)

    def get_coefficients(self) -> pd.DataFrame:
        return self.coef_df

