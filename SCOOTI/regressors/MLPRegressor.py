import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import KFold
import pandas as pd
import numpy as np
from tqdm import tqdm

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import KFold
import pandas as pd
import numpy as np

class GOFluxDataset(Dataset):
    def __init__(self, X: pd.DataFrame, y_col: np.ndarray):
        self.X = torch.tensor(X.values, dtype=torch.float32)
        self.y = torch.tensor(y_col.reshape(-1, 1), dtype=torch.float32)

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

class BetaRegressor(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.fc1 = nn.Linear(input_dim, 128)
        self.fc_alpha = nn.Linear(128, input_dim)
        self.fc_beta = nn.Linear(128, input_dim)

    def forward(self, X):
        h = F.relu(self.fc1(X))
        alpha = F.softplus(self.fc_alpha(h)) + 1e-6  # Ensure positivity
        beta_param = F.softplus(self.fc_beta(h)) + 1e-6
        gamma_dist = torch.distributions.Gamma(alpha, beta_param)
        beta = gamma_dist.rsample()
        return beta, alpha, beta_param

def train_model(model, dataloader, optimizer, epochs=100, l1_weight=1e-4, entropy_weight=1e-4, kl_weight=1e-3):
    model.train()
    for epoch in range(epochs):
        for X_batch, y_batch in dataloader:
            optimizer.zero_grad()
            beta, alpha, beta_param = model(X_batch)
            y_hat = torch.sum(X_batch * beta, dim=1, keepdim=True)  # predicted flux: (batch_size, 1)
            mse_loss = F.mse_loss(y_hat, y_batch)
            l1_loss = beta.abs().mean()
            entropy = - (F.softmax(beta, dim=1) * F.log_softmax(beta, dim=1)).sum(dim=1).mean()
            kl_loss = torch.mean(alpha - torch.log(beta_param) + torch.lgamma(alpha) + (1 - alpha) * torch.digamma(alpha))
            loss = mse_loss + l1_weight * l1_loss + entropy_weight * entropy + kl_weight * kl_loss
            loss.backward()
            optimizer.step()

def infer_beta(model, X_tensor):
    model.eval()
    with torch.no_grad():
        beta, _, _ = model(X_tensor)
        # In Gamma model, third return is beta_param
        return beta.mean(dim=0).cpu().numpy()

class TrainValWrapper:
    def __init__(self, epochs=50, batch_size=32, patience=10, verbose=True):
        self.epochs = epochs
        self.batch_size = batch_size
        self.patience = patience
        self.verbose = verbose
        self.coef_df = None

    def train_model(self, model, train_loader, val_loader, optimizer):
        best_loss = float('inf')
        best_model_state = None
        counter = 0

        for epoch in range(self.epochs):
            model.train()
            for X_batch, y_batch in train_loader:
                optimizer.zero_grad()
                beta, alpha, beta_param = model(X_batch)
                y_hat = torch.sum(X_batch * beta, dim=1, keepdim=True)
                mse_loss = F.mse_loss(y_hat, y_batch)
                l1_loss = beta.abs().mean()
                entropy = -(F.softmax(beta, dim=1) * F.log_softmax(beta, dim=1)).sum(dim=1).mean()
                kl_loss = torch.mean(alpha - torch.log(beta_param) + torch.lgamma(alpha) + (1 - alpha) * torch.digamma(alpha))
                loss = mse_loss + 1e-4 * l1_loss + 1e-4 * entropy + 1e-3 * kl_loss
                loss.backward()
                optimizer.step()

            # Validation
            model.eval()
            val_losses = []
            with torch.no_grad():
                for X_val, y_val in val_loader:
                    beta, alpha, beta_param = model(X_val)
                    y_hat = torch.sum(X_val * beta, dim=1, keepdim=True)
                    val_loss = F.mse_loss(y_hat, y_val)
                    val_losses.append(val_loss.item())
            avg_val_loss = np.mean(val_losses)

            if self.verbose:
                print(f"Epoch {epoch+1}: Validation Loss = {avg_val_loss:.6f}")

            if avg_val_loss < best_loss - 1e-4:
                best_loss = avg_val_loss
                best_model_state = {k: v.clone() for k, v in model.state_dict().items()}
                counter = 0
            else:
                counter += 1
                if counter >= self.patience:
                    print(f"Final Epoch {epoch+1}: Validation Loss = {avg_val_loss:.6f}")
                    print("Early stopping triggered!")
                    break

        if best_model_state:
            model.load_state_dict(best_model_state)

    def run_train_val(self, X: pd.DataFrame, y: pd.DataFrame):
        from sklearn.model_selection import train_test_split

        self.feature_names = X.columns.tolist()
        self.target_names = y.columns.tolist()
        all_coefs = []

        for j in tqdm(range(y.shape[1])):
            #if self.verbose and j % 100 == 0:
            #    print(f"Fitting target {j+1}/{y.shape[1]}")
            y_col = y.iloc[:, j].values
            X_train_df, X_val_df, y_train, y_val = train_test_split(X, y_col, test_size=0.2, random_state=42)
            train_loader = DataLoader(GOFluxDataset(X_train_df, y_train), batch_size=self.batch_size, shuffle=True)
            val_loader = DataLoader(GOFluxDataset(X_val_df, y_val), batch_size=self.batch_size, shuffle=False)
            model = BetaRegressor(input_dim=X.shape[1])
            optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
            self.train_model(model, train_loader, val_loader, optimizer)
            X_tensor = torch.tensor(X.values, dtype=torch.float32)
            coef = infer_beta(model, X_tensor)
            all_coefs.append(coef)

        coef_matrix = np.stack(all_coefs, axis=1)
        self.coef_df = pd.DataFrame(coef_matrix, index=self.feature_names, columns=self.target_names)
        return self.coef_df

