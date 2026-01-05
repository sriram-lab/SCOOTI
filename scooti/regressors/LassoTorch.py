import torch
import torch.nn as nn
import torch.optim as optim
from tqdm import tqdm
import pandas as pd
import numpy as np

class LassoRegressor:
    def __init__(self, X, alpha=1e-3, lr=1e-2, epochs=1000):
        self.objectives = X.columns
        self.X = torch.tensor(X.values, dtype=torch.float32).cuda()  # shape: (n_reactions, n_objectives)
        self.alpha = alpha
        self.lr = lr
        self.epochs = epochs

    def fit(self, y: pd.DataFrame):
        y_tensor = torch.tensor(y.values, dtype=torch.float32).cuda()  # shape: (n_samples, n_reactions)
        n_samples, n_reactions = y_tensor.shape
        n_objectives = self.X.shape[1]

        beta = nn.Parameter(torch.rand(n_samples, n_objectives, device='cuda'))  # β_i for each sample
        optimizer = optim.Adam([beta], lr=self.lr)

        for epoch in tqdm(range(self.epochs), desc="GPU Lasso"):
            optimizer.zero_grad()

            # Optional: Clamp β to be non-negative
            with torch.no_grad():
                beta.clamp_(min=0.0)

            y_hat = beta @ self.X.T  # (n_samples, n_reactions)
            loss = nn.functional.mse_loss(y_hat, y_tensor)
            l1_penalty = self.alpha * beta.abs().mean()
            total_loss = loss + l1_penalty
            total_loss.backward()
            optimizer.step()


        self.beta = beta.detach().cpu().numpy()
        self.samples = y.index
        self.reactions = y.columns

    def get_beta_df(self):
        return pd.DataFrame(self.beta.T, index=self.objectives, columns=self.samples)

    def predict(self):
        y_hat = self.beta @ self.X.cpu().numpy().T
        return pd.DataFrame(y_hat.T, index=self.reactions, columns=self.samples)

