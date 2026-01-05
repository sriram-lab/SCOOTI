"""
analysis_minimal.py
Lightweight analysis utilities to interpret metabolic objective coefficients and
optionally associated flux solutions, with minimal dependencies.
"""
from __future__ import annotations
import os
from pathlib import Path
from typing import Dict, Tuple, Optional, List, Union

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit
from scipy.stats import ttest_ind, mannwhitneyu

# Optional imports from legacy utilities when available
try:
    from scooti.regressionAnalyzer import read_sampling_objFlux as _legacy_read_sampling_objFlux
    from scooti.regressionAnalyzer import find_pareto as _legacy_find_pareto
    _HAS_LEGACY_PARETO = True
except Exception:
    _HAS_LEGACY_PARETO = False


def load_coefficients(coef_paths: Dict[str, str]) -> Tuple[pd.DataFrame, pd.Series]:
    """Load coefficient CSV(s) into a single DataFrame (rows: objectives, cols: samples).

    Expects each path to be a CSV where index or first column are objective IDs and
    columns are samples/conditions. Concatenates along columns and returns a label Series
    mapping each column to its group key (dict key).
    """
    dfs = []
    labels = []
    for key, path in coef_paths.items():
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Coefficient file not found: {p}")
        df = pd.read_csv(p, index_col=0)
        # normalize: columns are samples; index are objectives
        # Some files might have extra unnamed columns; drop empties
        df = df.loc[:, ~df.columns.str.fullmatch("")]
        dfs.append(df)
        labels.extend([(c, key) for c in df.columns])
    if not dfs:
        raise ValueError("No coefficient tables loaded")
    coef_df = pd.concat(dfs, axis=1)
    label_map = pd.Series({col: grp for col, grp in labels})
    # ensure alignment with DataFrame
    label_series = label_map.reindex(coef_df.columns).fillna("")
    return coef_df, label_series


def save_heatmap(coef_df: pd.DataFrame, save_path: Path, prefix: str, allocation_norm: bool = False) -> None:
    try:
        import seaborn as sns
        # Dynamic sizing to reduce label overlap
        X = coef_df.fillna(0)
        if allocation_norm:
            X = X.div(X.sum(axis=0), axis=1).fillna(0)
        nrows, ncols = X.shape
        width = min(max(10, ncols * 0.25), 40)
        height = min(max(8, nrows * 0.25), 40)
        g = sns.clustermap(
            X, cmap="viridis", figsize=(width, height),
            row_cluster=False, col_cluster=False, xticklabels=1, yticklabels=1
        )
        # Improve tick readability: show ALL labels with small font, rotate x
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, ha='center', fontsize=6)
        plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=6)
        g.figure.subplots_adjust(bottom=0.2, left=0.2)
        out = save_path / f"{prefix}_coef_clustermap.png"
        g.savefig(out, dpi=200)
        plt.close(g.figure)
    except Exception as e:
        # fallback: simple heatmap
        X = coef_df.fillna(0)
        if allocation_norm:
            X = X.div(X.sum(axis=0), axis=1).fillna(0)
        nrows, ncols = X.shape
        width = min(max(10, ncols * 0.25), 40)
        height = min(max(8, nrows * 0.25), 40)
        plt.figure(figsize=(width, height))
        plt.imshow(X, aspect='auto', cmap='viridis')
        plt.colorbar(label='Coefficient')
        plt.yticks(range(len(X.index)), X.index, fontsize=6)
        plt.xticks(range(len(X.columns)), X.columns, fontsize=6, rotation=90)
        out = save_path / f"{prefix}_coef_heatmap.png"
        plt.tight_layout()
        plt.savefig(out, dpi=200)
        plt.close()


def compute_entropy(coef_df: pd.DataFrame) -> pd.Series:
    X = coef_df.fillna(0).clip(lower=0)
    s = X.sum(axis=0)
    s[s == 0] = 1.0
    P = X / s
    with np.errstate(divide='ignore', invalid='ignore'):
        ent = -(P * np.log(P.where(P > 0, 1))).sum(axis=0)
    return ent


def _biomass_coefficients_from_gem(gem_path: str) -> Optional[pd.Series]:
    try:
        import cobra
        model = cobra.io.load_matlab_model(gem_path)
        r = model.reactions.get_by_id('biomass_objective')
        # Aggregate coefficients across compartments; strip bracket or underscore suffixes.
        agg: Dict[str, float] = {}
        for met, val in r.metabolites.items():
            raw = met.id
            # strip bracket-style compartment, e.g., atp[c]
            base = raw.split('[')[0]
            # also strip common underscore-style suffixes if present (e.g., atp_c)
            for suf in ('_c', '_m', '_n', '_r', '_g', '_l', '_x', '_p', '_e'):
                if base.endswith(suf):
                    base = base[: -len(suf)]
                    break
            if base == 'h2o':
                continue
            w = max(0.0, -float(val))
            if w <= 0:
                continue
            agg[base] = agg.get(base, 0.0) + w
        if not agg:
            return None
        return pd.Series(agg, name='coefficient')
    except Exception:
        return None


def compute_distance_to_biomass(coef_df: pd.DataFrame, gem_path: Optional[str] = None) -> Optional[pd.Series]:
    # Implement legacy-style distance using biomass composition and gh breakdown
    # 1) Build biomass composition vector from GEM if available
    biomass = None
    if gem_path:
        biomass = _biomass_coefficients_from_gem(gem_path)
        # biomass already excludes 'h2o' and zeros
    # 2) Prepare working coefficient table and drop 'h2o' only for distance
    work = coef_df.copy()
    if 'h2o' in work.index:
        work = work.drop(index='h2o')
    # 3) If we have a biomass composition, redistribute gh into biomass components
    if isinstance(biomass, pd.Series) and not biomass.empty:
        # Keep only biomass components for distance
        idx = [m for m in biomass.index if m in work.index]
        if len(idx) == 0:
            return None
        # Add gh contribution to each biomass component per sample, then drop gh
        if 'gh' in work.index:
            gh_row = work.loc['gh']
            # Broadcast: component_vector (len idx) outer-product gh_row (1 x n_samples)
            add_mat = pd.DataFrame(
                data=np.outer(biomass.loc[idx].to_numpy(), gh_row.to_numpy()),
                index=idx,
                columns=work.columns,
            )
            work.loc[idx, :] = work.loc[idx, :].add(add_mat, fill_value=0.0)
            work = work.drop(index='gh')
        # Reference is biomass coefficients over the same idx
        ref = biomass.loc[idx]
        # Compute distance across biomass components only, with normalization
        sub = work.loc[idx, :].fillna(0.0)
        # Column-normalize inferred coefficients (like legacy norm=True)
        colsum = sub.sum(axis=0).replace(0, 1.0)
        sub = sub.div(colsum, axis=1)
        # Normalize biomass reference to sum to 1 to match scale
        rs = float(ref.sum())
        if rs > 0:
            ref = ref / rs
        D = ((sub.subtract(ref, axis=0)) ** 2).sum(axis=0) ** 0.5
        return D
    # 4) Fallback: use a biomass-like row if present (gh/biomass), no breakdown
    names = work.index.astype(str)
    lower = names.str.lower()
    for s in ['biomass', 'biomass_objective', 'gh', 'growth']:
        if s in names.values:
            ref_row = s
            break
        elif s in lower.values:
            ref_row = names[lower == s].iloc[0]
            break
    else:
        ref_row = None
    if ref_row is None:
        return None
    # Fallback: normalize columns before distance to reduce scale effects
    colsum = work.sum(axis=0).replace(0, 1.0)
    work = work.div(colsum, axis=1)
    ref_vec = work.loc[ref_row]
    D = ((work.subtract(ref_vec, axis=1)) ** 2).sum(axis=0) ** 0.5
    return D


def _umap_safe_embed(X: np.ndarray, n_neighbors: int = 15, n_components: int = 2, random_state: int = 0, min_dist: float = 0.1) -> np.ndarray:
    """UMAP embedding with safety guard: cap components to N-1 to satisfy eigsh k < N.

    Parameters
    - X: samples x features
    - n_neighbors: will be capped to min(n_neighbors, N-1) where N = X.shape[0]
    - n_components: will be capped to min(max(2, n_components), max(2, N-1))
    """
    import umap
    N = int(X.shape[0])
    nn_req = int(n_neighbors)
    nc_req = int(n_components)
    nN = max(2, min(nn_req, max(2, N - 1)))
    nC = max(2, min(nc_req, max(2, N - 1)))
    init = 'random' if nC >= (N - 1) else 'spectral'
    if nC != nc_req or init == 'random':
        print(f"[UMAP] Adjusted n_components from {nc_req} to {nC} (N={N}); init='{init}'.", flush=True)
    if nN != nn_req:
        print(f"[UMAP] Adjusted n_neighbors from {nn_req} to {nN} (N={N}).", flush=True)
    reducer = umap.UMAP(n_components=nC, n_neighbors=nN, init=init, random_state=random_state, min_dist=min_dist)
    return reducer.fit_transform(X)


def reduce_scatter(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    save_path: Path,
    prefix: str,
    reduction: str = "auto",
    umap_params: Optional[Dict] = None,
    cluster: bool = False,
    return_emb: bool = False,
) -> Optional[pd.DataFrame]:
    X = coef_df.fillna(0).T
    emb = None
    title = ""
    red = (reduction or "auto").lower()
    if red == "umap" or red == "auto":
        try:
            params = umap_params or {}
            n_neighbors = int(params.get('n_neighbors', 15))
            n_components = int(params.get('n_components', 2))
            min_dist = float(params.get('min_dist', 0.1))
            random_state = int(params.get('random_state', 0))
            emb = _umap_safe_embed(X, n_neighbors=n_neighbors, n_components=n_components, random_state=random_state, min_dist=min_dist)
            title = "UMAP"
        except Exception:
            if red == "umap":
                # If forced umap but not available, fall back to PCA explicitly
                pass
    if emb is None:
        try:
            from sklearn.decomposition import PCA
            pca = PCA(n_components=2, random_state=0)
            emb = pca.fit_transform(X)
            title = "PCA"
        except Exception:
            # nothing to do
            return
    clusters = None
    if cluster:
        try:
            import hdbscan
            clusterer = hdbscan.HDBSCAN(min_cluster_size=5)
            clusters = clusterer.fit_predict(emb)
        except Exception:
            clusters = None

    # If embedding has more than 2 dimensions (e.g., n_components>2), use first two for plotting
    if emb.ndim == 2 and emb.shape[1] > 2:
        emb_plot = emb[:, :2]
    else:
        emb_plot = emb

    # Style to match regressionAnalyzer: larger markers, black edges, clean spines, no ticks
    sns.set_context("notebook", font_scale=2.0)
    sns.set_style("white")
    df = pd.DataFrame(emb_plot, columns=[f"{title}1", f"{title}2"])  # for consistent API
    df["label"] = labels.values

    # Choose palette similar to regressionAnalyzer
    unique_labels = np.unique(labels)
    palette = 'gist_rainbow' if len(unique_labels) > 5 else 'Set2'

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax = sns.scatterplot(
        data=df,
        x=f"{title}1",
        y=f"{title}2",
        hue="label",
        alpha=0.5,
        palette=palette,
        s=100,
        legend=True,
        linewidth=1,
        edgecolor='black',
        ax=ax,
    )
    for pos in ['right', 'top', 'bottom', 'left']:
        plt.gca().spines[pos].set_visible(False)
    plt.xticks([]); plt.yticks([])
    handles, lab = ax.get_legend_handles_labels()
    # Enlarge legend marker sizes
    for h in handles:
        try:
            h.set_sizes(h.get_sizes() * 10)
        except Exception:
            pass
    if len(unique_labels) > 5:
        plt.legend(frameon=False, bbox_to_anchor=(1, 1.1), loc='best', fontsize=20, handles=handles)
    else:
        plt.legend(frameon=False, fontsize=20, loc='best', handles=handles)
    plt.tight_layout()
    out = save_path / f"{prefix}_{title.lower()}_scatter.png"
    plt.savefig(out, dpi=200)
    plt.close()

    # Save embedding as CSV
    df = pd.DataFrame({
        f"{title}1": emb_plot[:, 0],
        f"{title}2": emb_plot[:, 1],
        "label": labels.values,
    })
    if clusters is not None:
        df["cluster"] = clusters
    df.to_csv(save_path / f"{prefix}_{title.lower()}_embedding.csv", index=False)

    if return_emb:
        return df
    return None


def compare_groups(coef_df: pd.DataFrame, labels: pd.Series, a: str, b: str) -> pd.DataFrame:
    from scipy.stats import mannwhitneyu
    A = coef_df.loc[:, labels == a]
    B = coef_df.loc[:, labels == b]
    rows = []
    for obj in coef_df.index:
        x = A.loc[obj].dropna().to_numpy()
        y = B.loc[obj].dropna().to_numpy()
        if len(x) >= 3 and len(y) >= 3:
            stat, p = mannwhitneyu(x, y, alternative='two-sided')
            rows.append((obj, np.mean(x), np.mean(y), p))
    return pd.DataFrame(rows, columns=['objective', f'mean_{a}', f'mean_{b}', 'p_value'])


def group_by_metabolite_type(coef_df: pd.DataFrame, mapping: Union[str, Path, Dict[str, str]]) -> pd.DataFrame:
    """Aggregate objective coefficients by metabolite type categories.

    mapping can be:
      - path to CSV with columns: objective, category (or first col index as objective, and a 'category' column)
      - dict: {objective_id: category}
    Returns a DataFrame indexed by category with summed coefficients per sample.
    """
    if isinstance(mapping, (str, Path)):
        mdf = pd.read_csv(mapping)
        if 'objective' not in mdf.columns:
            mdf.rename(columns={mdf.columns[0]: 'objective'}, inplace=True)
        if 'category' not in mdf.columns:
            raise ValueError("Mapping CSV must include a 'category' column")
        mapper = dict(zip(mdf['objective'].astype(str), mdf['category'].astype(str)))
    else:
        mapper = {str(k): str(v) for k, v in mapping.items()}

    cats = pd.Series(coef_df.index.astype(str)).map(mapper)
    # Drop objectives without mapping
    df = coef_df.copy()
    df['__category__'] = cats.values
    df = df.dropna(subset=['__category__'])
    grouped = df.groupby('__category__').sum(numeric_only=True)
    grouped.index.name = 'category'
    return grouped


def find_tradeoffs(coef_df: pd.DataFrame, method: str = 'spearman', top_k: int = 100, corr_threshold: float = -0.5) -> pd.DataFrame:
    """Identify pairs of objectives with strong negative correlation across samples.

    Returns a DataFrame with columns: obj_i, obj_j, correlation; limited to top_k most negative or below threshold.
    """
    corr = coef_df.T.corr(method=method)
    # Extract lower triangle without diagonal
    pairs = []
    idx = corr.index.tolist()
    for i in range(len(idx)):
        for j in range(i + 1, len(idx)):
            c = corr.iat[i, j]
            if c <= corr_threshold:
                pairs.append((idx[i], idx[j], c))
    # If none meet threshold, take top_k most negative
    if not pairs:
        for i in range(len(idx)):
            for j in range(i + 1, len(idx)):
                pairs.append((idx[i], idx[j], corr.iat[i, j]))
        pairs.sort(key=lambda x: x[2])
        pairs = pairs[:top_k]
    else:
        pairs.sort(key=lambda x: x[2])
    return pd.DataFrame(pairs, columns=['objective_i', 'objective_j', 'correlation'])


def plot_corr_heatmap(coef_df: pd.DataFrame, save_path: Path, prefix: str, title: str = 'Objective Correlation Heatmap', cmap: str = 'RdBu_r') -> None:
    corr = coef_df.T.corr(method='pearson')
    n = corr.shape[0]
    width = min(max(8, n * 0.3), 40)
    height = width
    plt.figure(figsize=(width, height))
    ax = sns.heatmap(corr, cmap=cmap, vmin=-1, vmax=1, square=True, cbar_kws={'label': 'Pearson r'})
    ax.set_title(title)
    plt.setp(ax.get_xticklabels(), rotation=90, ha='center', fontsize=6)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=6)
    plt.tight_layout()
    plt.savefig(save_path / f"{prefix}_corr_heatmap.png", dpi=200)
    plt.close()


def _read_sampling_objFlux_minimal(path: str, medium: str = '', allocation_norm: bool = True) -> Optional[pd.DataFrame]:
    if _HAS_LEGACY_PARETO and path:
        try:
            res = _legacy_read_sampling_objFlux(path=path, medium=medium or 'DMEMF12', allocation_norm=allocation_norm)
            return res
        except Exception:
            pass
    # Fallback: return None (no Pareto overlay)
    return None


def _pareto_frontier_2d(df: pd.DataFrame, x: str, y: str) -> pd.DataFrame:
    # Compute Pareto frontier for maximizing both x and y
    pts = df[[x, y]].to_numpy()
    keep = []
    for i in range(len(pts)):
        xi, yi = pts[i]
        dominated = False
        for j in range(len(pts)):
            if j == i:
                continue
            xj, yj = pts[j]
            if (xj >= xi and yj >= yi) and (xj > xi or yj > yi):
                dominated = True
                break
        if not dominated:
            keep.append(i)
    return df.iloc[keep]


def pareto_scatter(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    save_path: Path,
    prefix: str,
    pareto_mets: List[str],
    samplingFlux_path: Optional[str] = None,
    medium: str = '',
    allocation_norm: bool = True,
    pareto_line: str = 'connected',
) -> None:
    X = coef_df.copy().T  # samples x objectives
    if allocation_norm:
        s = X.sum(axis=1).replace({0: 1.0})
        X = X.div(s, axis=0)
    # Load reference feasible/pareto points if available
    res = _read_sampling_objFlux_minimal(samplingFlux_path or '', medium=medium, allocation_norm=allocation_norm)

    sns.set_context("notebook", font_scale=2.0)
    sns.set_style("white")
    unique_labels = np.unique(labels)
    palette = 'gist_rainbow' if len(unique_labels) > 5 else 'Set2'

    # Resolve human-friendly names to objective IDs present in X
    def _resolve_name(name: str) -> Optional[str]:
        n = str(name).strip().lower()
        # direct match
        for c in X.columns:
            if str(c).lower() == n:
                return c
        # synonyms
        syn = {
            'biomass': ['gh', 'biomass', 'biomass_objective', 'growth'],
            'cholesterol': ['chsterol', 'cholesterol'],
            'glutathione': ['gthrd', 'gthox', 'glutathione'],
        }
        if n in syn:
            for cand in syn[n]:
                for c in X.columns:
                    if str(c).lower() == cand:
                        return c
        return None

    resolved_targets: List[str] = []
    for y in pareto_mets:
        yy = _resolve_name(y)
        if yy is None:
            print(f"[tradeoff] Skipping target '{y}': not found in coefficients.")
            continue
        resolved_targets.append(yy)

    for y in resolved_targets:
        candidates = [m for m in X.columns if m != y]
        for met in candidates:
            if met not in X.columns or y not in X.columns:
                continue
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(1, 1, 1)
            ax = sns.scatterplot(
                x=X[met], y=X[y], hue=labels.values, palette=palette,
                s=100, alpha=0.5, linewidth=1, edgecolor='black', ax=ax
            )

            # Overlay Pareto frontier if reference available
            if res is not None and met in res.columns and y in res.columns:
                # Use legacy find_pareto when possible; fallback to simple frontier
                pareto_df = None
                if _HAS_LEGACY_PARETO:
                    try:
                        pareto_df, _ = _legacy_find_pareto(res, met, y)
                    except Exception:
                        pareto_df = None
                if pareto_df is None:
                    pareto_df = _pareto_frontier_2d(res[[met, y]].copy(), met, y)
                    pareto_df['cellType'] = 'Pareto'

                pareto_df = pareto_df.sort_values(by=[met])
                ax.plot(pareto_df[met], pareto_df[y], '--', color='k', label='Est. Pareto', linewidth=3.0)

                # Optional curve fit of Pareto
                if pareto_line == 'curvefit' and len(pareto_df) >= 2:
                    try:
                        def _quad(px, a, b, c):
                            return a * px + b * (px**2) + c
                        popt, _ = curve_fit(_quad, pareto_df[met], pareto_df[y], maxfev=5000)
                        xs = np.linspace(pareto_df[met].min(), pareto_df[met].max(), 100)
                        ys = _quad(xs, *popt)
                        ax.plot(xs, ys, '-', color='k', alpha=0.6, linewidth=2.0)
                    except Exception:
                        pass

            for pos in ['right', 'top', 'bottom', 'left']:
                plt.gca().spines[pos].set_visible(False)
            plt.xticks([]); plt.yticks([])
            handles, lab = ax.get_legend_handles_labels()
            for h in handles:
                try:
                    h.set_sizes(h.get_sizes() * 10)
                except Exception:
                    pass
            if len(unique_labels) > 5:
                plt.legend(frameon=False, bbox_to_anchor=(1, 1.1), loc='best', fontsize=20, handles=handles)
            else:
                plt.legend(frameon=False, fontsize=20, loc='best', handles=handles)
            plt.xlabel(met); plt.ylabel(y)
            plt.tight_layout()
            plt.savefig(save_path / f"{prefix}_pareto_{met}_vs_{y}.png", dpi=200)
            plt.close()


def trait_analysis_minimal(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    samplingFlux_path: Optional[str],
    save_path: Path,
    prefix: str,
    medium: str = '',
    allocation_norm: bool = True,
    n_pc: int = 2,
) -> None:
    # Load reference (controls)
    res = _read_sampling_objFlux_minimal(samplingFlux_path or '', medium=medium, allocation_norm=allocation_norm)
    if res is None:
        return
    # Prepare normalized matrices
    arch = res.copy().T  # objectives x samples
    arch = arch.iloc[:-1, :] if 'cellType' in arch.index else arch
    sel = coef_df.copy()
    if allocation_norm:
        sel = sel.div(sel.sum(axis=0), axis=1)
        arch = arch.div(arch.sum(axis=0), axis=1)
    # Align objectives
    common = sel.index.intersection(arch.index)
    arch = arch.loc[common]
    sel = sel.loc[common]
    # Build combined matrix: controls then each label mean
    ctrl = arch.copy()
    grp = sel.T.groupby(labels.values).mean().T
    merge = pd.concat([ctrl, grp], axis=1, join='outer').fillna(0)
    # Heatmap of archetype vs objectives
    plt.figure(figsize=(max(8, merge.shape[1]*0.6), 6))
    sns.heatmap(merge, cmap='Blues', cbar_kws={'shrink': 0.6})
    plt.title('Trait Analysis (Control vs Group Means)')
    plt.tight_layout()
    plt.savefig(save_path / f"{prefix}_trait_heatmap.png", dpi=200)
    plt.close()


# -----------------------------
# Additional polished plots
# -----------------------------

def _canvas_style(ax):
    for pos in ['right', 'top', 'bottom', 'left']:
        ax.spines[pos].set_visible(False)
    return ax


def plot_coeff_strip(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    save_path: Path,
    prefix: str,
    top_k: int = 12,
    normalized: bool = False,
    pv_thresh: float = 0.05,
    fc_thresh: float = 1.0,
) -> None:
    # Prepare matrix
    X = coef_df.copy()
    if normalized:
        # Legacy behavior: normalize each metabolite (row) by its max across samples
        X = X.div(X.max(axis=1), axis=0).fillna(0)
    # Select top_k objectives by overall mean
    means_all = X.mean(axis=1).sort_values(ascending=False)
    objs = list(means_all.index[:max(1, top_k)])
    # Compute per-objective p-values and fold-changes for two groups
    groups = list(pd.unique(labels))
    if len(groups) < 2:
        groups = groups + groups
    col1, col2 = groups[0], groups[1]
    def _stars(p):
        return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'n.s.'
    pvals, fcs, colors = [], [], []
    for obj in objs:
        x1 = X.loc[obj, labels == col1].astype(float)
        x2 = X.loc[obj, labels == col2].astype(float)
        try:
            # Nonparametric (legacy): Mann–Whitney U, two-sided
            _, p = mannwhitneyu(x2, x1, alternative='two-sided')
        except Exception:
            p = 1.0
        pvals.append(p)
        fc = (x2.mean() / max(1e-12, x1.mean())) if x1.mean() != 0 else np.inf
        fcs.append(fc)
        if p < pv_thresh and fc > fc_thresh:
            colors.append(sns.color_palette("Dark2")[1])
        elif p < pv_thresh and fc < 1.0 / max(1e-12, fc_thresh):
            colors.append(sns.color_palette("Dark2")[0])
        else:
            colors.append('k')
    # Filter to only significant metabolites (legacy behavior). If none, fall back to top_k selection.
    sig_mask = [(p < pv_thresh) and (fc > fc_thresh or fc < (1.0 / max(1e-12, fc_thresh))) for p, fc in zip(pvals, fcs)]
    if any(sig_mask):
        objs = [o for o, keep in zip(objs, sig_mask) if keep]
        pvals = [p for p, keep in zip(pvals, sig_mask) if keep]
        colors = [c for c, keep in zip(colors, sig_mask) if keep]

    # Build plotting dataframe and label objectives with stars
    obj_labels = [f"{o} \n {_stars(p)}" for o, p in zip(objs, pvals)]
    tidy = (
        X.loc[objs]
        .T.assign(label=labels.values)
        .melt(id_vars='label', var_name='objective', value_name='coefficient')
    )
    # Map objective names to labeled versions with stars
    label_map = dict(zip(objs, obj_labels))
    tidy['objective'] = tidy['objective'].map(label_map)
    # Plot: strip with Pastel2 and pointplot with Dark2 (means)
    sns.set_context("notebook", font_scale=2.0)
    sns.set_style("white")
    fig, ax = plt.subplots(1, 1, figsize=(max(10, len(objs)), 8))
    pal_strip = {k: c for k, c in zip(groups[:2], sns.color_palette("Pastel2")[:2])}
    pal_point = {k: c for k, c in zip(groups[:2], sns.color_palette("Dark2")[:2])}
    sns.stripplot(
        data=tidy, x='objective', y='coefficient', hue='label', dodge=True,
        jitter=0.25, s=6, palette=pal_strip, alpha=0.5, ax=ax, hue_order=groups[:2]
    )
    sns.pointplot(
        data=tidy, x='objective', y='coefficient', hue='label', dodge=0.35,
        join=False, markers='X', scale=1.0, ci=None, palette=pal_point,
        ax=ax, hue_order=groups[:2]
    )
    ax.set_xlabel('')
    ax.set_ylabel('Proportion' if normalized else 'Coefficient')
    # Color-code x tick labels per significance direction
    xt = ax.get_xticklabels()
    for i, tl in enumerate(xt):
        tl.set_color(colors[i] if i < len(colors) else 'k')
    plt.setp(ax.get_xticklabels(), rotation=90, ha='center')
    # Legend cleanup: keep one set
    handles, labs = ax.get_legend_handles_labels()
    bylab = dict(zip(labs, handles))
    ax.legend(bylab.values(), bylab.keys(), frameon=False, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10)
    plt.tight_layout()
    plt.savefig(save_path / f"{prefix}_coef_strip{'_norm' if normalized else ''}.png", dpi=200)
    plt.close(fig)


def plot_proportional_bars(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    save_path: Path,
    prefix: str,
    top_k: int = 12,
) -> None:
    # Normalize per sample, then average per group
    X = coef_df.div(coef_df.sum(axis=0), axis=1).fillna(0)
    groups = labels.unique()
    group_means = {g: X.loc[:, labels == g].mean(axis=1) for g in groups}
    mean_all = X.mean(axis=1).sort_values(ascending=False)
    objs = list(mean_all.index[:max(1, top_k)])
    # Add 'Other'
    palette = sns.color_palette('tab20', n_colors=len(objs) + 1)
    fig, ax = plt.subplots(1, 1, figsize=(max(8, len(groups) * 1.4), 6))
    bottoms = np.zeros(len(groups))
    for i, obj in enumerate(objs):
        vals = np.array([group_means[g].get(obj, 0.0) for g in groups])
        ax.bar(groups, vals, bottom=bottoms, color=palette[i], edgecolor='black', linewidth=0.5, label=obj)
        bottoms += vals
    # Other slice
    other = np.clip(1.0 - bottoms, 0, 1)
    ax.bar(groups, other, bottom=bottoms, color=palette[-1], edgecolor='black', linewidth=0.5, label='Other')
    ax.set_ylabel('Mean Proportion')
    _canvas_style(ax)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False, fontsize=8, ncol=1)
    plt.tight_layout()
    plt.savefig(save_path / f"{prefix}_coef_proportions.png", dpi=200)
    plt.close()


def plot_distance_to_biomass_fig(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    save_path: Path,
    prefix: str,
    gem_path: Optional[str] = None,
) -> None:
    # Reuse compute_distance_to_biomass to get distances
    dist = compute_distance_to_biomass(coef_df, gem_path=gem_path)
    if dist is None:
        return
    df = pd.DataFrame({'distance': dist.values, 'label': labels.values})
    sns.set_context("notebook", font_scale=1.8)
    sns.set_style("white")
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    sns.boxplot(data=df, x='label', y='distance', ax=ax, fliersize=0, color='white', linewidth=1, whis=1.5)
    sns.stripplot(data=df, x='label', y='distance', ax=ax, color='k', size=3, alpha=0.5, jitter=0.2)
    ax.set_xlabel('Group'); ax.set_ylabel('Distance to Biomass')
    ax.tick_params(axis='x', rotation=45, labelsize=8)
    _canvas_style(ax)
    plt.tight_layout()
    plt.savefig(save_path / f"{prefix}_distance_to_biomass.png", dpi=200)
    plt.close()

def plot_distance_with_hist_pvalue(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    save_path: Path,
    prefix: str,
    gem_path: Optional[str] = None,
    bins: int = 30,
) -> None:
    dist = compute_distance_to_biomass(coef_df, gem_path=gem_path)
    if dist is None:
        return
    df = pd.DataFrame({'distance': dist.values, 'label': labels.values})
    groups = list(pd.unique(df['label']))
    if len(groups) < 2:
        groups = groups + groups
    g1, g2 = groups[0], groups[1]
    x1 = df.loc[df['label'] == g1, 'distance']
    x2 = df.loc[df['label'] == g2, 'distance']
    try:
        _, pval = ttest_ind(x1, x2, equal_var=False, nan_policy='omit')
    except Exception:
        pval = np.nan

    sns.set_context("notebook", font_scale=1.8)
    sns.set_style("white")
    # Stacked panels: histogram (top, taller) + scatter (bottom, short & wide)
    fig = plt.figure(figsize=(10, 7))
    gs = fig.add_gridspec(2, 1, height_ratios=[2, 1], hspace=0.05)
    ax_hist = fig.add_subplot(gs[0])
    ax_scat = fig.add_subplot(gs[1], sharex=ax_hist)

    # Legacy-like: density hist with Pastel2 palette
    tmp_df = pd.DataFrame({'x': pd.concat([x1, x2], ignore_index=True),
                           'label': [str(g1)] * len(x1) + [str(g2)] * len(x2)})
    sns.histplot(data=tmp_df, x='x', hue='label', stat='density', kde=True,
                 palette='Pastel2', common_norm=False, bins=100, element='step', ax=ax_hist)
    ax_hist.set_ylabel('Density')
    ax_hist.grid(False)

    # Two narrow jitter bands for bottom scatter
    y1 = np.zeros_like(x1, dtype=float) + 0.1 + (np.random.rand(len(x1)) - 0.5) * 0.02
    y2 = np.zeros_like(x2, dtype=float) - 0.1 + (np.random.rand(len(x2)) - 0.5) * 0.02
    ax_scat.scatter(x1, y1, s=30, color='slateblue', edgecolor='black', alpha=0.8, label=str(g1))
    ax_scat.scatter(x2, y2, s=30, color='tomato', edgecolor='black', alpha=0.8, label=str(g2))
    ax_scat.set_yticks([])
    ax_scat.set_xlabel('Distance to Biomass Objective')
    ax_scat.grid(False)
    if not np.isnan(pval):
        ax_hist.set_title(f'P-value:{pval:.5f}', fontsize=18)
    # Keep frames for scatter; soften histogram spines
    for pos in ['right', 'top', 'bottom', 'left']:
        ax_hist.spines[pos].set_visible(False)
    if not np.isnan(pval):
        fig.suptitle(f'P-value: {pval:.5f}')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(save_path / f"{prefix}_distance_hist_scatter.png", dpi=200)
    plt.close(fig)


def plot_allocation_heatmap(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    save_path: Path,
    prefix: str,
    top_k: int = 25,
) -> None:
    # Normalize per sample and aggregate per group (mean)
    X = coef_df.div(coef_df.sum(axis=0), axis=1).fillna(0)
    grp_mean = X.T.groupby(labels.values).mean().T
    mean_all = grp_mean.mean(axis=1).sort_values(ascending=False)
    objs = list(mean_all.index[:max(1, top_k)])
    M = grp_mean.loc[objs]
    # Heatmap with readable ticks
    nrows, ncols = M.shape
    width = min(max(8, ncols * 0.6), 40)
    height = min(max(8, nrows * 0.35), 40)
    plt.figure(figsize=(width, height))
    ax = sns.heatmap(M, cmap='Blues', vmin=0, vmax=min(1.0, M.values.max() or 1.0), cbar_kws={'shrink': 0.6})
    plt.setp(ax.get_xticklabels(), rotation=90, ha='center', fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path / f"{prefix}_allocation_heatmap.png", dpi=200)
    plt.close()

def plot_lollipop_proportions(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    save_path: Path,
    prefix: str,
    threshold: float = 0.0,
    top_k: Optional[int] = None,
) -> None:
    # Legacy-like lollipop: vertical vlines and Dark2 palette
    X = coef_df.fillna(0)
    groups = list(pd.unique(labels))
    prop_by_group = {g: (X.loc[:, labels == g] > threshold).mean(axis=1) for g in groups}
    prop_df = pd.DataFrame(prop_by_group).fillna(0.0)
    prop_df['Objective'] = prop_df.index
    if groups:
        prop_df = prop_df.sort_values(by=[groups[0]])
    if top_k:
        prop_df = prop_df.tail(top_k)
    plot_df = prop_df.melt(id_vars=['Objective'], var_name='cellType', value_name='Proportion')
    mins = plot_df.groupby('Objective')['Proportion'].min()
    maxs = plot_df.groupby('Objective')['Proportion'].max()
    sns.set_context("notebook", font_scale=2.0)
    sns.set_style("white")
    fig, ax = plt.subplots(1, 1, figsize=(max(8, len(prop_df) / 2), 6))
    ax.vlines(x=prop_df['Objective'], ymin=mins[prop_df['Objective']].values, ymax=maxs[prop_df['Objective']].values,
              linewidth=3, color='grey')
    palette = sns.color_palette("Dark2", n_colors=len(groups) or 2)
    sns.scatterplot(data=plot_df, x='Objective', y='Proportion', hue='cellType', palette=palette,
                    s=200, zorder=7, ax=ax, hue_order=groups)
    plt.xticks(rotation=90)
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel('')
    for pos in ['right', 'top', 'bottom', 'left']:
        ax.spines[pos].set_visible(False)
    plt.tight_layout()
    plt.savefig(save_path / f"{prefix}_proportion_lollipop.png", dpi=200)
    plt.close(fig)

def plot_allocation_bar_distribution(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    save_path: Path,
    prefix: str,
    mets: Optional[List[str]] = None,
    top_k: int = 10,
) -> None:
    X = coef_df.div(coef_df.sum(axis=0), axis=1).fillna(0)
    if mets is None:
        mean_all = X.mean(axis=1).sort_values(ascending=False)
        mets = list(mean_all.index[:max(1, top_k)])
    bp_df = X.loc[mets].copy()
    bp_df.columns = labels.values
    bp_df['metabolites'] = bp_df.index
    tidy = bp_df.melt(id_vars=['metabolites'])
    tidy.columns = ['metabolites', 'cellType', 'allocation']
    sns.set_context("notebook", font_scale=2.0)
    fig, ax = plt.subplots(1, 1, figsize=(max(8, len(mets) / 2), 6))
    colors = ['teal', 'tomato', 'slateblue', 'lightgrey']
    cp = sns.color_palette(colors[:len(pd.unique(labels))])
    sns.barplot(data=tidy, x='metabolites', y='allocation', hue='cellType', palette=cp, alpha=0.5,
                ax=ax, hue_order=list(pd.unique(labels)))
    plt.xticks(rotation=90)
    for pos in ['right', 'top', 'bottom', 'left']:
        ax.spines[pos].set_visible(False)
    plt.tight_layout()
    plt.savefig(save_path / f"{prefix}_allocation_bar_distribution.png", dpi=200)
    plt.close(fig)


def plot_tradeoff_pairs(
    coef_df: pd.DataFrame,
    labels: pd.Series,
    mets: List[str],
    save_path: Path,
    prefix: str,
    normalize: bool = True,
) -> None:
    """Plot pairwise scatter for the given metabolite objectives.

    If normalize is True, scales each sample's objective vector to sum to 1 before plotting.
    """
    X = coef_df.copy().T  # samples x objectives
    if normalize:
        s = X.sum(axis=1).replace({0: 1.0})
        X = X.div(s, axis=0)
    # subset available metabolites
    mets = [m for m in mets if m in X.columns]
    if len(mets) < 2:
        return
    palette = 'gist_rainbow' if len(np.unique(labels)) > 5 else 'Set2'
    for i in range(len(mets)):
        for j in range(i + 1, len(mets)):
            xi, xj = mets[i], mets[j]
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(1, 1, 1)
            ax = sns.scatterplot(
                x=X[xi], y=X[xj], hue=labels.values, palette=palette,
                s=100, alpha=0.5, linewidth=1, edgecolor='black', ax=ax
            )
            for pos in ['right', 'top', 'bottom', 'left']:
                plt.gca().spines[pos].set_visible(False)
            plt.xticks([]); plt.yticks([])
            handles, lab = ax.get_legend_handles_labels()
            for h in handles:
                try:
                    h.set_sizes(h.get_sizes() * 10)
                except Exception:
                    pass
            if len(np.unique(labels)) > 5:
                plt.legend(frameon=False, bbox_to_anchor=(1, 1.1), loc='best', fontsize=20, handles=handles)
            else:
                plt.legend(frameon=False, fontsize=20, loc='best', handles=handles)
            plt.xlabel(xi); plt.ylabel(xj)
            plt.tight_layout()
            plt.savefig(save_path / f"{prefix}_tradeoff_{xi}_vs_{xj}.png", dpi=200)
            plt.close()


def run_minimal_analysis(cfg: dict) -> None:
    save_root = Path(cfg.get('save_root_path', './examples/analyze_demo/out/'))
    save_root.mkdir(parents=True, exist_ok=True)
    prefix = cfg.get('prefix', 'analyze_demo')
    reduction = (cfg.get('reduction') or 'auto')
    # Flags similar to coefAnalysis API
    do_clustering = bool(cfg.get('clustering', True))
    do_entropy = bool(cfg.get('entropy', True))
    do_distance = bool(cfg.get('distance', True))
    do_tradeoff = bool(cfg.get('tradeoff', False))
    do_metType_cluster = bool(cfg.get('metType_cluster', False))
    cluster_umap = bool(cfg.get('cluster_umap', False))

    coef_df, labels = load_coefficients(cfg.get('coef_paths', {}))
    # Drop columns with all-zero coefficients to avoid degenerate plots/metrics
    nz_mask = coef_df.fillna(0).any(axis=0)
    if (~nz_mask).any():
        coef_df = coef_df.loc[:, nz_mask]
        labels = labels.loc[nz_mask]
    # Optional label overrides via config
    labels_cfg = cfg.get('labels', {}) or {}
    mode = str(labels_cfg.get('mode', '')).lower()
    if mode == 'column':
        labels = pd.Series(coef_df.columns, index=coef_df.columns)
    elif mode == 'regex':
        import re
        pattern = labels_cfg.get('regex', '(.*)')
        group_idx = int(labels_cfg.get('regex_group', 1))
        rx = re.compile(pattern)
        def _rxlab(c):
            m = rx.search(str(c))
            if not m:
                return str(c)
            try:
                return m.group(group_idx)
            except Exception:
                return m.group(0)
        labels = pd.Series([_rxlab(c) for c in coef_df.columns], index=coef_df.columns)
    elif mode == 'split':
        delim = labels_cfg.get('split_delim', '_')
        idx = int(labels_cfg.get('split_index', 0))
        def _splitlab(c):
            parts = str(c).split(delim)
            return parts[idx] if 0 <= idx < len(parts) else str(c)
        labels = pd.Series([_splitlab(c) for c in coef_df.columns], index=coef_df.columns)
    elif mode == 'contains':
        pat = str(labels_cfg.get('contains_pattern', ''))
        true_label = str(labels_cfg.get('true_label', 'True'))
        false_label = str(labels_cfg.get('false_label', 'False'))
        def _contains(c):
            return true_label if pat and (pat in str(c)) else false_label
        labels = pd.Series([_contains(c) for c in coef_df.columns], index=coef_df.columns)
    # Optional group remap for group_key-derived labels
    group_map = labels_cfg.get('group_map') or {}
    if isinstance(group_map, dict) and len(group_map) > 0:
        labels = labels.map(lambda x: group_map.get(str(x), x))
    (save_root / f"{prefix}_coef_table.csv").write_text(coef_df.to_csv())

    # Heatmap
    if do_clustering:
        save_heatmap(
            coef_df,
            save_root,
            prefix,
            allocation_norm=bool(cfg.get('clustergram_allocation', False))
        )

    # Optional correlation heatmap (tradeoff view)
    if bool(cfg.get('corr_heatmap', False)) or bool(cfg.get('tradeoff', False)):
        plot_corr_heatmap(coef_df, save_root, prefix)

    # Entropy
    if do_entropy:
        ent = compute_entropy(coef_df)
        ent.to_csv(save_root / f"{prefix}_entropy.csv")

    # Distance
    if do_distance:
        dist = compute_distance_to_biomass(coef_df, gem_path=str(cfg.get('GEM_path') or ''))
        if dist is not None:
            dist.to_csv(save_root / f"{prefix}_distance_to_biomass.csv")

    # Reduction (auto tries umap; falls back to pca)
    umap_para = cfg.get('umap_para')
    umap_params = None
    if isinstance(umap_para, (list, tuple)) and len(umap_para) >= 2:
        umap_params = {'n_neighbors': umap_para[0], 'n_components': umap_para[1]}
    elif isinstance(cfg.get('umap_params'), dict):
        umap_params = cfg.get('umap_params')
    reduce_scatter(coef_df, labels, save_root, prefix, reduction=reduction, umap_params=umap_params, cluster=cluster_umap)

    # Optional comparison
    ca = cfg.get('coef_analysis', {})
    ref = ca.get('ref_col')
    tgt = ca.get('tgt_col')
    if ref and tgt and ref in labels.values and tgt in labels.values:
        cmp_df = compare_groups(coef_df, labels, ref, tgt)
        cmp_df.to_csv(save_root / f"{prefix}_compare_{ref}_vs_{tgt}.csv", index=False)

    # Built-in optional plots (polished minimal versions)
    plots_cfg = cfg.get('plots', {}) or {}
    if bool(plots_cfg.get('coef_strip', False)):
        plot_coeff_strip(coef_df, labels, save_root, prefix, top_k=int(plots_cfg.get('top_k', 12)), normalized=False)
    if bool(plots_cfg.get('coef_strip_norm', False)):
        plot_coeff_strip(coef_df, labels, save_root, prefix, top_k=int(plots_cfg.get('top_k', 12)), normalized=True)
    if bool(plots_cfg.get('proportional', False)):
        plot_proportional_bars(coef_df, labels, save_root, prefix, top_k=int(plots_cfg.get('top_k', 12)))
    if bool(plots_cfg.get('distance_plot', False)):
        plot_distance_to_biomass_fig(coef_df, labels, save_root, prefix, gem_path=str(cfg.get('GEM_path') or ''))
        # legacy-style overlay: histogram + scatter with p-value
        plot_distance_with_hist_pvalue(coef_df, labels, save_root, prefix, gem_path=str(cfg.get('GEM_path') or ''))
    if bool(plots_cfg.get('allocation_heatmap', False)):
        plot_allocation_heatmap(coef_df, labels, save_root, prefix, top_k=int(plots_cfg.get('top_k', 25)))
    if bool(plots_cfg.get('proportion_lollipop', False)):
        plot_lollipop_proportions(coef_df, labels, save_root, prefix, threshold=0.0, top_k=int(plots_cfg.get('top_k', 25)))
    if bool(plots_cfg.get('allocation_bar_distribution', False)):
        plot_allocation_bar_distribution(coef_df, labels, save_root, prefix, top_k=int(plots_cfg.get('top_k', 10)))

    # Optional: metabolite type aggregation
    if do_metType_cluster and cfg.get('metType_map'):
        grouped = group_by_metabolite_type(coef_df, cfg.get('metType_map'))
        grouped.to_csv(save_root / f"{prefix}_metType_grouped.csv")
        # heatmap for grouped
        try:
            import seaborn as sns
            g = sns.clustermap(grouped.fillna(0), cmap="viridis", figsize=(10, 8))
            g.savefig(save_root / f"{prefix}_metType_clustermap.png")
        except Exception:
            plt.figure(figsize=(10, 6))
            plt.imshow(grouped.fillna(0), aspect='auto', cmap='viridis')
            plt.colorbar(label='Coefficient')
            plt.tight_layout()
            plt.savefig(save_root / f"{prefix}_metType_heatmap.png", dpi=200)
            plt.close()

    # Optional: trade-off analysis and Pareto overlays
    if do_tradeoff:
        trade = find_tradeoffs(coef_df, method='spearman', top_k=int(cfg.get('tradeoff_top_k', 100)), corr_threshold=float(cfg.get('tradeoff_corr_threshold', -0.5)))
        trade.to_csv(save_root / f"{prefix}_tradeoffs.csv", index=False)
        pareto_mets = cfg.get('pareto_mets') or []
        if isinstance(pareto_mets, (list, tuple)) and len(pareto_mets) > 0:
            pareto_scatter(
                coef_df,
                labels,
                save_path=save_root,
                prefix=prefix,
                pareto_mets=list(pareto_mets),
                samplingFlux_path=cfg.get('samplingFlux_path'),
                medium=str(cfg.get('medium', '')),
                allocation_norm=True,
                pareto_line=str(cfg.get('pareto_line', 'connected')),
            )

    # Optional: trait analysis (control vs group means heatmap)
    if bool(cfg.get('traitAnalysis', False)):
        trait_analysis_minimal(
            coef_df,
            labels,
            samplingFlux_path=cfg.get('samplingFlux_path'),
            save_path=save_root,
            prefix=prefix,
            medium=str(cfg.get('medium', '')),
            allocation_norm=True,
            n_pc=int(cfg.get('trait_pc', 2)),
        )
