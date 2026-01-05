"""
analysis_runner.py
Configuration-driven wrapper around metObjAnalyzer to analyze metabolic objectives.
"""
import json
import argparse
from pathlib import Path
from typing import Dict, Any

import pandas as pd
from scooti.metObjAnalyzer import metObjAnalyzer
from scooti.analysis_minimal import run_minimal_analysis


def _load_json(p: str) -> Dict[str, Any]:
    with open(p, "r") as f:
        return json.load(f)


def _make_label_func(labels_cfg: Dict[str, Any]):
    mode = str(labels_cfg.get('mode', '')).lower()
    group_map = labels_cfg.get('group_map') or {}
    def _apply_group_map(series: pd.Series) -> pd.Series:
        if isinstance(group_map, dict) and len(group_map) > 0:
            return series.map(lambda x: group_map.get(str(x), x))
        return series
    if mode == 'column':
        def _lf(df: pd.DataFrame) -> pd.Series:
            s = pd.Series(df.columns, index=df.columns)
            return _apply_group_map(s)
        return _lf
    if mode == 'regex':
        import re
        pattern = labels_cfg.get('regex', '(.*)')
        group_idx = int(labels_cfg.get('regex_group', 1))
        rx = re.compile(pattern)
        def _lf(df: pd.DataFrame) -> pd.Series:
            cols = list(df.columns)
            out = []
            for c in cols:
                m = rx.search(str(c))
                if not m:
                    out.append(str(c))
                else:
                    try:
                        out.append(m.group(group_idx))
                    except Exception:
                        out.append(m.group(0))
            s = pd.Series(out, index=cols)
            return _apply_group_map(s)
        return _lf
    if mode == 'split':
        delim = labels_cfg.get('split_delim', '_')
        idx = int(labels_cfg.get('split_index', 0))
        def _lf(df: pd.DataFrame) -> pd.Series:
            cols = list(df.columns)
            out = []
            for c in cols:
                parts = str(c).split(delim)
                out.append(parts[idx] if 0 <= idx < len(parts) else str(c))
            s = pd.Series(out, index=cols)
            return _apply_group_map(s)
        return _lf
    if mode == 'contains':
        pat = str(labels_cfg.get('contains_pattern', ''))
        true_label = str(labels_cfg.get('true_label', 'True'))
        false_label = str(labels_cfg.get('false_label', 'False'))
        def _lf(df: pd.DataFrame) -> pd.Series:
            cols = list(df.columns)
            out = [true_label if pat and (pat in str(c)) else false_label for c in cols]
            s = pd.Series(out, index=cols)
            return _apply_group_map(s)
        return _lf
    return None


def main():
    ap = argparse.ArgumentParser(description="Run metabolic objective analyses from a JSON config")
    ap.add_argument("--config", required=True, help="Path to analyze config JSON")
    args = ap.parse_args()

    cfg = _load_json(args.config)

    flux_paths = cfg.get("flux_paths", {})
    coef_paths = cfg.get("coef_paths", {})
    save_root_path = cfg.get("save_root_path", "./examples/analyze_demo/out/")
    GEM_path = cfg.get("GEM_path", "./scooti/metabolicModel/GEMs/Shen2019.mat")
    uncon_model_path = cfg.get("uncon_model_path", "")
    col_map = cfg.get("col_map", {})
    # Default: let analyzers derive labels from group keys (coef_paths) when label_func is None
    label_func = None
    # Optional: build a label_func from config for legacy engine
    labels_cfg = cfg.get('labels', {}) or {}
    if isinstance(labels_cfg, dict) and labels_cfg.get('mode'):
        lf = _make_label_func(labels_cfg)
        if lf is not None:
            label_func = lf
    samplingFlux_path = cfg.get("samplingFlux_path", "")
    sel_para = cfg.get("sel_para", "")
    prefix = cfg.get("prefix", "analyze_demo")
    medium = cfg.get("medium", "DMEMF12")

    # Resolve coef_paths entries: allow directory input by auto-selecting a CSV
    resolved_coef_paths = {}
    for name, path in coef_paths.items():
        p = Path(path)
        if p.is_dir():
            # Prefer flux_sl_*.csv then any *.csv
            cand = sorted(p.glob("flux_sl_*.csv")) or sorted(p.glob("*.csv"))
            if cand:
                resolved_coef_paths[name] = str(cand[0])
                print(f"[analysis] coef_paths[{name}] -> {resolved_coef_paths[name]}")
            else:
                resolved_coef_paths[name] = str(p)
                print(f"[analysis][WARN] No CSV found under {p}; keeping directory path (may fail): {p}")
        else:
            resolved_coef_paths[name] = str(p)
            print(f"[analysis] coef_paths[{name}] = {resolved_coef_paths[name]}")

    # Ensure save path exists
    Path(save_root_path).mkdir(parents=True, exist_ok=True)
    print(f"[analysis] init: save_root_path={save_root_path} medium={medium} prefix={prefix}")

    # Persist resolved coef paths back into config so minimal engine can use them
    if resolved_coef_paths:
        cfg["coef_paths"] = resolved_coef_paths

    # Engine selection (minimal avoids optional plotting deps)
    engine = cfg.get("engine", "minimal").lower()
    if engine == "minimal":
        run_minimal_analysis(cfg)
        return
    moa = metObjAnalyzer(
        flux_paths=flux_paths,
        coef_paths=resolved_coef_paths,
        save_root_path=save_root_path,
        GEM_path=GEM_path,
        uncon_model_path=uncon_model_path,
        col_map=col_map,
        label_func=label_func,
        samplingFlux_path=samplingFlux_path,
        sel_para=sel_para,
        prefix=prefix,
        medium=medium,
    )

    # Optional steps
    if "get_flux" in cfg:
        gf = cfg["get_flux"]
        print("[analysis] get_flux")
        moa.get_flux(
            kappa=gf.get("kappa", [10, 1, 0.1]),
            rho=gf.get("rho", [10, 1, 0.1]),
            rank=bool(gf.get("rank", False)),
        )

    if "flux_analysis" in cfg:
        fa = cfg["flux_analysis"]
        print("[analysis] fluxAnalysis")
        moa.fluxAnalysis(
            kappa_arr=fa.get("kappaArr", [10, 1, 0.1]),
            rho_arr=fa.get("rhoArr", [10, 1, 0.1]),
        )

    if "get_coef" in cfg:
        gc = cfg["get_coef"]
        print("[analysis] get_coef")
        moa.get_coef(metType_cluster=bool(gc.get("metType_cluster", False)))
        # If no label_func was built, allow post-hoc override for backward compatibility
        if label_func is None and isinstance(labels_cfg, dict) and labels_cfg.get('mode'):
            mode = str(labels_cfg.get('mode', '')).lower()
            if mode in ('column', 'regex', 'split', 'contains'):
                import re
                cols = list(moa.coef_df.columns)
                if mode == 'column':
                    moa.labels = pd.Series(cols, index=cols)
                elif mode == 'regex':
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
                    moa.labels = pd.Series([_rxlab(c) for c in cols], index=cols)
                elif mode == 'split':
                    delim = labels_cfg.get('split_delim', '_')
                    idx = int(labels_cfg.get('split_index', 0))
                    def _splitlab(c):
                        parts = str(c).split(delim)
                        return parts[idx] if 0 <= idx < len(parts) else str(c)
                    moa.labels = pd.Series([_splitlab(c) for c in cols], index=cols)
                elif mode == 'contains':
                    pat = str(labels_cfg.get('contains_pattern', ''))
                    true_label = str(labels_cfg.get('true_label', 'True'))
                    false_label = str(labels_cfg.get('false_label', 'False'))
                    def _contains(c):
                        return true_label if pat and (pat in str(c)) else false_label
                    moa.labels = pd.Series([_contains(c) for c in cols], index=cols)
            group_map = labels_cfg.get('group_map') or {}
            if isinstance(group_map, dict) and len(group_map) > 0:
                moa.labels = moa.labels.map(lambda x: group_map.get(str(x), x))
        # Log labels summary for clarity
        try:
            uniq = pd.unique(moa.labels)
            print(f"[analysis] labels (unique): {list(map(str, uniq))}")
        except Exception:
            pass

    if "coef_analysis" in cfg:
        ca = dict(cfg["coef_analysis"])  # shallow copy
        print("[analysis] coefAnalysis")
        # Adjust UMAP parameters to avoid k >= N errors for small sample counts
        try:
            n_samples = int(moa.coef_df.shape[1]) if hasattr(moa, 'coef_df') and isinstance(moa.coef_df, pd.DataFrame) else 0
        except Exception:
            n_samples = 0
        up = ca.get("umap_para", [5, 50])
        if not isinstance(up, (list, tuple)) or len(up) != 2:
            up = [5, 50]
        # neighbors and components should be < n_samples, and at least 2
        if n_samples >= 3:
            max_k = max(2, n_samples - 1)
            neighbors = int(min(max_k, max(2, int(up[0]))))
            components = int(min(max_k, max(2, int(up[1]))))
            if components != up[1] or neighbors != up[0]:
                print(f"[analysis] Adjusting umap_para from {up} to [{neighbors}, {components}] for n={n_samples}")
            up = [neighbors, components]
        else:
            # too few samples; disable clustering to avoid errors
            ca["clustering"] = False
        moa.coefAnalysis(
            unknown_clustering=bool(ca.get("unknown_clustering", False)),
            clustering=bool(ca.get("clustering", False)),
            entropy=bool(ca.get("entropy", False)),
            distance=bool(ca.get("distance", True)),
            compare=bool(ca.get("compare", True)),
            umap_para=up,
            method=ca.get("method", "average"),
            ref_col=ca.get("ref_col", None),
        )

    # Optional: tradeoff analysis (legacy engine)
    if engine != "minimal" and "tradeoff_analysis" in cfg:
        ta = cfg["tradeoff_analysis"]
        # Skip sampling-dependent analyses if samplingFlux_path is empty or invalid
        try:
            sp = cfg.get("samplingFlux_path", "") or ''
        except Exception:
            sp = ''
        try:
            from pathlib import Path as _P
            if not sp or not _P(sp).exists():
                if 'pareto2DAnalysis' in ta and ta.get('pareto2DAnalysis'):
                    print('[analysis] samplingFlux_path missing; disabling pareto2DAnalysis')
                if 'traitAnalysis' in ta and ta.get('traitAnalysis'):
                    print('[analysis] samplingFlux_path missing; disabling traitAnalysis')
                ta['pareto2DAnalysis'] = False
                ta['traitAnalysis'] = False
        except Exception:
            pass
        print("[analysis] tradeoff_analysis")
        # Ensure coefficients are loaded for correlation and Pareto analyses
        try:
            has_coef = hasattr(moa, 'coef_df') and isinstance(moa.coef_df, pd.DataFrame) and moa.coef_df.shape[1] > 0
        except Exception:
            has_coef = False
        if not has_coef:
            try:
                moa.get_coef(metType_cluster=False)
            except Exception:
                pass
        moa.tradeoff_analysis(
            input_type=ta.get("table", "coef"),
            corr_th=ta.get("corr_th", 0.0),
            allocation_norm=bool(ta.get("normalize", True)),
            pareto_mets=ta.get("pareto_mets", ["gh"]),
            pareto3D=bool(ta.get("pareto3D", False)),
            pareto2DAnalysis=bool(ta.get("pareto2DAnalysis", False)),
            # Trait-level summaries are enabled via archetype analysis flags in this API
            archetypeAnalysis_w_control=bool(ta.get("traitAnalysis", False)),
        )
        # Optional correlation heatmap of objectives used in tradeoff analysis
        try:
            if bool(ta.get("corr_heatmap", True)):
                import numpy as np
                import seaborn as sns
                import matplotlib.pyplot as plt
                if ta.get("table", "coef") == "coef" and hasattr(moa, 'coef_df'):
                    df = moa.coef_df.copy()
                elif ta.get("table", "coef") == "fluxRecon" and hasattr(moa, 'fluxRecon_df'):
                    df = moa.fluxRecon_df.copy()
                else:
                    df = None
                if df is not None and isinstance(df, pd.DataFrame) and df.shape[0] > 1:
                    corr = df.T.corr(method='pearson')
                    n = corr.shape[0]
                    width = min(max(8, n * 0.3), 40)
                    height = width
                    plt.figure(figsize=(width, height))
                    ax = sns.heatmap(corr, cmap='RdBu_r', vmin=-1, vmax=1, square=True,
                                     cbar_kws={'label': 'Pearson r'})
                    ax.set_title('Objective Correlation Heatmap')
                    plt.setp(ax.get_xticklabels(), rotation=90, ha='center', fontsize=6)
                    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=6)
                    plt.tight_layout()
                    Path(save_root_path).mkdir(parents=True, exist_ok=True)
                    plt.savefig(Path(save_root_path) / f"{prefix}_tradeoff_corr_heatmap.png", dpi=200)
                    plt.close()
        except Exception:
            pass


if __name__ == "__main__":
    main()
