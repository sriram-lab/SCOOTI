import json, os, re, sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from statsmodels.stats.multitest import fdrcorrection
from scooti.utils.findSigGenes import findSigGenes
from pathlib import Path

# Utilities

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def ttest_and_fdr(exp_vals: np.ndarray, ctrl_vals: np.ndarray):
    import scipy.stats as ss
    pvals = []
    fcs = []
    for i in range(exp_vals.shape[0]):
        _, p = ss.ttest_ind(exp_vals[i], ctrl_vals[i])
        if np.isnan(p):
            p = 1.0
        pvals.append(p)
        ctrl_mean = np.mean(ctrl_vals[i])
        exp_mean = np.mean(exp_vals[i])
        fc = exp_mean / ctrl_mean if ctrl_mean != 0 else np.inf
        fcs.append(fc)
    reject, pvals_fdr = fdrcorrection(pvals)
    # Fix inf/nan fold-changes
    fcs = np.array(fcs, dtype=float)
    if np.any(np.isinf(fcs)):
        finite = fcs[np.isfinite(fcs)]
        if finite.size:
            fcs[np.isinf(fcs)] = np.max(finite)
        else:
            fcs[np.isinf(fcs)] = 1.0
    fcs[np.isnan(fcs)] = 1.0
    return reject, pvals_fdr, fcs


def save_up_dw(index, up_mask, dw_mask, out_dir, prefix):
    ensure_dir(out_dir)
    up_genes = (index[up_mask] if hasattr(index, '__getitem__') else np.array(index)[up_mask])
    dw_genes = (index[dw_mask] if hasattr(index, '__getitem__') else np.array(index)[dw_mask])
    pd.DataFrame({"upgenes": up_genes}).to_csv(os.path.join(out_dir, f"{prefix}_upgenes.csv"), index=False)
    pd.DataFrame({"dwgenes": dw_genes}).to_csv(os.path.join(out_dir, f"{prefix}_dwgenes.csv"), index=False)
    # Reverse (P)
    pd.DataFrame({"dwgenes": up_genes}).to_csv(os.path.join(out_dir, f"{prefix}_P_dwgenes.csv"), index=False)
    pd.DataFrame({"upgenes": dw_genes}).to_csv(os.path.join(out_dir, f"{prefix}_P_upgenes.csv"), index=False)


def run_sharma(cfg):
    path = cfg["sharma_csv"]
    out_dir = cfg["out_dir"]
    alpha = cfg.get("fdr_alpha", 0.05)
    df = pd.read_csv(path, index_col=0)
    # Columns like c1.. sort by numeric suffix
    cols_sorted = sorted(df.columns, key=lambda x: int(str(x).split('c')[-1]))
    df = df[cols_sorted]
    # Coerce to numeric and drop non-numeric duplicate header rows if present
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all")
    from tqdm import tqdm
    for label, groups in tqdm(list(cfg["sharma_groups"].items()), desc="Sharma_21 groups"):
        # Select by explicit sample names (e.g., c2, c6) to avoid position mismatches
        exp_names = [f"c{i}" for i in groups["exp"] if f"c{i}" in df.columns]
        ctrl_names = [f"c{i}" for i in groups["ctrl"] if f"c{i}" in df.columns]
        genedf = df[exp_names + ctrl_names]
        exp = genedf[exp_names].values
        ctrl = genedf[ctrl_names].values
        reject, pvals_fdr, fcs = ttest_and_fdr(exp, ctrl)
        # Masks
        up = (pvals_fdr < alpha) & (fcs > cfg.get("fc_threshold", 1.0))
        dw = (pvals_fdr < alpha) & (fcs < 1.0/cfg.get("fc_threshold", 1.0))
        save_up_dw(genedf.index, up, dw, out_dir, f"Sharma_21_{label}")


def run_min(cfg):
    path = cfg["min_csv"]
    out_dir = cfg["out_dir"]
    alpha = cfg.get("fdr_alpha", 0.05)
    df = pd.read_csv(path)
    df.index = df["gene"]
    from tqdm import tqdm
    for item in tqdm(cfg["min_pairs"], desc="Min_18 pairs"):
        label = item["label"]
        regex = item["regex"]
        cols = df.columns[pd.Series(df.columns).str.contains(regex, regex=True)]
        genedf = df[cols].T
        # Split evenly: first half experimental, second half control per original logic
        half = len(genedf)//2
        exp = genedf.iloc[:half, :].values.T
        ctrl = genedf.iloc[half:, :].values.T
        reject, pvals_fdr, fcs = ttest_and_fdr(exp, ctrl)
        up = (pvals_fdr < alpha) & (fcs > cfg.get("fc_threshold", 1.0))
        dw = (pvals_fdr < alpha) & (fcs < 1.0/cfg.get("fc_threshold", 1.0))
        save_up_dw(genedf.columns, up, dw, out_dir, f"min_18_GSE122927_{label}")


def run_qp(cfg):
    path = cfg["qp_csv"]
    out_dir = cfg["out_dir"]
    alpha = cfg.get("fdr_alpha", 0.05)
    df = pd.read_csv(path)
    df.index = df["gene"]
    # Use first ctrl_n columns as control, next exp_n as experimental by default
    ctrl_n = int(cfg.get("qp_split", {}).get("ctrl_n", 3))
    exp_n = int(cfg.get("qp_split", {}).get("exp_n", 3))
    cols = df.columns[1:1+ctrl_n+exp_n]
    genedf = df[cols]
    ctrl = genedf.iloc[:, :ctrl_n].values.T
    exp = genedf.iloc[:, ctrl_n:ctrl_n+exp_n].values.T
    reject, pvals_fdr, fcs = ttest_and_fdr(exp, ctrl)
    up = (pvals_fdr < alpha) & (fcs > cfg.get("fc_threshold", 1.0))
    dw = (pvals_fdr < alpha) & (fcs < 1.0/cfg.get("fc_threshold", 1.0))
    save_up_dw(genedf.index, up, dw, out_dir, "johnson_18_GSE117444")


def run_single_cell(cfg):
    """Identify DEGs for single-cell embryogenesis transitions using findSigGenes.
    Config keys:
      sc_path: path to 10x folder (contains matrix.mtx, genes.tsv, barcodes.tsv)
      sc_label_split: delimiter to parse stage from cell names (default "_")
      sc_label_index: index used after split to extract stage (default 0)
      sc_transitions: list of {label, ref_regex, exp_regex}
      sc_method: "AVGSTD" or "CI"; sc_std_num: default 2
    """
    out_dir = cfg["out_dir"]
    sc_path = cfg.get("sc_path", "")
    if sc_path and not os.path.isabs(sc_path):
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        sc_path = os.path.abspath(os.path.join(repo_root, sc_path))
    if not sc_path:
      return
    label_split = cfg.get("sc_label_split", "_")
    label_index = int(cfg.get("sc_label_index", 0))
    sc_method = cfg.get("sc_method", "AVGSTD")
    sc_std = float(cfg.get("sc_std_num", 2))
    sc_correction = cfg.get("sc_correction", True)

    fr = findSigGenes(sc_path)
    fr.read_scRNAseq()
    genedf = fr.get_genedf(transpose=False)
    # derive labels from cell IDs
    labels = pd.Series(genedf.index).apply(lambda x: str(x).split(label_split)[label_index])

    for tr in cfg.get("sc_transitions", []):
        tlabel = tr.get("label", "sc_transition")
        ref_pat = re.compile(tr.get("ref_regex", ".*"))
        exp_pat = re.compile(tr.get("exp_regex", ".*"))
        ref_cells = labels.apply(lambda s: bool(ref_pat.search(s))).to_numpy()
        exp_cells = labels.apply(lambda s: bool(exp_pat.search(s))).to_numpy()
        print(f"[debug] transition {tlabel}: ref={ref_cells.sum()} exp={exp_cells.sum()} of {len(labels)}")
        print(f"[debug] transition {tlabel}: ref={ref_cells.sum()} exp={exp_cells.sum()} of {len(labels)}")
        updf, dwdf = fr.get_transition_genes(ref_cells, exp_cells, split_str=label_split, method=sc_method, std_num=sc_std, correction=sc_correction, save_files=False)
        # Write per-exp-cell files matching example_sigGenes naming
        for cell in updf.index:
            if 'sc1C2C' in tlabel:
                base = f'1C_to_{cell}'
            elif 'sc2CBC' in tlabel:
                base = f'2C_to_{cell}'
            else:
                base = cell
            up_genes = updf.loc[cell]
            dw_genes = dwdf.loc[cell]
            up_list = up_genes[up_genes == True].index.to_numpy()
            dw_list = dw_genes[dw_genes == True].index.to_numpy()
            pd.DataFrame({"upgenes": up_list}).to_csv(os.path.join(out_dir, f"{base}_upgenes.csv"), index=False)
            pd.DataFrame({"dwgenes": dw_list}).to_csv(os.path.join(out_dir, f"{base}_dwgenes.csv"), index=False)



def extract_first_col(df: pd.DataFrame):
    for c in df.columns:
        if str(c).strip().lower() == 'gene':
            series = df[c]
            return [str(x).strip() for x in series.dropna().astype(str) if str(x).strip()]
    try:
        series = df.iloc[:, -1]
    except Exception:
        series = df.iloc[:, 0]
    return [str(x).strip() for x in series.dropna().astype(str) if str(x).strip()]


def run_scembryo_excels(cfg):
    """Produce per-file upgenes/dwgenes CSVs from scEmbryo Excel files.
    Follows RegulatorsFromDatasets mapping: for names like '1C_to_2cell_...'
      - use the '2C' sheet as upgenes and '1C' as dwgenes (last column or 'gene' column)
      - for '2C_to_32cell_...' use '32C' as upgenes and '2C' as dwgenes.
    """
    indir = Path(cfg.get('sc_xlsx_dir', ''))
    outdir = Path(cfg.get('out_dir', indir))
    if not indir or not indir.exists():
        return
    outdir.mkdir(parents=True, exist_ok=True)
    xfiles = sorted(indir.glob('*.xlsx')); print('[scembryo-excel] files found:', len(xfiles))
    for xf in xfiles:
        try:
            xls = pd.ExcelFile(xf)
        except Exception as e:
            print(f'[scembryo-excel][SKIP] {xf.name}: {e}', file=sys.stderr)
            continue
        sheets = {s.lower(): s for s in xls.sheet_names}
        base = xf.stem
        name = base.lower()
        up_sheet = None
        dw_sheet = None
        if '1c_to_2cell' in name:
            up_sheet = sheets.get('2c')
            dw_sheet = sheets.get('1c')
        elif '2c_to_32cell' in name:
            # allow '32c' or '32cell' in sheet
            up_sheet = sheets.get('32c') or sheets.get('32cell')
            dw_sheet = sheets.get('2c')
        else:
            # Fallback: try any two sheets, map second as up, first as down
            if xls.sheet_names:
                dw_sheet = xls.sheet_names[0]
            if len(xls.sheet_names) > 1:
                up_sheet = xls.sheet_names[1]
        up_list, dw_list = [], []
        if up_sheet:
            try:
                dfu = xls.parse(up_sheet)
                up_list = extract_first_col(dfu)
            except Exception as e:
                print(f'[scembryo-excel][WARN] {xf.name}: up {up_sheet} parse error: {e}', file=sys.stderr)
        if dw_sheet:
            try:
                dfd = xls.parse(dw_sheet)
                dw_list = extract_first_col(dfd)
            except Exception as e:
                print(f'[scembryo-excel][WARN] {xf.name}: dw {dw_sheet} parse error: {e}', file=sys.stderr)
        # Write per-file outputs into out_dir
        pd.DataFrame({'upgenes': up_list}).to_csv(outdir / f'{base}_upgenes.csv', index=False)
        pd.DataFrame({'dwgenes': dw_list}).to_csv(outdir / f'{base}_dwgenes.csv', index=False)
        print(f'[scembryo-excel] {xf.name} -> {base}_upgenes.csv ({len(up_list)}), {base}_dwgenes.csv ({len(dw_list)})')



def main():
    print("[debug] demo_identify_siggenes starting")
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else "identify_siggenes_config.json"
    with open(cfg_path, "r") as f:
        cfg = json.load(f); print("[debug] cfg keys:", list(cfg.keys()))
    # Resolve out_dir relative to repo root if needed
    out_dir = cfg.get("out_dir", "./examples/identifySigGenes_demo/out/")
    if not os.path.isabs(out_dir):
        # Assume run from repo root; adjust if run from elsewhere
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        out_dir = os.path.abspath(os.path.join(repo_root, out_dir))
    os.makedirs(out_dir, exist_ok=True)
    cfg["out_dir"] = out_dir

    # Bulk
    if all(k in cfg for k in ("sharma_csv", "sharma_groups")):
        run_sharma(cfg)
    if all(k in cfg for k in ("min_csv", "min_pairs")):
        run_min(cfg)
    if "qp_csv" in cfg:
        run_qp(cfg)
    # Single-cell (optional)
    if "sc_xlsx_dir" in cfg:
        print("[debug] running excel on", cfg.get("sc_xlsx_dir")); run_scembryo_excels(cfg)
    if "sc_path" in cfg:
        try:
            run_single_cell(cfg)
        except Exception as e:
            print(f"[warn] single-cell processing skipped: {e}")

    print(f"[identifySigGenes_demo] Done. Outputs in: {out_dir}")

if __name__ == "__main__":
    main()
