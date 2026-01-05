#!/usr/bin/env python3
import argparse, sys
from pathlib import Path
import pandas as pd
import re

def extract_first_col(df: pd.DataFrame):
    # Prefer a column named 'gene' if present
    for c in df.columns:
        if str(c).strip().lower() == 'gene':
            series = df[c]
            return [str(x).strip() for x in series.dropna().astype(str) if str(x).strip()]
    # If there are at least 2 columns, the last column often contains gene symbols
    try:
        series = df.iloc[:, -1]
    except Exception:
        series = df.iloc[:, 0]
    return [str(x).strip() for x in series.dropna().astype(str) if str(x).strip()]


def main():
    ap = argparse.ArgumentParser(description='Per-file conversion: XLSX with upgenes/dwgenes sheets -> *_upgenes.csv and *_dwgenes.csv')
    ap.add_argument('indir', type=Path, help='Directory with .xlsx files')
    ap.add_argument('--outdir', type=Path, default=None, help='Output directory (default: same as input)')
    args = ap.parse_args()

    indir: Path = args.indir
    outdir: Path = args.outdir or indir
    outdir.mkdir(parents=True, exist_ok=True)

    xfiles = sorted(indir.glob('*.xlsx'))
    if not xfiles:
        print(f'[convert] No .xlsx files under {indir}', file=sys.stderr)
        sys.exit(2)


    for xf in xfiles:
        try:
            xls = pd.ExcelFile(xf)
        except Exception as e:
            print(f'[convert][SKIP] Cannot open {xf}: {e}', file=sys.stderr)
            continue
        sheets = {s.lower(): s for s in xls.sheet_names}
        # Convert either explicit upgenes/dwgenes sheets OR infer ref/exp sheets from filename
        base = xf.stem
        up_list, dw_list = [], []
        if 'upgenes' in sheets or 'dwgenes' in sheets:
            if 'upgenes' in sheets:
                try:
                    dfu = xls.parse(sheets['upgenes'])
                    up_list = extract_first_col(dfu)
                except Exception as e:
                    print(f'[convert][WARN] {xf.name}: upgenes sheet parse error: {e}', file=sys.stderr)
            if 'dwgenes' in sheets:
                try:
                    dfd = xls.parse(sheets['dwgenes'])
                    dw_list = extract_first_col(dfd)
                except Exception as e:
                    print(f'[convert][WARN] {xf.name}: dwgenes sheet parse error: {e}', file=sys.stderr)
        else:
            # Infer sheets from filename like '1C_to_2cell_*'
            name = base.lower()
            m = re.search(r'(\w+)_to_(\w+)', name)
            ref_tok = exp_tok = None
            if m:
                ref_tok, exp_tok = m.group(1), m.group(2)
            def norm(s: str) -> str:
                s = s.lower()
                s = re.sub(r'[^a-z0-9]', '', s)
                s = s.replace('cell', 'c')
                return s
            find = {norm(k): v for k, v in sheets.items()}
            ref_sheet = find.get(norm(ref_tok)) if ref_tok else None
            exp_sheet = find.get(norm(exp_tok)) if exp_tok else None
            if not exp_sheet and exp_tok:
                for nk, v in find.items():
                    if nk.startswith(norm(exp_tok)[:2]):
                        exp_sheet = v; break
            if not ref_sheet and ref_tok:
                for nk, v in find.items():
                    if nk.startswith(norm(ref_tok)[:2]):
                        ref_sheet = v; break
            if exp_sheet:
                try:
                    up_list = extract_first_col(xls.parse(exp_sheet))
                except Exception as e:
                    print(f'[convert][WARN] {xf.name}: exp sheet parse error: {e}', file=sys.stderr)
            if ref_sheet:
                try:
                    dw_list = extract_first_col(xls.parse(ref_sheet))
                except Exception as e:
                    print(f'[convert][WARN] {xf.name}: ref sheet parse error: {e}', file=sys.stderr)
        # Write outputs
        up_path = outdir / f'{base}_upgenes.csv'
        dw_path = outdir / f'{base}_dwgenes.csv'
        pd.DataFrame({'upgenes': up_list}).to_csv(up_path, index=False)
        pd.DataFrame({'dwgenes': dw_list}).to_csv(dw_path, index=False)
        print(f"[convert] {xf.name} -> {up_path.name} ({len(up_list)}), {dw_path.name} ({len(dw_list)})")

if __name__ == '__main__':
    main()
