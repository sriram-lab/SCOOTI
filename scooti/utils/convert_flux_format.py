import argparse
import gzip
import os
import shutil


def convert_file(src_path: str, to_ext: str, keep: bool = False):
    if to_ext not in {"csv", "csv.gz"}:
        raise ValueError("to_ext must be 'csv' or 'csv.gz'")

    if src_path.endswith("_fluxes.csv") and to_ext == "csv.gz":
        dst_path = src_path + ".gz"
        with open(src_path, "rb") as f_in, gzip.open(dst_path, "wb", compresslevel=6) as f_out:
            shutil.copyfileobj(f_in, f_out)
        if not keep:
            os.remove(src_path)
        return dst_path

    if src_path.endswith("_fluxes.csv.gz") and to_ext == "csv":
        dst_path = src_path[:-3]
        with gzip.open(src_path, "rb") as f_in, open(dst_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        if not keep:
            os.remove(src_path)
        return dst_path

    return None


def main():
    ap = argparse.ArgumentParser(description="Convert SCOOTI flux files between .csv and .csv.gz")
    ap.add_argument("--root", required=True, help="Root directory to scan for *_fluxes.csv[.gz]")
    ap.add_argument("--to", choices=["csv", "csv.gz"], required=True, help="Target format")
    ap.add_argument("--keep", action="store_true", help="Keep source files (default deletes after convert)")
    args = ap.parse_args()

    converted = 0
    for dirpath, _dirnames, filenames in os.walk(args.root):
        for fn in filenames:
            if fn.endswith("_fluxes.csv") or fn.endswith("_fluxes.csv.gz"):
                src = os.path.join(dirpath, fn)
                out = convert_file(src, args.to, keep=args.keep)
                if out:
                    converted += 1
                    print(f"Converted: {src} -> {out}")

    print(f"Done. Converted {converted} file(s).")


if __name__ == "__main__":
    main()

