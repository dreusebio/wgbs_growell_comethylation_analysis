#!/usr/bin/env python3
# this script will help rename some of my cytosine report
# Usage
# ./00_rename_growell_files.py \
#   --xlsx /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/metadata/Sample_names_GROWELL.xlsx \
#   --dir /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/data/processed/08_cytosine_reports \
#   --dry-run

# === DRY RUN ===
# Excel: /path/to/your_mapping.xlsx
# Dir:   /path/to/your_files_directory

# 92A_merged_name_sorted.deduplicated.bismark.cov.gz.CpG_report.txt.gz  ->  758053_merged_name_sorted.deduplicated.bismark.cov.gz.CpG_report.txt.gz
# 92A_merged_name_sorted.deduplicated.bismark.cov.gz.cytosine_context_summary.txt  ->  758053_merged_name_sorted.deduplicated.bismark.cov.gz.cytosine_context_summary.txt
# ...

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd


def load_mapping(xlsx_path: Path, sheet: str | None) -> dict[str, str]:
    # Read Excel (requires openpyxl, which you already have in your environment)
    if sheet:
        df = pd.read_excel(xlsx_path, sheet_name=sheet)
    else:
        df = pd.read_excel(xlsx_path)   # automatically reads first sheet

    required = {"Samples_Name_on_Tube", "Unique_Aliquot_ID"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Excel is missing required columns: {sorted(missing)}")

    # Build mapping: "92A" -> "758053"
    m = {}
    for _, row in df.iterrows():
        s = str(row["Samples_Name_on_Tube"]).strip()
        u = str(row["Unique_Aliquot_ID"]).strip()

        # Skip blanks/NaN-like
        if s.lower() in {"nan", ""} or u.lower() in {"nan", ""}:
            continue

        # Force 6-digit formatting if it's numeric-ish
        # (If IDs are already strings, this will keep them)
        try:
            u_int = int(float(u))
            u = f"{u_int:06d}"
        except Exception:
            # if not numeric, leave as-is
            pass

        m[s] = u

    return m


def plan_renames(target_dir: Path, mapping: dict[str, str]) -> list[tuple[Path, Path, str]]:
    """
    Only rename files that start with '<Samples_Name_on_Tube>_'.
    Example:
      92A_merged_name_sorted... -> 758053_merged_name_sorted...
    """
    plans = []
    # Matches: prefix before first underscore
    prefix_re = re.compile(r"^([^_]+)_(.+)$")

    for p in sorted(target_dir.iterdir()):
        if not p.is_file():
            continue

        m = prefix_re.match(p.name)
        if not m:
            continue

        prefix = m.group(1)
        rest = m.group(2)

        if prefix not in mapping:
            continue

        new_prefix = mapping[prefix]
        new_name = f"{new_prefix}_{rest}"
        new_path = p.with_name(new_name)

        action = f"{p.name}  ->  {new_name}"
        plans.append((p, new_path, action))

    return plans


def main():
    ap = argparse.ArgumentParser(
        description="Rename files that start with Samples_Name_on_Tube_ to Unique_Aliquot_ID_."
    )
    ap.add_argument("--xlsx", required=True, help="Path to Excel file containing mapping.")
    ap.add_argument("--sheet", default=None, help="Excel sheet name (optional).")
    ap.add_argument("--dir", required=True, help="Directory containing the files to rename.")
    ap.add_argument("--dry-run", action="store_true", help="Print planned renames; do not rename.")
    ap.add_argument("--log", default="rename_changes.log", help="Log file (default: rename_changes.log).")

    args = ap.parse_args()
    xlsx_path = Path(args.xlsx).expanduser().resolve()
    target_dir = Path(args.dir).expanduser().resolve()
    log_path = Path(args.log).expanduser().resolve()

    if not xlsx_path.exists():
        print(f"ERROR: Excel file not found: {xlsx_path}", file=sys.stderr)
        sys.exit(1)
    if not target_dir.exists():
        print(f"ERROR: Directory not found: {target_dir}", file=sys.stderr)
        sys.exit(1)

    mapping = load_mapping(xlsx_path, args.sheet)
    plans = plan_renames(target_dir, mapping)

    if not plans:
        print("No matching files found to rename (prefix not in mapping or no '<prefix>_' files).")
        sys.exit(0)

    # Collision check
    collisions = [action for oldp, newp, action in plans if newp.exists() and newp != oldp]
    if collisions:
        print("ERROR: The following renames would overwrite existing files:", file=sys.stderr)
        for c in collisions:
            print("  " + c, file=sys.stderr)
        print("Resolve collisions first (move/rename existing target files), then re-run.", file=sys.stderr)
        sys.exit(2)

    # Write log header + planned actions
    mode = "DRY RUN" if args.dry_run else "EXECUTE"
    lines = [f"=== {mode} ===", f"Excel: {xlsx_path}", f"Dir:   {target_dir}", ""]
    lines += [action for _, _, action in plans]
    lines.append("")

    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Print (cat-like) what would change
    print("\n".join(lines))

    if args.dry_run:
        print(f"[dry-run] Wrote plan to: {log_path}")
        return

    # Execute renames
    for oldp, newp, _ in plans:
        os.rename(oldp, newp)

    print(f"[done] Renamed {len(plans)} files.")
    print(f"[done] Log saved to: {log_path}")


if __name__ == "__main__":
    main()