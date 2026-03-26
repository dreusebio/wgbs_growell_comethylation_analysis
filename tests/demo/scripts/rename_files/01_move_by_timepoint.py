#!/usr/bin/env python3
# Usage
# ./01_move_by_timepoint.py   --xlsx /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/metadata/Sample_names_GROWELL.xlsx   --dir /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/data/processed/08_cytosine_reports

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd


def load_id_to_timepoint(xlsx_path: Path, sheet: str | None) -> dict[str, str]:
    if sheet:
        df = pd.read_excel(xlsx_path, sheet_name=sheet)
    else:
        df = pd.read_excel(xlsx_path)   # automatically reads first sheet

    required = {"Unique_Aliquot_ID", "POD_Timepoint"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Excel is missing required columns: {sorted(missing)}")

    id2tp: dict[str, str] = {}

    for _, row in df.iterrows():
        uid = str(row["Unique_Aliquot_ID"]).strip()
        tp = str(row["POD_Timepoint"]).strip()

        if uid.lower() in {"nan", ""} or tp.lower() in {"nan", ""}:
            continue

        # normalize ID to 6 digits when numeric
        try:
            uid_int = int(float(uid))
            uid = f"{uid_int:06d}"
        except Exception:
            pass

        id2tp[uid] = tp

    return id2tp


def safe_dirname(name: str) -> str:
    # Make a filesystem-safe directory name
    # Keep letters/numbers/underscore/dash; map spaces and slashes to underscores
    name = name.strip()
    name = re.sub(r"[ /]+", "_", name)
    name = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    return name


def plan_moves(target_dir: Path, id2tp: dict[str, str], out_root: Path) -> list[tuple[Path, Path, str]]:
    """
    Move files with prefix '<Unique_Aliquot_ID>_' into '<out_root>/<Timepoint>/'.
    """
    plans = []
    prefix_re = re.compile(r"^(\d{6})_(.+)$")

    for p in sorted(target_dir.iterdir()):
        if not p.is_file():
            continue

        m = prefix_re.match(p.name)
        if not m:
            continue

        uid = m.group(1)
        if uid not in id2tp:
            continue

        tp = safe_dirname(id2tp[uid])
        dest_dir = out_root / tp
        dest_path = dest_dir / p.name

        action = f"mkdir -p {dest_dir}\nmove {p}  ->  {dest_path}"
        plans.append((p, dest_path, action))

    return plans


def main():
    ap = argparse.ArgumentParser(
        description="Move files into timepoint directories using Unique_Aliquot_ID -> POD_Timepoint from Excel."
    )
    ap.add_argument("--xlsx", required=True, help="Path to Excel file containing mapping.")
    ap.add_argument("--sheet", default=None, help="Excel sheet name (optional).")
    ap.add_argument("--dir", required=True, help="Directory containing the files to move.")
    ap.add_argument("--out-root", default=None,
                    help="Root directory where timepoint folders will be created (default: same as --dir).")
    ap.add_argument("--dry-run", action="store_true", help="Print planned moves; do not move anything.")
    ap.add_argument("--log", default="move_by_timepoint.log", help="Log file (default: move_by_timepoint.log).")
    ap.add_argument("--copy", action="store_true", help="Copy instead of move (default is move).")

    args = ap.parse_args()
    xlsx_path = Path(args.xlsx).expanduser().resolve()
    target_dir = Path(args.dir).expanduser().resolve()
    out_root = Path(args.out_root).expanduser().resolve() if args.out_root else target_dir
    log_path = Path(args.log).expanduser().resolve()

    if not xlsx_path.exists():
        print(f"ERROR: Excel file not found: {xlsx_path}", file=sys.stderr)
        sys.exit(1)
    if not target_dir.exists():
        print(f"ERROR: Directory not found: {target_dir}", file=sys.stderr)
        sys.exit(1)
    if not out_root.exists():
        print(f"ERROR: out-root not found: {out_root}", file=sys.stderr)
        sys.exit(1)

    id2tp = load_id_to_timepoint(xlsx_path, args.sheet)
    plans = plan_moves(target_dir, id2tp, out_root)

    if not plans:
        print("No files matched the pattern '<6-digit-ID>_' with an ID found in the Excel mapping.")
        sys.exit(0)

    # Collision check (destination exists)
    collisions = [act for src, dst, act in plans if dst.exists() and dst != src]
    if collisions:
        print("ERROR: Some destination files already exist (would overwrite):", file=sys.stderr)
        for c in collisions[:50]:
            print("  " + c.replace("\n", " | "), file=sys.stderr)
        if len(collisions) > 50:
            print(f"  ... and {len(collisions)-50} more", file=sys.stderr)
        print("Resolve collisions first (move/rename existing dest files), then re-run.", file=sys.stderr)
        sys.exit(2)

    mode = "DRY RUN" if args.dry_run else ("COPY" if args.copy else "MOVE")
    header = [
        f"=== {mode} BY TIMEPOINT ===",
        f"Excel:    {xlsx_path}",
        f"From dir: {target_dir}",
        f"Out root: {out_root}",
        "",
    ]

    # Print a readable plan
    actions = []
    for _, _, act in plans:
        actions.append(act)
        actions.append("")  # blank line between items

    text = "\n".join(header + actions).rstrip() + "\n"
    log_path.write_text(text, encoding="utf-8")

    # "cat" plan to stdout
    print(text)

    if args.dry_run:
        print(f"[dry-run] Plan written to: {log_path}")
        return

    # Execute
    for src, dst, _ in plans:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if args.copy:
            # copy then keep original
            import shutil
            shutil.copy2(src, dst)
        else:
            os.rename(src, dst)

    print(f"[done] {mode.lower()}d {len(plans)} files.")
    print(f"[done] Log saved to: {log_path}")


if __name__ == "__main__":
    main()