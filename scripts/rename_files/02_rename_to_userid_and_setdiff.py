#!/usr/bin/env python3

# Usage


# ./02_rename_to_userid_and_setdiff.py \
#   --xlsx /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/metadata/Sample_names_GROWELL.xlsx \
#   --root /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/data/processed/08_cytosine_reports \
#   --folders Baseline 36-38wk Postpartum \
#   --expected-n 49 \
#   --dry-run \

#   --xlsx /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/metadata/Sample_names_GROWELL.xlsx \
#   --dir /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/data/processed/08_cytosine_reports \
#   --dry-run

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd


def safe_dirname(name: str) -> str:
    name = str(name).strip()
    name = re.sub(r"[ /]+", "_", name)
    name = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    return name


def load_uid_to_userid(xlsx_path: Path, sheet: str | None) -> dict[str, str]:
    if sheet:
        df = pd.read_excel(xlsx_path, sheet_name=sheet)
    else:
        df = pd.read_excel(xlsx_path)   # automatically reads first sheet

    required = {"Unique_Aliquot_ID", "User_ID"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Excel is missing required columns: {sorted(missing)}")

    m: dict[str, str] = {}

    for _, row in df.iterrows():
        uid = str(row["Unique_Aliquot_ID"]).strip()
        user_id = str(row["User_ID"]).strip()

        if uid.lower() in {"nan", ""} or user_id.lower() in {"nan", ""}:
            continue

        # normalize aliquot ID to 6 digits when numeric
        try:
            uid_int = int(float(uid))
            uid = f"{uid_int:06d}"
        except Exception:
            pass

        # normalize user_id to integer-like string when numeric (e.g., "92.0" -> "92")
        try:
            user_int = int(float(user_id))
            user_id = str(user_int)
        except Exception:
            user_id = user_id.strip()

        m[uid] = user_id

    return m


def list_timepoint_dirs(root: Path, names: list[str] | None) -> list[Path]:
    if names:
        return [root / safe_dirname(n) for n in names]

    # Default: all subdirs under root
    return sorted([p for p in root.iterdir() if p.is_dir()])


def sample_ids_in_dir(tp_dir: Path, mode: str) -> set[str]:
    """
    mode:
      - 'aliquot' expects filenames like 758053_...
      - 'user'   expects filenames like 92_...
    """
    ids: set[str] = set()
    if mode == "aliquot":
        pat = re.compile(r"^(\d{6})_")
    elif mode == "user":
        pat = re.compile(r"^([^_]+)_")
    else:
        raise ValueError("mode must be 'aliquot' or 'user'")

    for p in tp_dir.iterdir():
        if not p.is_file():
            continue
        m = pat.match(p.name)
        if m:
            ids.add(m.group(1))
    return ids


def plan_renames_in_dir(tp_dir: Path, uid2user: dict[str, str]) -> list[tuple[Path, Path, str]]:
    plans: list[tuple[Path, Path, str]] = []
    pat = re.compile(r"^(\d{6})_(.+)$")

    for p in sorted(tp_dir.iterdir()):
        if not p.is_file():
            continue
        m = pat.match(p.name)
        if not m:
            continue

        uid = m.group(1)
        rest = m.group(2)

        if uid not in uid2user:
            continue

        new_prefix = uid2user[uid]
        new_name = f"{new_prefix}_{rest}"
        new_path = p.with_name(new_name)
        action = f"{tp_dir.name}: {p.name}  ->  {new_name}"
        plans.append((p, new_path, action))

    return plans


def collision_check(plans: list[tuple[Path, Path, str]]) -> list[str]:
    collisions = []
    for src, dst, action in plans:
        if dst.exists() and dst != src:
            collisions.append(action)
    return collisions


def setdiff_report(tp_dirs: list[Path], expected_n: int) -> str:
    # collect per-folder sample User_IDs (assumes already renamed to user prefix)
    per = {}
    for d in tp_dirs:
        if not d.exists():
            per[d.name] = set()
            continue
        per[d.name] = sample_ids_in_dir(d, mode="user")

    folder_names = [d.name for d in tp_dirs]
    all_union = set().union(*per.values()) if per else set()

    lines = []
    lines.append("=== SET DIFF REPORT (by User_ID) ===")
    lines.append("Folders: " + ", ".join(folder_names))
    lines.append(f"Union size across all folders: {len(all_union)}")
    lines.append("")

    # counts + missing-from-each vs union
    for name in folder_names:
        s = per.get(name, set())
        lines.append(f"[{name}] unique samples: {len(s)} (expected {expected_n})")
        missing = sorted(all_union - s, key=lambda x: (len(x), x))
        if missing:
            lines.append(f"  Missing vs union ({len(missing)}): " + ", ".join(missing))
        else:
            lines.append("  Missing vs union: none")
        lines.append("")

    # pairwise diffs
    lines.append("=== PAIRWISE DIFFERENCES ===")
    for i in range(len(folder_names)):
        for j in range(i + 1, len(folder_names)):
            a = folder_names[i]
            b = folder_names[j]
            A = per[a]
            B = per[b]
            only_a = sorted(A - B, key=lambda x: (len(x), x))
            only_b = sorted(B - A, key=lambda x: (len(x), x))
            lines.append(f"{a} \\ {b} (in {a} not {b}): {len(only_a)}")
            if only_a:
                lines.append("  " + ", ".join(only_a))
            lines.append(f"{b} \\ {a} (in {b} not {a}): {len(only_b)}")
            if only_b:
                lines.append("  " + ", ".join(only_b))
            lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def main():
    ap = argparse.ArgumentParser(
        description="Rename 6-digit Unique_Aliquot_ID prefixes to User_ID within timepoint folders, then compute set diffs."
    )
    ap.add_argument("--xlsx", required=True, help="Excel file with Unique_Aliquot_ID, User_ID, POD_Timepoint.")
    ap.add_argument("--sheet", default=None, help="Excel sheet name (optional).")
    ap.add_argument("--root", required=True, help="Root folder containing timepoint directories (Baseline, 36-38wk, Postpartum).")
    ap.add_argument("--folders", nargs="*", default=None,
                    help="Optional explicit folder names (e.g., Baseline 36-38wk Postpartum). If omitted, all subdirs are used.")
    ap.add_argument("--expected-n", type=int, default=49, help="Expected unique sample count per folder (default: 49).")
    ap.add_argument("--dry-run", action="store_true", help="Print planned renames + setdiff; do not rename.")
    ap.add_argument("--log", default="rename_to_userid_and_setdiff.log", help="Log file (default: rename_to_userid_and_setdiff.log).")

    args = ap.parse_args()
    xlsx_path = Path(args.xlsx).expanduser().resolve()
    root = Path(args.root).expanduser().resolve()
    log_path = Path(args.log).expanduser().resolve()

    if not xlsx_path.exists():
        print(f"ERROR: Excel not found: {xlsx_path}", file=sys.stderr)
        sys.exit(1)
    if not root.exists():
        print(f"ERROR: Root folder not found: {root}", file=sys.stderr)
        sys.exit(1)

    uid2user = load_uid_to_userid(xlsx_path, args.sheet)
    tp_dirs = list_timepoint_dirs(root, args.folders)

    # Build rename plans for each folder
    all_plans: list[tuple[Path, Path, str]] = []
    missing_dirs = []
    for d in tp_dirs:
        if not d.exists():
            missing_dirs.append(d)
            continue
        all_plans.extend(plan_renames_in_dir(d, uid2user))

    mode = "DRY RUN" if args.dry_run else "EXECUTE"

    header = []
    header.append(f"=== RENAME Unique_Aliquot_ID -> User_ID ({mode}) ===")
    header.append(f"Excel: {xlsx_path}")
    header.append(f"Root:  {root}")
    header.append("Folders considered: " + ", ".join([d.name for d in tp_dirs]))
    if missing_dirs:
        header.append("WARNING missing folders: " + ", ".join([d.name for d in missing_dirs]))
    header.append("")

    # Collision check before doing anything
    collisions = collision_check(all_plans)
    if collisions:
        msg = ["ERROR: collisions detected (would overwrite existing files):"]
        msg += ["  " + c for c in collisions[:100]]
        if len(collisions) > 100:
            msg.append(f"  ... and {len(collisions)-100} more")
        msg.append("Resolve collisions (rename/move existing dest files) then re-run.")
        text = "\n".join(header + msg) + "\n"
        log_path.write_text(text, encoding="utf-8")
        print(text, end="")
        sys.exit(2)

    # Render planned rename actions
    rename_lines = header[:]
    if not all_plans:
        rename_lines.append("No files found matching '^\\d{6}_' that also exist in mapping.")
        rename_lines.append("")
    else:
        rename_lines.append(f"Planned renames: {len(all_plans)} files")
        rename_lines.append("")
        for _, _, action in all_plans:
            rename_lines.append(action)
        rename_lines.append("")

    # If executing, do the renames now
    if not args.dry_run:
        for src, dst, _ in all_plans:
            os.rename(src, dst)
        rename_lines.append(f"[done] Renamed {len(all_plans)} files.")
        rename_lines.append("")

    # After (or during dry-run), compute setdiff report
    # Assumes filenames are now User_ID_... (or will be after rename).
    # In dry-run, we compute using current state; that's usually OK if you already renamed,
    # but if you haven't renamed yet, setdiff will reflect aliquot-prefixes.
    # So: in dry-run we ALSO compute expected per-folder IDs from the planned operations.
    if args.dry_run:
        # Predict per-folder user IDs after rename by applying plans virtually
        per_after = {}
        for d in tp_dirs:
            if not d.exists():
                per_after[d.name] = set()
                continue
            # start with existing user IDs
            existing_user_ids = sample_ids_in_dir(d, mode="user")
            # add user IDs that would be created from aliquot-prefix renames
            planned_for_dir = [pl for pl in all_plans if pl[0].parent == d]
            planned_user_ids = set()
            for src, _, _ in planned_for_dir:
                m = re.match(r"^(\d{6})_", src.name)
                if m:
                    uid = m.group(1)
                    if uid in uid2user:
                        planned_user_ids.add(uid2user[uid])
            # after rename, aliquot-prefixed files become user-prefixed
            # union is safe approximation for set diff
            per_after[d.name] = existing_user_ids.union(planned_user_ids)

        # build setdiff text from per_after
        folder_names = [d.name for d in tp_dirs]
        all_union = set().union(*per_after.values()) if per_after else set()

        lines = []
        lines.append("=== SET DIFF REPORT (predicted after rename; by User_ID) ===")
        lines.append("Folders: " + ", ".join(folder_names))
        lines.append(f"Union size across all folders: {len(all_union)}")
        lines.append("")
        for name in folder_names:
            s = per_after.get(name, set())
            lines.append(f"[{name}] unique samples: {len(s)} (expected {args.expected_n})")
            missing = sorted(all_union - s, key=lambda x: (len(x), x))
            if missing:
                lines.append(f"  Missing vs union ({len(missing)}): " + ", ".join(missing))
            else:
                lines.append("  Missing vs union: none")
            lines.append("")
        lines.append("=== PAIRWISE DIFFERENCES ===")
        for i in range(len(folder_names)):
            for j in range(i + 1, len(folder_names)):
                a = folder_names[i]
                b = folder_names[j]
                A = per_after[a]
                B = per_after[b]
                only_a = sorted(A - B, key=lambda x: (len(x), x))
                only_b = sorted(B - A, key=lambda x: (len(x), x))
                lines.append(f"{a} \\ {b} (in {a} not {b}): {len(only_a)}")
                if only_a:
                    lines.append("  " + ", ".join(only_a))
                lines.append(f"{b} \\ {a} (in {b} not {a}): {len(only_b)}")
                if only_b:
                    lines.append("  " + ", ".join(only_b))
                lines.append("")
        setdiff_text = "\n".join(lines).rstrip() + "\n"
    else:
        setdiff_text = setdiff_report(tp_dirs, args.expected_n)

    out_text = "\n".join(rename_lines).rstrip() + "\n\n" + setdiff_text
    log_path.write_text(out_text, encoding="utf-8")
    print(out_text, end="")
    print(f"[log] {log_path}")


if __name__ == "__main__":
    main()