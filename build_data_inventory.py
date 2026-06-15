#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_data_inventory.py -- inventory the DATA that lives on the D: drive (it is
NOT pushed to GitHub, per UKB DUA) and curate only the small AGGREGATE result /
traceability tables into ./results_public/ for the repo.

Why: participant-level data cannot leave the secure environment. The repo ships
a *map* of where every data file sits on D:, plus DUA-safe aggregate results, so
the work is fully findable and the manuscript stays traceable -- without exposing
any patient record.

Outputs (repo root):
    DATA_INVENTORY_ON_D.csv   every data file: D: location, size, rows, DUA class
    DATA_INVENTORY_ON_D.md    human summary (counts + GB by folder + class)
    results_public/           copied AGGREGATE/traceability tables only

Classification rule (robust): a CSV with > 300 data rows, or a non-CSV data file,
or a participant-named file, is PARTICIPANT (never copied). Aggregate result
tables (odds ratios, coefficients, AUC, Table1/2, *_results, *_TRACEABILITY,
*_LEDGER, subgroup, LOCO-CV) with <= 300 rows are AGGREGATE (eligible to copy).
"""
import csv
import os
import re
import shutil

ROOT = os.path.dirname(os.path.abspath(__file__))
INV_ROOTS = [os.path.join(ROOT, d) for d in ("data", "projects", "assets")]
RESULTS_PUBLIC = os.path.join(ROOT, "results_public")

DATA_EXT = {".csv", ".parquet", ".tsv", ".rds", ".dta", ".xlsx", ".xls", ".feather", ".sav"}
RESERVED = {"nul", "con", "prn", "aux"} | {"com%d" % i for i in range(1, 10)} \
    | {"lpt%d" % i for i in range(1, 10)}

# Eligible-to-publish names (aggregate results / traceability):
ALLOW = re.compile(
    r"(odds_ratio|coefficient|_results|TRACEAB|TRACE_LEDGER|_LEDGER|Table\d|"
    r"_summary|LOCO_CV|Subgroup|reproduced_|statistics|metrics|calibration)",
    re.I)
# Hard participant markers (never publish even if small):
DANGER = re.compile(
    r"(analysis_dataset|DRAGON|combined_|ukb_|carrier|_nmr|_prs|master|patient|"
    r"individual|_raw|icd10|gp_|cohort|FULL_MASTER|soretd|baseline_char)", re.I)

ROW_CAP = 302  # read at most this many lines to decide ">300 rows"


def human(n):
    n = float(n)
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return ("%.0f %s" % (n, u)) if u == "B" else ("%.1f %s" % (n, u))
        n /= 1024
    return "%.1f PB" % n


def lines_capped(path, cap=ROW_CAP):
    n = 0
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for _ in f:
                n += 1
                if n >= cap:
                    break
    except OSError:
        return -1
    return n


def classify(path, ext):
    """Return (dua_class, rows_est, hit_cap)."""
    if ext != ".csv":
        return "PARTICIPANT", None, False           # cannot inspect -> conservative
    n = lines_capped(path)
    if n < 0:
        return "PARTICIPANT", None, False
    hit = n >= ROW_CAP
    rows = max(0, n - 1)                              # minus header
    if hit or rows > 300:
        return "PARTICIPANT", (">300" if hit else rows), hit
    return "AGGREGATE", rows, hit


def main():
    os.makedirs(RESULTS_PUBLIC, exist_ok=True)
    rows_out = []
    by_folder = {}
    copied, skipped_big, skipped_danger = 0, 0, 0

    for base in INV_ROOTS:
        if not os.path.isdir(base):
            continue
        top = os.path.basename(base)
        for dirpath, dirnames, filenames in os.walk(base):
            for fn in filenames:
                stem = fn.split(".")[0].lower()
                if stem in RESERVED:
                    continue
                ext = os.path.splitext(fn)[1].lower()
                if ext not in DATA_EXT:
                    continue
                full = os.path.join(dirpath, fn)
                try:
                    sz = os.path.getsize(full)
                except (OSError, ValueError):
                    continue
                dua, rows_est, _ = classify(full, ext)
                rel = os.path.relpath(full, ROOT).replace("\\", "/")

                # curate: only from data/, aggregate, allowlisted, not danger, small
                curated = "no"
                if (top == "data" and dua == "AGGREGATE" and ALLOW.search(fn)
                        and not DANGER.search(fn) and sz <= 1_048_576):
                    try:
                        shutil.copy2(full, os.path.join(RESULTS_PUBLIC, fn))
                        curated = "yes"
                        copied += 1
                    except OSError:
                        pass
                elif top == "data" and ALLOW.search(fn) and dua == "PARTICIPANT":
                    skipped_big += 1
                elif top == "data" and DANGER.search(fn):
                    skipped_danger += 1

                rows_out.append({
                    "folder": top,
                    "rel_path": rel,
                    "abs_path_on_D": full.replace("\\", "/"),
                    "size_bytes": sz,
                    "size_h": human(sz),
                    "rows_est": rows_est if rows_est is not None else "",
                    "dua_class": dua,
                    "in_repo": curated,
                })
                agg = by_folder.setdefault(top, {"n": 0, "bytes": 0, "part": 0, "aggr": 0})
                agg["n"] += 1
                agg["bytes"] += sz
                agg["part"] += 1 if dua == "PARTICIPANT" else 0
                agg["aggr"] += 1 if dua == "AGGREGATE" else 0

    rows_out.sort(key=lambda r: (r["folder"], -r["size_bytes"]))
    inv_csv = os.path.join(ROOT, "DATA_INVENTORY_ON_D.csv")
    with open(inv_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "folder", "rel_path", "abs_path_on_D", "size_bytes", "size_h",
            "rows_est", "dua_class", "in_repo"])
        w.writeheader()
        w.writerows(rows_out)

    total_n = len(rows_out)
    total_b = sum(r["size_bytes"] for r in rows_out)
    n_part = sum(1 for r in rows_out if r["dua_class"] == "PARTICIPANT")
    n_aggr = total_n - n_part
    with open(os.path.join(ROOT, "DATA_INVENTORY_ON_D.md"), "w", encoding="utf-8") as fh:
        fh.write("# Data inventory -- on the D: drive (NOT in this repo)\n\n")
        fh.write("UK Biobank / Welsh-registry participant data stays on D: under the "
                 "UKB DUA. This is the map of where it lives; the repo ships only the "
                 "aggregate result tables in `results_public/`.\n\n")
        fh.write("| Metric | Value |\n|---|---|\n")
        fh.write("| Data files inventoried | **%d** |\n" % total_n)
        fh.write("| Total size on D: | **%s** |\n" % human(total_b))
        fh.write("| PARTICIPANT (restricted, on D: only) | **%d** |\n" % n_part)
        fh.write("| AGGREGATE (safe) | **%d** |\n" % n_aggr)
        fh.write("| Copied into results_public/ | **%d** |\n\n" % copied)
        fh.write("## By folder\n\n| Folder | Files | Size | Participant | Aggregate |\n")
        fh.write("|---|---|---|---|---|\n")
        for k in sorted(by_folder):
            a = by_folder[k]
            fh.write("| `%s/` | %d | %s | %d | %d |\n"
                     % (k, a["n"], human(a["bytes"]), a["part"], a["aggr"]))
        fh.write("\nFull per-file map: `DATA_INVENTORY_ON_D.csv` "
                 "(column `abs_path_on_D` = exact location on the D: drive).\n")

    print("[inventory] %d data files, %s on D:  (PARTICIPANT=%d, AGGREGATE=%d)"
          % (total_n, human(total_b), n_part, n_aggr))
    print("[inventory] results_public/: copied %d aggregate tables; "
          "skipped %d big + %d participant-named" % (copied, skipped_big, skipped_danger))
    print("[inventory] wrote DATA_INVENTORY_ON_D.csv + .md")


if __name__ == "__main__":
    main()
