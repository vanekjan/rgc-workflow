#!/usr/bin/env python3
"""
compare_pair_scan_cut_profiles.py

Compare two QADB-enriched qe_pair_scan pair-summary CSV files, typically
cut_profile kin0 versus kin1, using the same broad kin0 skim input.

The merge key is run1, run2, and target_pair_class_name. Outputs:
  <prefix>_pair_comparison.csv
  <prefix>_score_comparison.png/.pdf
  <prefix>_yield_ratio_comparison.png/.pdf
  <prefix>_R2Q_screen_comparison.png/.pdf
  <prefix>_rank_shift_topN.png/.pdf
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from typing import Dict, List, Tuple

try:
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit(
        "Error: matplotlib is required. Try loading a Python environment with matplotlib.\n"
        f"Original import error: {exc}"
    )

Row = Dict[str, str]
Key = Tuple[int, int, str]

VERSION = "compare_pair_scan_cut_profiles_v1_2026_08_13"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare kin0 and kin1 QE pair-scan rankings.")
    p.add_argument("--kin0-csv", required=True, help="QADB-enriched kin0 pair CSV")
    p.add_argument("--kin1-csv", required=True, help="QADB-enriched kin1 pair CSV")
    p.add_argument("--outdir", default="outputs/pair_scan/comparison", help="Output directory")
    p.add_argument("--prefix", default="kin0_vs_kin1", help="Output filename prefix")
    p.add_argument("--top-n", type=int, default=50, help="Top-N pairs for rank-shift bar plot")
    p.add_argument("--pair-class", default="pos_neg", choices=["all", "pos_neg", "pos_pos", "neg_neg"], help="Pair class to compare")
    return p.parse_args()


def read_rows(path: str) -> List[Row]:
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        if "qadb_clean_rank" not in r and "qadb_rank" in r:
            r["qadb_clean_rank"] = r.get("qadb_rank", "")
    return rows


def ival(row: Row, key: str, default: int = 0) -> int:
    try:
        v = row.get(key, "")
        return default if v in {"", None} else int(float(v))
    except Exception:
        return default


def fval(row: Row, key: str, default: float = math.nan) -> float:
    try:
        v = row.get(key, "")
        return default if v in {"", None} else float(v)
    except Exception:
        return default


def sval(row: Row, key: str, default: str = "") -> str:
    v = row.get(key, default)
    return default if v is None else str(v)


def make_key(row: Row) -> Key:
    return (ival(row, "run1"), ival(row, "run2"), sval(row, "target_pair_class_name"))


def index_rows(rows: List[Row], pair_class: str) -> Dict[Key, Row]:
    out: Dict[Key, Row] = {}
    for r in rows:
        if pair_class != "all" and sval(r, "target_pair_class_name") != pair_class:
            continue
        out[make_key(r)] = r
    return out


def savefig(outdir: str, prefix: str, stem: str) -> None:
    png = os.path.join(outdir, f"{prefix}_{stem}.png")
    pdf = os.path.join(outdir, f"{prefix}_{stem}.pdf")
    plt.tight_layout()
    plt.savefig(png, dpi=180)
    plt.savefig(pdf)
    plt.close()
    print(f"wrote {png}")
    print(f"wrote {pdf}")


def scatter(rows: List[Row], outdir: str, prefix: str, xkey: str, ykey: str, xlabel: str, ylabel: str, title: str, stem: str, line_one: bool = False) -> None:
    xs = [fval(r, xkey) for r in rows]
    ys = [fval(r, ykey) for r in rows]
    pts = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if not pts:
        print(f"No finite points for {stem}")
        return
    xs, ys = zip(*pts)
    plt.figure(figsize=(7.2, 6.2))
    plt.scatter(xs, ys, s=18, alpha=0.75)
    if line_one:
        lo = min(min(xs), min(ys))
        hi = max(max(xs), max(ys))
        plt.plot([lo, hi], [lo, hi], linewidth=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.25)
    savefig(outdir, prefix, stem)


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    kin0_rows = read_rows(args.kin0_csv)
    kin1_rows = read_rows(args.kin1_csv)
    i0 = index_rows(kin0_rows, args.pair_class)
    i1 = index_rows(kin1_rows, args.pair_class)

    common = sorted(set(i0).intersection(i1), key=lambda k: (ival(i1[k], "rank", 999999), ival(i0[k], "rank", 999999)))
    out_rows: List[Row] = []

    for key in common:
        a = i0[key]
        b = i1[key]
        run1, run2, cls = key
        r: Row = {
            "run1": str(run1),
            "run2": str(run2),
            "target_pair_class_name": cls,
            "kin0_rank": sval(a, "rank"),
            "kin1_rank": sval(b, "rank"),
            "rank_shift_kin1_minus_kin0": str(ival(b, "rank", 999999) - ival(a, "rank", 999999)),
            "kin0_qadb_category": sval(a, "qadb_pair_category"),
            "kin1_qadb_category": sval(b, "qadb_pair_category"),
            "kin0_score_qadb": sval(a, "pair_quality_score_qadb"),
            "kin1_score_qadb": sval(b, "pair_quality_score_qadb"),
            "score_ratio_kin1_over_kin0": f"{fval(b, 'pair_quality_score_qadb') / fval(a, 'pair_quality_score_qadb'):.6g}" if fval(a, "pair_quality_score_qadb", 0.0) > 0 else "",
            "kin0_yield_ratio_12": sval(a, "yield_ratio_12"),
            "kin1_yield_ratio_12": sval(b, "yield_ratio_12"),
            "delta_yield_ratio_kin1_minus_kin0": f"{fval(b, 'yield_ratio_12') - fval(a, 'yield_ratio_12'):.6g}",
            "kin0_R2Q_screen": sval(a, "R2Q_screen"),
            "kin1_R2Q_screen": sval(b, "R2Q_screen"),
            "delta_R2Q_kin1_minus_kin0": f"{fval(b, 'R2Q_screen') - fval(a, 'R2Q_screen'):.6g}",
            "kin0_abs_deltaQ": sval(a, "abs_deltaQ"),
            "kin1_abs_deltaQ": sval(b, "abs_deltaQ"),
        }
        out_rows.append(r)

    out_csv = os.path.join(args.outdir, f"{args.prefix}_pair_comparison.csv")
    if out_rows:
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
            w.writeheader()
            w.writerows(out_rows)
        print(f"wrote {out_csv}")
    else:
        print("No common pairs found.")
        return

    scatter(out_rows, args.outdir, args.prefix, "kin0_score_qadb", "kin1_score_qadb", "kin0 QADB-adjusted score", "kin1 QADB-adjusted score", "Pair quality score: kin0 vs kin1", "score_comparison", line_one=True)
    scatter(out_rows, args.outdir, args.prefix, "kin0_yield_ratio_12", "kin1_yield_ratio_12", "kin0 Y1/Y2", "kin1 Y1/Y2", "Final yield ratio: kin0 vs kin1", "yield_ratio_comparison", line_one=True)
    scatter(out_rows, args.outdir, args.prefix, "kin0_R2Q_screen", "kin1_R2Q_screen", "kin0 R2Q screen", "kin1 R2Q screen", "Integrated R2Q screen: kin0 vs kin1", "R2Q_screen_comparison", line_one=True)

    top = out_rows[:args.top_n]
    labels = [f"{r['run1']}/{r['run2']}" for r in reversed(top)]
    shifts = [float(r["rank_shift_kin1_minus_kin0"] or 0) for r in reversed(top)]
    height = max(6.0, 0.23 * len(labels) + 2.0)
    plt.figure(figsize=(9.5, height))
    plt.barh(labels, shifts)
    plt.axvline(0, linewidth=1)
    plt.xlabel("Rank shift: kin1 rank - kin0 rank")
    plt.ylabel("Run pair")
    plt.title(f"Top {len(labels)} common pairs: rank shift from kin0 to kin1")
    plt.grid(True, axis="x", alpha=0.25)
    savefig(args.outdir, args.prefix, "rank_shift_topN")

    print("\nSummary")
    print("-------")
    print(f"kin0 pairs:     {len(kin0_rows)}")
    print(f"kin1 pairs:     {len(kin1_rows)}")
    print(f"common pairs:   {len(common)}")
    print(f"output:         {out_csv}")


if __name__ == "__main__":
    main()

