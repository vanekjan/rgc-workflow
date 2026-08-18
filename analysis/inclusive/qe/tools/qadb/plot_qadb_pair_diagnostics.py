#!/usr/bin/env python3
"""
plot_qadb_pair_diagnostics.py

Make diagnostic plots from a QADB-enriched qe_pair_scan pair-summary CSV.

Typical use from analysis/qe:

  python3 tools/qadb/plot_qadb_pair_diagnostics.py \
    --pair-csv outputs/qe_pair_scan_spring23_kin0_pair_summary_qadb.csv \
    --outdir outputs/pair_diagnostics_kin0 \
    --prefix kin0

The script writes PNG/PDF plots and a sorted candidate table:

  <prefix>_top_pairs_qadb.csv
  <prefix>_yield_ratio_vs_abs_deltaQ.png/.pdf
  <prefix>_R2Q_vs_abs_deltaQ.png/.pdf
  <prefix>_raw_vs_final_yield_ratio.png/.pdf
  <prefix>_theta_survival_vs_final_yield_ratio.png/.pdf
  <prefix>_score_vs_abs_deltaQ.png/.pdf
  <prefix>_top50_qadb_score.png/.pdf
  <prefix>_qadb_pair_category_counts.png/.pdf

No ROOT is required. It uses Python stdlib plus matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Optional, Tuple

try:
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "Error: matplotlib is required for plotting. Try loading a Python environment with matplotlib.\n"
        f"Original import error: {exc}"
    )

Row = Dict[str, str]

PLOT_QADB_PAIR_DIAGNOSTICS_VERSION = "plot_qadb_pair_diagnostics_v2_qadb_clean_rank_2026_08_12"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot QADB-enriched QE pair-scan diagnostics."
    )
    p.add_argument("--pair-csv", required=True, help="QADB-enriched pair summary CSV")
    p.add_argument("--outdir", default="outputs/pair_diagnostics", help="Output directory")
    p.add_argument("--prefix", default="pairscan", help="Output filename prefix")
    p.add_argument("--top-n", type=int, default=50, help="Number of top pairs for ranked bar plot and top CSV")
    p.add_argument("--include-rejected", type=int, default=0, choices=[0, 1], help="Include QADB-rejected pairs in top-pair ranking plots/tables")
    p.add_argument("--pair-class", default="all", choices=["all", "pos_neg", "pos_pos", "neg_neg"], help="Restrict ranked/top plots to a target pair class")
    p.add_argument("--min-abs-deltaq", type=float, default=0.0, help="Minimum |Delta Q| for ranked/top plots")
    p.add_argument("--max-abs-r2q-plot", type=float, default=5.0, help="Maximum |R2Q_screen| shown in R2Q scatter plot")
    p.add_argument("--yield-ratio-min", type=float, default=0.0, help="Optional minimum yield_ratio_12 for ranked/top plots. 0 disables.")
    p.add_argument("--yield-ratio-max", type=float, default=0.0, help="Optional maximum yield_ratio_12 for ranked/top plots. 0 disables.")
    return p.parse_args()


def read_csv(path: str) -> List[Row]:
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    # Backward compatibility: older QADB postprocessor versions used qadb_rank.
    # The current name is qadb_clean_rank.
    for r in rows:
        if "qadb_clean_rank" not in r and "qadb_rank" in r:
            r["qadb_clean_rank"] = r.get("qadb_rank", "")
    return rows


def fval(row: Row, key: str, default: float = math.nan) -> float:
    try:
        s = row.get(key, "")
        if s is None or s == "":
            return default
        return float(s)
    except Exception:
        return default


def ival(row: Row, key: str, default: int = 0) -> int:
    try:
        s = row.get(key, "")
        if s is None or s == "":
            return default
        return int(float(s))
    except Exception:
        return default


def sval(row: Row, key: str, default: str = "") -> str:
    s = row.get(key, default)
    return default if s is None else str(s)


def is_finite(x: float) -> bool:
    return math.isfinite(x)


def is_rejected(row: Row) -> bool:
    rec = sval(row, "qadb_pair_recommendation").lower()
    cat = sval(row, "qadb_pair_category").lower()
    score = fval(row, "pair_quality_score_qadb", 0.0)
    pen = fval(row, "qadb_pair_penalty", 1.0)
    return rec.startswith("reject") or cat in {"junk_or_tuning", "empty_or_special", "helicity_warning"} or score <= 0.0 or pen <= 0.0


def pair_label(row: Row) -> str:
    return f"{ival(row, 'run1')}/{ival(row, 'run2')}"


def sort_rows_by_qadb_score(rows: Iterable[Row]) -> List[Row]:
    return sorted(
        rows,
        key=lambda r: (
            fval(r, "pair_quality_score_qadb", -1.0),
            fval(r, "pair_quality_score", -1.0),
            fval(r, "abs_deltaQ", -1.0),
        ),
        reverse=True,
    )


def filtered_for_top(rows: List[Row], args: argparse.Namespace) -> List[Row]:
    out: List[Row] = []
    for r in rows:
        if not args.include_rejected and is_rejected(r):
            continue
        if args.pair_class != "all" and sval(r, "target_pair_class_name") != args.pair_class:
            continue
        if fval(r, "abs_deltaQ", -1.0) < args.min_abs_deltaq:
            continue
        yr = fval(r, "yield_ratio_12")
        if args.yield_ratio_min > 0.0 and (not is_finite(yr) or yr < args.yield_ratio_min):
            continue
        if args.yield_ratio_max > 0.0 and (not is_finite(yr) or yr > args.yield_ratio_max):
            continue
        out.append(r)
    return sort_rows_by_qadb_score(out)


def savefig(outdir: str, prefix: str, stem: str) -> None:
    png = os.path.join(outdir, f"{prefix}_{stem}.png")
    pdf = os.path.join(outdir, f"{prefix}_{stem}.pdf")
    plt.tight_layout()
    plt.savefig(png, dpi=180)
    plt.savefig(pdf)
    plt.close()
    print(f"wrote {png}")
    print(f"wrote {pdf}")


def grouped_scatter(
    rows: List[Row],
    outdir: str,
    prefix: str,
    stem: str,
    title: str,
    xkey: str,
    ykey: str,
    xlabel: str,
    ylabel: str,
    category_key: str = "qadb_pair_category",
    yline: Optional[float] = None,
    xline: Optional[float] = None,
    ylim: Optional[Tuple[float, float]] = None,
) -> None:
    groups: Dict[str, List[Tuple[float, float]]] = defaultdict(list)
    for r in rows:
        x = fval(r, xkey)
        y = fval(r, ykey)
        if not (is_finite(x) and is_finite(y)):
            continue
        if ylim is not None and (y < ylim[0] or y > ylim[1]):
            continue
        cat = sval(r, category_key, "unknown") or "unknown"
        groups[cat].append((x, y))

    plt.figure(figsize=(8.5, 6.2))
    for cat in sorted(groups.keys()):
        pts = groups[cat]
        if not pts:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        plt.scatter(xs, ys, s=18, alpha=0.75, label=cat)
    if yline is not None:
        plt.axhline(yline, linewidth=1)
    if xline is not None:
        plt.axvline(xline, linewidth=1)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8, loc="best")
    savefig(outdir, prefix, stem)


def write_top_csv(rows: List[Row], outdir: str, prefix: str, top_n: int) -> str:
    keep_cols = [
        "rank", "qadb_clean_rank", "pair_scan_rank",
        "run1", "run2", "target_pair_class_name",
        "P1", "P2", "Q1", "Q2", "abs_deltaQ",
        "N1", "N2", "Y1", "Y2", "yield_ratio_12", "R2Q_screen",
        "raw_yield_ratio_12", "inclusive_yield_ratio_12", "theta_survival_ratio_12",
        "dominant_trigger_bit_all_1", "dominant_trigger_bit_all_2",
        "dominant_trigger_bit_selected_1", "dominant_trigger_bit_selected_2",
        "same_dominant_trigger_bit_all", "same_dominant_trigger_bit_selected",
        "pair_quality_score", "qadb_pair_category", "qadb_pair_penalty",
        "qadb_pair_recommendation", "pair_quality_score_qadb",
        "qadb1_qadb_category", "qadb1_misc_comments",
        "qadb2_qadb_category", "qadb2_misc_comments",
    ]
    path = os.path.join(outdir, f"{prefix}_top_pairs_qadb.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keep_cols, extrasaction="ignore")
        w.writeheader()
        for r in rows[:top_n]:
            w.writerow(r)
    print(f"wrote {path}")
    return path


def plot_top_scores(rows: List[Row], outdir: str, prefix: str, top_n: int) -> None:
    top = rows[:top_n]
    if not top:
        print("No rows available for top score plot after filtering.")
        return
    # Reverse so best is at top in a horizontal bar chart.
    plot_rows = list(reversed(top))
    labels = [pair_label(r) for r in plot_rows]
    scores = [fval(r, "pair_quality_score_qadb", 0.0) for r in plot_rows]

    height = max(6.0, 0.22 * len(plot_rows) + 2.0)
    plt.figure(figsize=(9.0, height))
    plt.barh(labels, scores)
    plt.xlabel("QADB-adjusted pair quality score")
    plt.ylabel("Run pair")
    plt.title(f"Top {len(plot_rows)} run pairs after QADB adjustment")
    plt.grid(True, axis="x", alpha=0.25)
    savefig(outdir, prefix, "top50_qadb_score")


def plot_category_counts(rows: List[Row], outdir: str, prefix: str) -> None:
    counts = Counter(sval(r, "qadb_pair_category", "unknown") or "unknown" for r in rows)
    items = sorted(counts.items(), key=lambda kv: kv[1])
    labels = [k for k, _ in items]
    values = [v for _, v in items]

    plt.figure(figsize=(8.5, 5.5))
    plt.barh(labels, values)
    plt.xlabel("Number of run pairs")
    plt.ylabel("QADB pair category")
    plt.title("Run-pair counts by QADB category")
    plt.grid(True, axis="x", alpha=0.25)
    savefig(outdir, prefix, "qadb_pair_category_counts")


def print_summary(rows: List[Row], top_rows: List[Row], args: argparse.Namespace) -> None:
    print("\nSummary")
    print("-------")
    print(f"input pairs: {len(rows)}")
    print(f"top-filter pairs: {len(top_rows)}")
    print("QADB pair categories:")
    for cat, n in sorted(Counter(sval(r, "qadb_pair_category", "unknown") or "unknown" for r in rows).items()):
        print(f"  {cat:22s} {n}")

    print("\nTop pairs after QADB adjustment:")
    cols = [
        "rank", "qadb_clean_rank", "pair_scan_rank",
        "run1", "run2", "target_pair_class_name", "abs_deltaQ",
        "yield_ratio_12", "R2Q_screen", "raw_yield_ratio_12",
        "pair_quality_score", "qadb_pair_category", "pair_quality_score_qadb",
    ]
    print(",".join(cols))
    for r in top_rows[:min(args.top_n, 20)]:
        print(",".join(sval(r, c) for c in cols))


def main() -> None:
    args = parse_args()
    rows = read_csv(args.pair_csv)
    os.makedirs(args.outdir, exist_ok=True)

    top_rows = filtered_for_top(rows, args)
    write_top_csv(top_rows, args.outdir, args.prefix, args.top_n)

    # All-pair diagnostic scatter plots. These keep rejected pairs visible so the
    # user can see how QADB categories populate the metric space.
    grouped_scatter(
        rows, args.outdir, args.prefix,
        "yield_ratio_vs_abs_deltaQ",
        "Final yield ratio vs tensor lever arm",
        "abs_deltaQ", "yield_ratio_12",
        "|Q1 - Q2|", "Y1/Y2",
        yline=1.0,
    )

    grouped_scatter(
        rows, args.outdir, args.prefix,
        "R2Q_vs_abs_deltaQ",
        "Integrated R2Q screen vs tensor lever arm",
        "abs_deltaQ", "R2Q_screen",
        "|Q1 - Q2|", "R2Q screen",
        yline=0.0,
        ylim=(-abs(args.max_abs_r2q_plot), abs(args.max_abs_r2q_plot)),
    )

    grouped_scatter(
        rows, args.outdir, args.prefix,
        "raw_vs_final_yield_ratio",
        "Raw yield ratio vs final yield ratio",
        "raw_yield_ratio_12", "yield_ratio_12",
        "Raw yield ratio", "Final yield ratio Y1/Y2",
        yline=1.0,
        xline=1.0,
    )

    grouped_scatter(
        rows, args.outdir, args.prefix,
        "theta_survival_vs_final_yield_ratio",
        "Theta-survival ratio vs final yield ratio",
        "theta_survival_ratio_12", "yield_ratio_12",
        "Theta-survival ratio", "Final yield ratio Y1/Y2",
        yline=1.0,
        xline=1.0,
    )

    grouped_scatter(
        rows, args.outdir, args.prefix,
        "score_vs_abs_deltaQ",
        "QADB-adjusted score vs tensor lever arm",
        "abs_deltaQ", "pair_quality_score_qadb",
        "|Q1 - Q2|", "QADB-adjusted pair quality score",
    )

    plot_top_scores(top_rows, args.outdir, args.prefix, args.top_n)
    plot_category_counts(rows, args.outdir, args.prefix)
    print_summary(rows, top_rows, args)


if __name__ == "__main__":
    main()
