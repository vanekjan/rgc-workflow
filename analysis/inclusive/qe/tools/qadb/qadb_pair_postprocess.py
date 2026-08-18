#!/usr/bin/env python3
"""
qadb_pair_postprocess.py

Post-process CLAS12 QADB metadata for RGC Spring 2023 and merge it into
qe_pair_scan pair-summary CSV files.

Primary inputs:
  - miscTable.md
  - qaTree.json.table
  - qe_pair_scan_*_pair_summary.csv

Outputs:
  - qadb_run_metadata.csv
  - qe_pair_scan_*_pair_summary_qadb.csv

This script intentionally works as a postprocessor: it does not modify the skim
or the two-run analysis output. It adds run-level and pair-level QA metadata for
pair selection/ranking.
"""

from __future__ import annotations

import argparse
import csv
import html
import os
import re
import sys
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

MISC_URL = "https://raw.githubusercontent.com/JeffersonLab/clas12-qadb/main/qadb/pass1/rgc_sp23/miscTable.md"
QA_TABLE_URL = "https://raw.githubusercontent.com/JeffersonLab/clas12-qadb/main/qadb/pass1/rgc_sp23/qaTree.json.table"

POSTPROCESS_VERSION = "qadb_pair_postprocess_v2_ranked_2026_08_12"

DEFECT_NAMES = {
    0: "TotalOutlier",
    1: "TerminalOutlier",
    2: "MarginalOutlier",
    3: "SectorLoss",
    4: "LowLiveTime",
    5: "Misc",
    6: "TotalOutlierFT",
    7: "TerminalOutlierFT",
    8: "MarginalOutlierFT",
    9: "LossFT",
    10: "BSAWrong",
    11: "BSAUnknown",
    12: "TSAWrong",
    13: "TSAUnknown",
    14: "DSAWrong",
    15: "DSAUnknown",
    16: "ChargeHigh",
    17: "ChargeNegative",
    18: "ChargeUnknown",
    19: "PossiblyNoBeam",
}

# Approximate severity scale for ranking pair metadata. This is intentionally
# conservative; it should guide review, not silently veto physics choices.
CATEGORY_SEVERITY = {
    "clean": 0,
    "cnd_warning": 1,
    "qadb_warning": 2,
    "target_pol_warning": 3,
    "helicity_warning": 4,
    "junk_or_tuning": 5,
    "empty_or_special": 6,
}

CATEGORY_PENALTY = {
    "clean": 1.00,
    "cnd_warning": 0.90,
    "qadb_warning": 0.70,
    "target_pol_warning": 0.30,
    "helicity_warning": 0.20,
    "junk_or_tuning": 0.00,
    "empty_or_special": 0.00,
}


@dataclass
class RunQA:
    run: int
    misc_defect: str = "unknown"
    misc_comments: str = ""
    n_qadb_bins: int = 0
    n_bins_with_any_defect: int = 0
    defect_bin_counts: Dict[str, int] = field(default_factory=dict)
    qa_comments: List[str] = field(default_factory=list)

    # classified comment/defect flags
    has_misc: bool = False
    has_target_pol_issue: bool = False
    has_cnd_issue: bool = False
    has_raster_issue: bool = False
    has_junk_or_tuning_issue: bool = False
    has_empty_or_special_issue: bool = False
    has_helicity_issue: bool = False
    has_detector_issue: bool = False
    has_sector_nf_issue: bool = False
    has_charge_issue: bool = False
    has_low_livetime_issue: bool = False
    has_no_beam_issue: bool = False
    category: str = "clean"

    def count(self, defect_name: str) -> int:
        return int(self.defect_bin_counts.get(defect_name, 0))

    def frac(self, defect_name: str) -> float:
        if self.n_qadb_bins <= 0:
            return 0.0
        return self.count(defect_name) / float(self.n_qadb_bins)

    @property
    def any_defect_fraction(self) -> float:
        if self.n_qadb_bins <= 0:
            return 0.0
        return self.n_bins_with_any_defect / float(self.n_qadb_bins)


def read_text(path_or_url: str, cache_dir: Optional[Path] = None) -> str:
    """Read a local file or URL. URLs are cached if cache_dir is provided."""
    if re.match(r"^https?://", path_or_url):
        if cache_dir is not None:
            cache_dir.mkdir(parents=True, exist_ok=True)
            name = path_or_url.rstrip("/").split("/")[-1]
            cached = cache_dir / name
            if cached.exists() and cached.stat().st_size > 0:
                return cached.read_text(errors="replace")
            with urllib.request.urlopen(path_or_url, timeout=60) as resp:
                data = resp.read().decode("utf-8", errors="replace")
            cached.write_text(data)
            return data
        with urllib.request.urlopen(path_or_url, timeout=60) as resp:
            return resp.read().decode("utf-8", errors="replace")
    return Path(path_or_url).read_text(errors="replace")


def clean_comment(text: str) -> str:
    text = html.unescape(text)
    text = re.sub(r"</?pre>", "", text)
    text = re.sub(r"<[^>]+>", "", text)
    text = text.replace("\r", " ").replace("\n", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def parse_misc_table(text: str) -> Dict[int, RunQA]:
    """Parse qadb/pass1/rgc_sp23/miscTable.md."""
    runs: Dict[int, RunQA] = {}
    row_re = re.compile(r"^\|\s*(\d+)\s*\|\s*(yes|no)\s*\|\s*(.*?)\s*\|\s*$", re.I)
    for line in text.splitlines():
        m = row_re.match(line)
        if not m:
            continue
        run = int(m.group(1))
        misc = m.group(2).lower()
        comments = clean_comment(m.group(3))
        runs[run] = RunQA(run=run, misc_defect=misc, misc_comments=comments)
    return runs


def parse_qa_table(text: str, runs: Optional[Dict[int, RunQA]] = None) -> Dict[int, RunQA]:
    """Parse qadb/pass1/rgc_sp23/qaTree.json.table.

    Example table line:
      17482 0 :: empty target :: 1-TerminalOutlier[all] 5-Misc[all]
    """
    if runs is None:
        runs = {}

    # Defect token: 5-Misc[all], 0-TotalOutlier[4], etc.
    defect_re = re.compile(r"(\d+)\-([A-Za-z0-9_]+)\[([^\]]+)\]")
    line_re = re.compile(r"^\s*(\d{4,6})\s+(\d+)\s+::\s*(.*?)\s*::\s*(.*)$")

    for line in text.splitlines():
        m = line_re.match(line)
        if not m:
            continue
        run = int(m.group(1))
        # bin_number = int(m.group(2))
        comment = clean_comment(m.group(3))
        defect_text = m.group(4)
        obj = runs.setdefault(run, RunQA(run=run))
        obj.n_qadb_bins += 1
        if comment and comment not in obj.qa_comments:
            obj.qa_comments.append(comment)

        defects = defect_re.findall(defect_text)
        if defects:
            obj.n_bins_with_any_defect += 1
        for bit_s, name, _scope in defects:
            bit = int(bit_s)
            canonical = DEFECT_NAMES.get(bit, name)
            obj.defect_bin_counts[canonical] = obj.defect_bin_counts.get(canonical, 0) + 1

    return runs


def classify_run(obj: RunQA, bad_bin_warning_threshold: float = 0.25) -> None:
    combined_comments = "; ".join(x for x in [obj.misc_comments] + obj.qa_comments if x).strip()
    low = combined_comments.lower()

    obj.has_misc = (obj.misc_defect == "yes") or obj.count("Misc") > 0
    obj.has_target_pol_issue = any(s in low for s in [
        "target polarization", "polarization dropped", "polarization issue", "pol. issue",
    ])
    obj.has_cnd_issue = "cnd" in low
    obj.has_raster_issue = "raster" in low
    obj.has_junk_or_tuning_issue = any(s in low for s in ["junk", "tuning", "beam tuning", "raster tuning"])
    obj.has_empty_or_special_issue = any(s in low for s in [
        "empty target", "low luminosity", "special run", "lh2", "ld2", "carbon", "dummy",
    ])
    obj.has_helicity_issue = (
        "helicity" in low or obj.count("BSAWrong") > 0 or obj.count("BSAUnknown") > 0
    )
    obj.has_detector_issue = any(s in low for s in [
        "cnd", "dc", "ft", "fd", "sector", "pmt", "tdc", "tof", "calorimeter", "ecal",
    ])
    obj.has_sector_nf_issue = ("sector" in low and "n/f" in low)
    obj.has_charge_issue = obj.count("ChargeHigh") > 0 or obj.count("ChargeNegative") > 0 or obj.count("ChargeUnknown") > 0
    obj.has_low_livetime_issue = obj.count("LowLiveTime") > 0
    obj.has_no_beam_issue = obj.count("PossiblyNoBeam") > 0

    if obj.has_empty_or_special_issue:
        obj.category = "empty_or_special"
    elif obj.has_junk_or_tuning_issue:
        obj.category = "junk_or_tuning"
    elif obj.has_target_pol_issue:
        obj.category = "target_pol_warning"
    elif obj.has_helicity_issue:
        obj.category = "helicity_warning"
    elif obj.has_cnd_issue:
        obj.category = "cnd_warning"
    elif obj.has_misc or obj.any_defect_fraction >= bad_bin_warning_threshold:
        obj.category = "qadb_warning"
    else:
        obj.category = "clean"


def classify_all(runs: Dict[int, RunQA], bad_bin_warning_threshold: float) -> Dict[int, RunQA]:
    for obj in runs.values():
        classify_run(obj, bad_bin_warning_threshold=bad_bin_warning_threshold)
    return runs


def run_to_row(obj: RunQA) -> Dict[str, object]:
    comments = obj.misc_comments or "; ".join(obj.qa_comments)
    row: Dict[str, object] = {
        "run": obj.run,
        "qadb_category": obj.category,
        "qadb_penalty": CATEGORY_PENALTY.get(obj.category, 0.0),
        "misc_defect": obj.misc_defect,
        "misc_comments": comments,
        "n_qadb_bins": obj.n_qadb_bins,
        "n_bins_with_any_defect": obj.n_bins_with_any_defect,
        "any_defect_fraction": f"{obj.any_defect_fraction:.6g}",
        "has_misc": int(obj.has_misc),
        "has_target_pol_issue": int(obj.has_target_pol_issue),
        "has_cnd_issue": int(obj.has_cnd_issue),
        "has_raster_issue": int(obj.has_raster_issue),
        "has_junk_or_tuning_issue": int(obj.has_junk_or_tuning_issue),
        "has_empty_or_special_issue": int(obj.has_empty_or_special_issue),
        "has_helicity_issue": int(obj.has_helicity_issue),
        "has_detector_issue": int(obj.has_detector_issue),
        "has_sector_nf_issue": int(obj.has_sector_nf_issue),
        "has_charge_issue": int(obj.has_charge_issue),
        "has_low_livetime_issue": int(obj.has_low_livetime_issue),
        "has_no_beam_issue": int(obj.has_no_beam_issue),
    }
    for name in DEFECT_NAMES.values():
        key = f"n_{name}_bins"
        row[key] = obj.count(name)
        row[f"frac_{name}_bins"] = f"{obj.frac(name):.6g}"
    return row


def write_run_csv(path: Path, runs: Dict[int, RunQA]) -> None:
    rows = [run_to_row(runs[k]) for k in sorted(runs)]
    if not rows:
        raise RuntimeError("No QADB run metadata rows were parsed")
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def prefixed_run_row(prefix: str, obj: Optional[RunQA]) -> Dict[str, object]:
    if obj is None:
        base = {
            "qadb_category": "missing",
            "qadb_penalty": 0.0,
            "misc_defect": "missing",
            "misc_comments": "missing from parsed QADB tables",
            "n_qadb_bins": 0,
            "any_defect_fraction": "0",
            "has_misc": 0,
            "has_target_pol_issue": 0,
            "has_cnd_issue": 0,
            "has_raster_issue": 0,
            "has_junk_or_tuning_issue": 0,
            "has_empty_or_special_issue": 0,
            "has_helicity_issue": 0,
            "has_detector_issue": 0,
            "has_sector_nf_issue": 0,
            "has_charge_issue": 0,
            "has_low_livetime_issue": 0,
            "has_no_beam_issue": 0,
        }
    else:
        base = run_to_row(obj)
        base.pop("run", None)
    return {f"{prefix}{k}": v for k, v in base.items()}


def pair_category(r1: Optional[RunQA], r2: Optional[RunQA]) -> str:
    cats = [(r1.category if r1 else "missing"), (r2.category if r2 else "missing")]
    if "missing" in cats:
        return "missing_qadb"
    return max(cats, key=lambda c: CATEGORY_SEVERITY.get(c, 99))


def recommendation_for_category(cat: str) -> str:
    if cat == "clean":
        return "keep_candidate"
    if cat == "cnd_warning":
        return "ok_for_inclusive_electron_review_for_hadron_or_cnd"
    if cat == "qadb_warning":
        return "review_qadb_bins"
    if cat == "target_pol_warning":
        return "avoid_for_final_or_split_by_time"
    if cat == "helicity_warning":
        return "review_or_fix_helicity_before_spin_use"
    if cat in {"junk_or_tuning", "empty_or_special"}:
        return "reject_for_final_physics_pair"
    return "review_missing_qadb"


def _as_float(value: object, default: float = 0.0) -> float:
    try:
        if value is None:
            return default
        text = str(value).strip()
        if text == "":
            return default
        return float(text)
    except Exception:
        return default


def _as_int(value: object, default: int = 0) -> int:
    try:
        if value is None:
            return default
        text = str(value).strip()
        if text == "":
            return default
        return int(float(text))
    except Exception:
        return default


def _qadb_only_sort_key(row: Dict[str, object]) -> Tuple[float, int, int, int, int, int, int, int, float, int]:
    """Sort key for QADB-only rank. Lower tuple is better."""
    penalty = _as_float(row.get("qadb_pair_penalty"), 0.0)
    cat = str(row.get("qadb_pair_category", "missing_qadb"))
    severity = CATEGORY_SEVERITY.get(cat, 99)
    max_defect_fraction = max(
        _as_float(row.get("qadb1_any_defect_fraction"), 0.0),
        _as_float(row.get("qadb2_any_defect_fraction"), 0.0),
    )
    return (
        -penalty,
        severity,
        _as_int(row.get("qadb_pair_has_junk_or_tuning_issue")),
        _as_int(row.get("qadb_pair_has_empty_or_special_issue")),
        _as_int(row.get("qadb_pair_has_helicity_issue")),
        _as_int(row.get("qadb_pair_has_target_pol_issue")),
        _as_int(row.get("qadb_pair_has_detector_issue")),
        _as_int(row.get("qadb_pair_has_misc")),
        max_defect_fraction,
        _as_int(row.get("pair_scan_rank"), 999999),
    )


def _main_sort_key(row: Dict[str, object]) -> Tuple[float, float, float, int]:
    """Sort key for the final/main rank. Lower tuple is better."""
    return (
        -_as_float(row.get("pair_quality_score_qadb"), 0.0),
        -_as_float(row.get("qadb_pair_penalty"), 0.0),
        -_as_float(row.get("pair_quality_score"), 0.0),
        _as_int(row.get("pair_scan_rank"), 999999),
    )


def merge_pair_csv(pair_csv: Path, out_csv: Path, runs: Dict[int, RunQA]) -> None:
    """Merge run-level QADB metadata into a qe_pair_scan pair CSV.

    Output ranking convention:
      rank           = final/main rank after QADB penalty, sorted by
                       pair_quality_score_qadb descending
      qadb_clean_rank      = QADB-only metadata rank; clean pairs rank ahead of
                       warning/rejected pairs independent of physics score
      pair_scan_rank = original rank from qe_pair_scan, or input row order
                       if the input CSV has no rank column

    The output CSV is automatically sorted by the final/main rank.
    """
    with pair_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        input_fieldnames = list(reader.fieldnames or [])

    if not rows:
        raise RuntimeError(f"Pair CSV has no rows: {pair_csv}")
    if "run1" not in input_fieldnames or "run2" not in input_fieldnames:
        raise RuntimeError("Pair CSV must contain run1 and run2 columns")

    # Preserve the scanner's original rank but free the name 'rank' for the
    # final QADB-adjusted rank. This makes the first column the one users
    # should normally sort/read by.
    data_fieldnames = [x for x in input_fieldnames if x != "rank"]

    new_rows: List[Dict[str, object]] = []
    added_fields: List[str] = []

    for input_index, row in enumerate(rows, start=1):
        run1 = int(float(row["run1"]))
        run2 = int(float(row["run2"]))
        qa1 = runs.get(run1)
        qa2 = runs.get(run2)
        pcat = pair_category(qa1, qa2)
        ppenalty = CATEGORY_PENALTY.get(pcat, 0.0)

        enriched: Dict[str, object] = {k: v for k, v in row.items() if k != "rank"}
        enriched["pair_scan_rank"] = row.get("rank", input_index)

        extras: Dict[str, object] = {}
        extras.update(prefixed_run_row("qadb1_", qa1))
        extras.update(prefixed_run_row("qadb2_", qa2))
        extras["qadb_pair_category"] = pcat
        extras["qadb_pair_penalty"] = ppenalty
        extras["qadb_pair_recommendation"] = recommendation_for_category(pcat)

        # Boolean pair flags from either run.
        flag_names = [
            "has_misc", "has_target_pol_issue", "has_cnd_issue", "has_raster_issue",
            "has_junk_or_tuning_issue", "has_empty_or_special_issue", "has_helicity_issue",
            "has_detector_issue", "has_sector_nf_issue", "has_charge_issue",
            "has_low_livetime_issue", "has_no_beam_issue",
        ]
        for flag in flag_names:
            v1 = int(getattr(qa1, flag, False)) if qa1 else 0
            v2 = int(getattr(qa2, flag, False)) if qa2 else 0
            extras[f"qadb_pair_{flag}"] = int(bool(v1 or v2))

        # Penalized score for convenience if the scanner score exists.
        score_val = None
        for score_key in ["pair_quality_score", "quality_score", "score"]:
            if score_key in row and row[score_key] not in {"", None}:
                try:
                    score_val = float(row[score_key])
                    break
                except ValueError:
                    pass
        extras["pair_quality_score_qadb"] = f"{(score_val * ppenalty):.6g}" if score_val is not None else ""

        for k, v in extras.items():
            enriched[k] = v
            if k not in added_fields:
                added_fields.append(k)
        new_rows.append(enriched)

    # Assign QADB-only rank before the final score sort. This keeps a separate
    # rank that answers: "which pairs are cleanest by QADB metadata alone?"
    qadb_sorted_indices = sorted(range(len(new_rows)), key=lambda i: _qadb_only_sort_key(new_rows[i]))
    for qadb_clean_rank, idx in enumerate(qadb_sorted_indices, start=1):
        new_rows[idx]["qadb_clean_rank"] = qadb_clean_rank

    # Main output order: QADB-adjusted score rank.
    new_rows.sort(key=_main_sort_key)
    for final_rank, row in enumerate(new_rows, start=1):
        row["rank"] = final_rank

    # Put the high-level ranks first, then the original scanner columns, then
    # the QADB metadata block.
    front_fields = ["rank", "qadb_clean_rank", "pair_scan_rank"]
    remaining_data_fields = [x for x in data_fieldnames if x not in front_fields]
    remaining_added_fields = [x for x in added_fields if x not in front_fields]
    fieldnames = front_fields + remaining_data_fields + remaining_added_fields

    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(new_rows)

def summarize(runs: Dict[int, RunQA], pair_csv: Optional[Path], out_pair_csv: Optional[Path]) -> None:
    counts: Dict[str, int] = {}
    for r in runs.values():
        counts[r.category] = counts.get(r.category, 0) + 1
    print("QADB run metadata summary:")
    for cat in sorted(counts, key=lambda c: CATEGORY_SEVERITY.get(c, 99)):
        print(f"  {cat:20s} {counts[cat]:5d}")
    if pair_csv and out_pair_csv:
        print(f"\nMerged pair CSV:")
        print(f"  input : {pair_csv}")
        print(f"  output: {out_pair_csv}")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Parse RGC Spring 2023 QADB metadata and merge it into qe_pair_scan pair-summary CSVs."
    )
    p.add_argument("--misc-table", default=None,
                   help="Path or URL to miscTable.md. Default: RGC Spring 2023 GitHub raw table.")
    p.add_argument("--qa-table", default=None,
                   help="Path or URL to qaTree.json.table. Default: RGC Spring 2023 GitHub raw table.")
    p.add_argument("--qadb-dir", default=None,
                   help="Local clas12-qadb directory. If supplied, uses qadb/pass1/rgc_sp23 tables inside it.")
    p.add_argument("--cache-dir", default="qadb_cache",
                   help="Directory for downloaded QADB files. Default: qadb_cache")
    p.add_argument("--pair-csv", default=None,
                   help="qe_pair_scan pair summary CSV to enrich.")
    p.add_argument("--out-run-csv", default="outputs/qadb_rgc_sp23_run_metadata.csv",
                   help="Output run metadata CSV.")
    p.add_argument("--out-pair-csv", default=None,
                   help="Output enriched pair CSV. Default: input name with _qadb before .csv")
    p.add_argument("--bad-bin-warning-threshold", type=float, default=0.25,
                   help="If any-defect bin fraction exceeds this value, category becomes qadb_warning unless a stronger category applies. Default 0.25")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    qadb_dir = Path(args.qadb_dir).expanduser() if args.qadb_dir else None
    if qadb_dir:
        default_misc = qadb_dir / "qadb" / "pass1" / "rgc_sp23" / "miscTable.md"
        default_qa = qadb_dir / "qadb" / "pass1" / "rgc_sp23" / "qaTree.json.table"
    else:
        default_misc = Path(args.misc_table) if args.misc_table else None
        default_qa = Path(args.qa_table) if args.qa_table else None

    misc_source = args.misc_table or (str(default_misc) if default_misc else MISC_URL)
    qa_source = args.qa_table or (str(default_qa) if default_qa else QA_TABLE_URL)
    cache_dir = Path(args.cache_dir) if args.cache_dir else None

    out_run_csv = Path(args.out_run_csv)
    out_run_csv.parent.mkdir(parents=True, exist_ok=True)

    print(f"QADB pair postprocessor version: {POSTPROCESS_VERSION}")
    print("Reading QADB inputs:")
    print(f"  misc table: {misc_source}")
    print(f"  qa table:   {qa_source}")

    misc_text = read_text(misc_source, cache_dir=cache_dir)
    qa_text = read_text(qa_source, cache_dir=cache_dir)

    runs = parse_misc_table(misc_text)
    runs = parse_qa_table(qa_text, runs=runs)
    classify_all(runs, bad_bin_warning_threshold=args.bad_bin_warning_threshold)

    write_run_csv(out_run_csv, runs)
    print(f"\nWrote run metadata CSV:\n  {out_run_csv}")

    pair_csv = Path(args.pair_csv) if args.pair_csv else None
    out_pair_csv = None
    if pair_csv:
        if args.out_pair_csv:
            out_pair_csv = Path(args.out_pair_csv)
        else:
            out_pair_csv = pair_csv.with_name(pair_csv.stem + "_qadb.csv")
        out_pair_csv.parent.mkdir(parents=True, exist_ok=True)
        merge_pair_csv(pair_csv, out_pair_csv, runs)
        print(f"  sorted by QADB-adjusted pair rank: {out_pair_csv}")

    summarize(runs, pair_csv, out_pair_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
