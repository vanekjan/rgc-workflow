#!/usr/bin/env python3

"""
submit_slurm_skim.py

Submit rgcskim jobs to the JLab Slurm farm as an array job.

Design:
  one Slurm array task = one HIPO file

This script:
  1. reads the normal configs/skim.yaml file,
  2. resolves the HIPO input files,
  3. applies run selection and target-map filtering at file-list level,
  4. writes a file list under logs/slurm_filelists/,
  5. submits scripts/batch/rgcskim_array.sbatch with one task per file.

The batch worker still calls scripts/run_all_hipo.sh for each file, so the
normal rgcskim output naming and ROOT/log folder organization are preserved.
"""

import argparse
import csv
import datetime as dt
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys


def load_yaml(path):
    try:
        import yaml
    except ImportError:
        print("Error: Python module 'yaml' was not found.")
        print("Try:")
        print("  python3 -m pip install --user PyYAML")
        sys.exit(1)

    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    if data is None:
        print(f"Error: YAML file is empty: {path}")
        sys.exit(1)

    if not isinstance(data, dict):
        print(f"Error: YAML top level must be a dictionary: {path}")
        sys.exit(1)

    return data


def get_section(cfg, name):
    value = cfg.get(name, {})
    if value is None:
        return {}
    if not isinstance(value, dict):
        print(f"Error: YAML section '{name}' must be a dictionary.")
        sys.exit(1)
    return value


def as_string(value):
    if value is None:
        return ""
    return str(value).strip()


def normalize_lower(value):
    return as_string(value).lower()


def sanitize_tag(raw):
    raw = as_string(raw)
    if raw == "":
        return "slurm"
    return re.sub(r"[^A-Za-z0-9._-]", "_", raw)


def repo_relative_or_absolute(repo_root, raw_path):
    raw_path = as_string(raw_path)

    if raw_path == "":
        return raw_path

    path = Path(os.path.expanduser(raw_path))

    if path.is_absolute():
        return str(path)

    return str((repo_root / path).resolve())


def path_for_command(repo_root, raw_path):
    """
    Keep relative paths relative for run_all_hipo.sh, because batch jobs cd into
    the repo before running. Absolute paths stay absolute.
    """
    raw_path = as_string(raw_path)

    if raw_path == "":
        return ""

    path = Path(os.path.expanduser(raw_path))

    if path.is_absolute():
        return str(path)

    return raw_path


def extract_run_number(path):
    base = Path(path).stem
    matches = re.findall(r"[0-9]{6}", base)

    if not matches:
        return None

    return int(matches[-1])


def parse_runs_value(value):
    runs = []

    if value is None:
        return runs

    if isinstance(value, list):
        items = value
    else:
        items = re.split(r"[,\s]+", str(value).strip())

    for item in items:
        item = as_string(item)

        if item == "":
            continue

        if not re.fullmatch(r"[0-9]+", item):
            print(f"Error: invalid run number in list: {item}")
            sys.exit(1)

        runs.append(int(item))

    return runs


def read_runlist(path):
    runs = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()

            if line == "":
                continue

            parts = re.split(r"[,\s]+", line)

            for part in parts:
                if part == "":
                    continue

                if not re.fullmatch(r"[0-9]+", part):
                    print(f"Error: invalid run number in runlist {path}: {part}")
                    sys.exit(1)

                runs.append(int(part))

    return runs


def read_target_map_runs(path):
    runs = set()

    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)

        for row in reader:
            if not row:
                continue

            found_run = None

            for item in row:
                item = item.strip()
                match = re.search(r"\b[0-9]{5,6}\b", item)

                if not match:
                    continue

                candidate = int(match.group(0))

                if 10000 <= candidate <= 999999:
                    found_run = candidate
                    break

            if found_run is not None:
                runs.add(found_run)

    return runs


def collect_hipo_files(input_path):
    path = Path(input_path)

    if path.is_file():
        if path.suffix != ".hipo":
            print(f"Error: input file is not a .hipo file: {path}")
            sys.exit(1)
        return [str(path)]

    if path.is_dir():
        files = sorted(str(p) for p in path.glob("*.hipo"))

        if not files:
            print(f"Error: no .hipo files found in input directory: {path}")
            sys.exit(1)

        return files

    print(f"Error: input.path does not exist: {path}")
    sys.exit(1)


def filter_by_target_map(files, target_map_runs):
    kept = []
    skipped_no_run = 0
    skipped_no_target = 0

    for path in files:
        run = extract_run_number(path)

        if run is None:
            skipped_no_run += 1
            continue

        if run not in target_map_runs:
            skipped_no_target += 1
            continue

        kept.append(path)

    return kept, skipped_no_run, skipped_no_target


def filter_by_run_selection(files, run_mode, runsel, repo_root):
    run_mode = normalize_lower(run_mode)

    if run_mode in ["", "none", "all", "test_one_run"]:
        kept = files
    elif run_mode in ["single", "run"]:
        run = as_string(runsel.get("run", ""))

        if run == "":
            print("Error: run_selection.mode single requires run_selection.run.")
            sys.exit(1)

        selected = {int(run)}
        kept = [p for p in files if extract_run_number(p) in selected]

    elif run_mode in ["list", "runs"]:
        selected = set(parse_runs_value(runsel.get("runs", [])))

        if not selected:
            print("Error: run_selection.mode list requires run_selection.runs.")
            sys.exit(1)

        kept = [p for p in files if extract_run_number(p) in selected]

    elif run_mode in ["range", "runrange", "run_range"]:
        start = as_string(runsel.get("start", ""))
        end = as_string(runsel.get("end", ""))

        if start == "" or end == "":
            runrange = as_string(runsel.get("runrange", ""))

            if runrange != "":
                parts = re.split(r"[:-]", runrange)

                if len(parts) == 2:
                    start = parts[0].strip()
                    end = parts[1].strip()

        if start == "" or end == "":
            print("Error: run_selection.mode range requires start/end or runrange.")
            sys.exit(1)

        start_i = int(start)
        end_i = int(end)

        if end_i < start_i:
            print("Error: run_selection.end is smaller than run_selection.start.")
            sys.exit(1)

        kept = []

        for path in files:
            run = extract_run_number(path)

            if run is not None and start_i <= run <= end_i:
                kept.append(path)

    elif run_mode in ["file", "runlist"]:
        runlist = as_string(runsel.get("runlist", ""))

        if runlist == "":
            print("Error: run_selection.mode file requires run_selection.runlist.")
            sys.exit(1)

        runlist_path = repo_relative_or_absolute(repo_root, runlist)

        if not os.path.isfile(runlist_path):
            print(f"Error: runlist does not exist: {runlist_path}")
            sys.exit(1)

        selected = set(read_runlist(runlist_path))
        kept = [p for p in files if extract_run_number(p) in selected]

    else:
        print(f"Error: unsupported run_selection.mode: {run_mode}")
        print("Allowed: all, test_one_run, single, range, list, file")
        sys.exit(1)

    if run_mode == "test_one_run":
        kept = kept[:1]

    return kept


def build_argument_parser():
    parser = argparse.ArgumentParser(
        description="Submit rgcskim YAML configuration as a Slurm array job."
    )

    parser.add_argument(
        "config",
        help="YAML config, usually configs/skim.yaml",
    )

    parser.add_argument(
        "--account",
        default="hallc",
        help="Slurm account. Default: hallc",
    )

    parser.add_argument(
        "--partition",
        default="production",
        help="Slurm partition. Default: production",
    )

    parser.add_argument(
        "--time",
        default="04:00:00",
        help="Slurm wall time. Default: 04:00:00",
    )

    parser.add_argument(
        "--mem",
        default="2G",
        help="Slurm memory per task. Default: 2G",
    )

    parser.add_argument(
        "--cpus-per-task",
        default="1",
        help="Slurm CPUs per task. Default: 1",
    )

    parser.add_argument(
        "--array-limit",
        type=int,
        default=20,
        help="Maximum simultaneous array tasks. Default: 20. Use 0 for no limit.",
    )

    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Limit number of HIPO files submitted. Useful for testing.",
    )

    parser.add_argument(
        "--jobtag",
        default="",
        help="Override output.jobtag from YAML.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Create file list and print sbatch command, but do not submit.",
    )

    return parser


def main():
    parser = build_argument_parser()
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[2]

    config_path = Path(os.path.expanduser(args.config))

    if not config_path.is_absolute():
        config_path = repo_root / config_path

    config_path = config_path.resolve()

    if not config_path.is_file():
        print(f"Error: config file does not exist: {config_path}")
        sys.exit(1)

    cfg = load_yaml(config_path)

    runner = get_section(cfg, "runner")
    input_cfg = get_section(cfg, "input")
    skim = get_section(cfg, "skim")
    cuts = get_section(cfg, "cuts")
    target = get_section(cfg, "target")
    runsel = get_section(cfg, "run_selection")
    output = get_section(cfg, "output")

    runner_script = path_for_command(repo_root, runner.get("script", "scripts/run_all_hipo.sh"))
    runner_check = repo_relative_or_absolute(repo_root, runner_script)

    if not os.path.isfile(runner_check):
        print(f"Error: runner script does not exist: {runner_check}")
        sys.exit(1)

    dataset = as_string(input_cfg.get("dataset", ""))
    input_path_raw = as_string(input_cfg.get("path", ""))
    input_path = repo_relative_or_absolute(repo_root, input_path_raw)
    period = normalize_lower(input_cfg.get("period", "none"))

    mode = normalize_lower(skim.get("mode", ""))
    region = normalize_lower(skim.get("region", "all"))
    pid = as_string(skim.get("pid", ""))

    electrontree = as_string(cuts.get("electrontree", 0))
    detpidcut = as_string(cuts.get("detpidcut", 0))
    kincut = as_string(cuts.get("kincut", 0))

    target_map = path_for_command(repo_root, target.get("map", ""))
    target_map_check = repo_relative_or_absolute(repo_root, target_map)
    polsource = normalize_lower(target.get("polsource", "offline"))
    targetfilter = normalize_lower(target.get("targetfilter", "all"))

    run_mode = normalize_lower(runsel.get("mode", "all"))
    output_format = normalize_lower(output.get("format", "root"))

    yaml_jobtag = as_string(output.get("jobtag", ""))
    jobtag = as_string(args.jobtag) if as_string(args.jobtag) != "" else yaml_jobtag
    if jobtag == "":
        jobtag = "slurm"

    if mode not in ["sidis", "inclusive", "dihadron"]:
        print(f"Error: skim.mode must be sidis, inclusive, or dihadron. Got: {mode}")
        sys.exit(1)

    if period not in ["none", "auto", "summer22", "fall22", "spring23"]:
        print(f"Error: input.period is invalid: {period}")
        sys.exit(1)

    if detpidcut not in ["0", "1"]:
        print(f"Error: cuts.detpidcut must be 0 or 1. Got: {detpidcut}")
        sys.exit(1)

    if kincut not in ["0", "1"]:
        print(f"Error: cuts.kincut must be 0 or 1. Got: {kincut}")
        sys.exit(1)

    if output_format not in ["root", "hipo", "both"]:
        print(f"Error: output.format must be root, hipo, or both. Got: {output_format}")
        sys.exit(1)

    if targetfilter not in ["all", "matched", "0", "1"]:
        print(f"Error: target.targetfilter must be all, matched, 0, or 1. Got: {targetfilter}")
        sys.exit(1)

    if mode == "sidis" and pid == "":
        print("Error: skim.mode sidis requires skim.pid.")
        sys.exit(1)

    if mode != "sidis":
        electrontree = "0"

    print("")
    print("Resolving HIPO input files...")
    files = collect_hipo_files(input_path)
    original_count = len(files)

    skipped_no_run = 0
    skipped_no_target = 0

    if targetfilter in ["matched", "1"]:
        if target_map == "":
            print("Error: targetfilter is matched, but target.map is empty.")
            sys.exit(1)

        if not os.path.isfile(target_map_check):
            print(f"Error: target map does not exist: {target_map_check}")
            sys.exit(1)

        target_runs = read_target_map_runs(target_map_check)

        if not target_runs:
            print(f"Error: no run numbers could be read from target map: {target_map_check}")
            sys.exit(1)

        files, skipped_no_run, skipped_no_target = filter_by_target_map(files, target_runs)

    files = filter_by_run_selection(files, run_mode, runsel, repo_root)

    if args.max_files > 0:
        files = files[: args.max_files]

    if not files:
        print("Error: no HIPO files remain after target/run/max-files filtering.")
        print(f"Original HIPO files:       {original_count}")
        print(f"Skipped, no run tag:       {skipped_no_run}")
        print(f"Skipped, not in targetmap: {skipped_no_target}")
        sys.exit(1)

    safe_jobtag = sanitize_tag(jobtag)
    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")

    filelist_dir = repo_root / "logs" / "slurm_filelists"
    slurm_log_dir = repo_root / "logs" / "slurm"

    filelist_dir.mkdir(parents=True, exist_ok=True)
    slurm_log_dir.mkdir(parents=True, exist_ok=True)

    filelist_path = filelist_dir / f"rgcskim_{safe_jobtag}_{timestamp}.txt"

    with open(filelist_path, "w", encoding="utf-8") as f:
        for path in files:
            f.write(path + "\n")

    sbatch_script = repo_root / "scripts" / "batch" / "rgcskim_array.sbatch"

    if not sbatch_script.is_file():
        print(f"Error: Slurm worker script does not exist: {sbatch_script}")
        sys.exit(1)

    job_name = sanitize_tag(f"rgcskim_{mode}_{region}_{safe_jobtag}")[:80]
    n_files = len(files)

    if n_files == 1:
        array_spec = "0"
    else:
        array_spec = f"0-{n_files - 1}"

    if args.array_limit > 0 and n_files > 1:
        array_spec = f"{array_spec}%{args.array_limit}"

    stdout_pattern = str(slurm_log_dir / f"{job_name}_%A_%a.out")
    stderr_pattern = str(slurm_log_dir / f"{job_name}_%A_%a.err")

    export_items = {
        "RGC_SKIM_REPO": str(repo_root),
        "RGC_SKIM_FILELIST": str(filelist_path),
        "RGC_SKIM_RUNNER": runner_script,
        "RGC_SKIM_DATASET": dataset,
        "RGC_SKIM_MODE": mode,
        "RGC_SKIM_REGION": region,
        "RGC_SKIM_PID": pid,
        "RGC_SKIM_ELECTRONTREE": electrontree,
        "RGC_SKIM_DETPIDCUT": detpidcut,
        "RGC_SKIM_KINCUT": kincut,
        "RGC_SKIM_TARGETMAP": target_map,
        "RGC_SKIM_POLSOURCE": polsource,
        "RGC_SKIM_TARGETFILTER": targetfilter,
        "RGC_SKIM_PERIOD": period,
        "RGC_SKIM_OUTFORMAT": output_format,
        "RGC_SKIM_JOBTAG": jobtag,
    }

    export_arg = "ALL," + ",".join(f"{key}={value}" for key, value in export_items.items())

    sbatch_cmd = [
        "sbatch",
        "--partition",
        args.partition,
        "--time",
        args.time,
        "--cpus-per-task",
        str(args.cpus_per_task),
        "--mem",
        args.mem,
        "--job-name",
        job_name,
        "--output",
        stdout_pattern,
        "--error",
        stderr_pattern,
        "--array",
        array_spec,
        "--export",
        export_arg,
    ]

    if as_string(args.account) != "":
        sbatch_cmd.extend(["--account", args.account])

    sbatch_cmd.append(str(sbatch_script))

    print("")
    print("Resolved Slurm skim submission")
    print("------------------------------------------------------------")
    print(f"Repo:                 {repo_root}")
    print(f"Config:               {config_path}")
    print(f"Runner:               {runner_script}")
    print(f"Input path:           {input_path}")
    print(f"Original HIPO files:  {original_count}")
    print(f"Selected HIPO files:  {n_files}")
    print(f"Skipped, no run tag:  {skipped_no_run}")
    print(f"Skipped, no target:   {skipped_no_target}")
    print(f"Mode:                 {mode}")
    print(f"Region:               {region if mode == 'inclusive' else '(not used)'}")
    print(f"PID:                  {pid if mode == 'sidis' else '(not used)'}")
    print(f"Detector/PID cut:     {detpidcut}")
    print(f"Kinematic cut:        {kincut}")
    print(f"Target filter:        {targetfilter}")
    print(f"Period:               {period}")
    print(f"Output format:        {output_format}")
    print(f"Job tag:              {jobtag}")
    print(f"Slurm account:        {args.account}")
    print(f"Slurm partition:      {args.partition}")
    print(f"Slurm time:           {args.time}")
    print(f"Slurm mem:            {args.mem}")
    print(f"Array spec:           {array_spec}")
    print(f"File list:            {filelist_path}")
    print(f"Slurm stdout:         {stdout_pattern}")
    print(f"Slurm stderr:         {stderr_pattern}")
    print("")
    print("sbatch command")
    print("------------------------------------------------------------")
    print(shlex.join(sbatch_cmd))
    print("")

    if args.dry_run:
        print("Dry run only. No Slurm job was submitted.")
        return

    result = subprocess.run(sbatch_cmd, check=False)

    if result.returncode != 0:
        print("")
        print(f"Error: sbatch failed with exit code {result.returncode}")
        sys.exit(result.returncode)


if __name__ == "__main__":
    main()
