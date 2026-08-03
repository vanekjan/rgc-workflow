#!/usr/bin/env python3

"""
submit_generic_hipo_array.py

Submit a generic command over HIPO files using a Slurm array job.

Design:
  one Slurm array task = one HIPO file

This tool is independent of the rgcskim-specific skim workflow.
"""

import argparse
import datetime as dt
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys


def as_string(value):
    if value is None:
        return ""
    return str(value).strip()


def sanitize_tag(raw):
    raw = as_string(raw)
    if raw == "":
        return "generic_hipo"

    cleaned = re.sub(r"[^A-Za-z0-9._-]", "_", raw)
    cleaned = cleaned.strip("._-")

    if cleaned == "":
        return "generic_hipo"

    return cleaned


def resolve_path(raw, base=None):
    raw = as_string(raw)

    if raw == "":
        return ""

    path = Path(os.path.expanduser(raw))

    if not path.is_absolute() and base is not None:
        path = Path(base) / path

    return str(path.resolve())


def extract_run_number(path):
    matches = re.findall(r"[0-9]{6}", str(path))

    if not matches:
        return None

    return int(matches[-1])


def parse_run_tokens(values):
    runs = []

    if values is None:
        return runs

    if isinstance(values, str):
        raw_items = re.split(r"[,\s]+", values.strip())
    else:
        raw_items = []

        for value in values:
            raw_items.extend(re.split(r"[,\s]+", str(value).strip()))

    for item in raw_items:
        item = item.strip()

        if item == "":
            continue

        if not re.fullmatch(r"[0-9]+", item):
            print(f"Error: invalid run number: {item}")
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


def parse_runrange(raw):
    raw = as_string(raw)

    if raw == "":
        return None

    parts = re.split(r"[:-]", raw)

    if len(parts) != 2:
        print(f"Error: invalid --runrange value: {raw}")
        print("Use START:END or START-END, for example 17491:17520")
        sys.exit(1)

    start = parts[0].strip()
    end = parts[1].strip()

    if not re.fullmatch(r"[0-9]+", start) or not re.fullmatch(r"[0-9]+", end):
        print(f"Error: invalid --runrange value: {raw}")
        sys.exit(1)

    start_i = int(start)
    end_i = int(end)

    if end_i < start_i:
        print("Error: --runrange end is smaller than start.")
        sys.exit(1)

    return (start_i, end_i)


def read_filelist(path):
    files = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()

            if line == "":
                continue

            files.append(resolve_path(line))

    return files


def collect_hipo_files(input_path, input_type, recursive):
    path = Path(input_path)

    if input_type == "auto":
        if path.is_dir():
            input_type = "dir"
        elif path.is_file() and path.suffix == ".hipo":
            input_type = "hipo"
        elif path.is_file():
            input_type = "filelist"
        else:
            print(f"Error: input does not exist: {path}")
            sys.exit(1)

    if input_type == "hipo":
        if not path.is_file():
            print(f"Error: HIPO input file does not exist: {path}")
            sys.exit(1)

        if path.suffix != ".hipo":
            print(f"Error: input file is not a .hipo file: {path}")
            sys.exit(1)

        return [str(path.resolve())]

    if input_type == "dir":
        if not path.is_dir():
            print(f"Error: input directory does not exist: {path}")
            sys.exit(1)

        if recursive:
            files = sorted(str(p.resolve()) for p in path.rglob("*.hipo"))
        else:
            files = sorted(str(p.resolve()) for p in path.glob("*.hipo"))

        if not files:
            print(f"Error: no .hipo files found in input directory: {path}")
            sys.exit(1)

        return files

    if input_type == "filelist":
        if not path.is_file():
            print(f"Error: input file list does not exist: {path}")
            sys.exit(1)

        files = read_filelist(path)

        if not files:
            print(f"Error: input file list is empty: {path}")
            sys.exit(1)

        return files

    print(f"Error: unsupported input type: {input_type}")
    sys.exit(1)


def filter_files(files, args):
    selected_runs = []

    if as_string(args.run) != "":
        selected_runs.append({int(args.run)})

    runs_from_runs = set(parse_run_tokens(args.runs))

    if runs_from_runs:
        selected_runs.append(runs_from_runs)

    if as_string(args.runlist) != "":
        runlist_path = resolve_path(args.runlist, base=args.work_dir)

        if not os.path.isfile(runlist_path):
            print(f"Error: runlist does not exist: {runlist_path}")
            sys.exit(1)

        selected_runs.append(set(read_runlist(runlist_path)))

    runrange = parse_runrange(args.runrange)

    active_filter = bool(selected_runs or runrange)

    if not active_filter:
        return files, 0

    kept = []
    skipped_no_run = 0

    for path in files:
        run = extract_run_number(path)

        if run is None:
            skipped_no_run += 1
            continue

        keep = True

        for run_set in selected_runs:
            if run not in run_set:
                keep = False
                break

        if keep and runrange is not None:
            start, end = runrange
            if not (start <= run <= end):
                keep = False

        if keep:
            kept.append(path)

    return kept, skipped_no_run


def write_filelist(path, files):
    with open(path, "w", encoding="utf-8") as f:
        for item in files:
            f.write(item + "\n")


def write_shell_env(path, values):
    with open(path, "w", encoding="utf-8") as f:
        f.write("# Auto-generated by submit_generic_hipo_array.py\n")
        f.write("# Do not edit by hand.\n\n")

        for key, value in values.items():
            f.write(f"export {key}={shlex.quote(str(value))}\n")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Submit a generic command over HIPO files as a Slurm array job."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input HIPO file, HIPO directory, or text file list.",
    )

    parser.add_argument(
        "--input-type",
        default="auto",
        choices=["auto", "hipo", "dir", "filelist"],
        help="How to interpret --input. Default: auto",
    )

    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recursively search for .hipo files when --input is a directory.",
    )

    parser.add_argument(
        "--command",
        required=True,
        help="Command template to run for each HIPO file.",
    )

    parser.add_argument(
        "--work-dir",
        default=os.getcwd(),
        help="Directory where the command should run. Default: current directory.",
    )

    parser.add_argument(
        "--output-dir",
        default="",
        help="Directory available as {output_dir}. Default: generic_slurm_output/<jobtag>_<timestamp>/",
    )

    parser.add_argument(
        "--state-dir",
        default="logs/generic_slurm",
        help="Directory for generated filelists, env files, and Slurm logs. Default: logs/generic_slurm",
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
        default="generic_hipo",
        help="Tag used in job names, logs, and default output folders.",
    )

    parser.add_argument(
        "--run",
        default="",
        help="Process only one run number.",
    )

    parser.add_argument(
        "--runs",
        nargs="*",
        default=[],
        help="Process selected runs. Accepts space- or comma-separated values.",
    )

    parser.add_argument(
        "--runrange",
        default="",
        help="Process inclusive run range, for example 17491:17520.",
    )

    parser.add_argument(
        "--runlist",
        default="",
        help="Text file containing run numbers.",
    )

    parser.add_argument(
        "--hipo-path",
        default="/group/clas12/packages/hipo/dev",
        help="HIPO path exported inside the Slurm worker. Default: /group/clas12/packages/hipo/dev",
    )

    parser.add_argument(
        "--no-clas12-env",
        action="store_true",
        help="Do not load the CLAS12 module environment in the Slurm worker.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Create helper files and print sbatch command, but do not submit.",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    work_dir = Path(resolve_path(args.work_dir))

    if not work_dir.is_dir():
        print(f"Error: --work-dir does not exist or is not a directory: {work_dir}")
        sys.exit(1)

    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_jobtag = sanitize_tag(args.jobtag)

    input_path = resolve_path(args.input, base=work_dir)

    files = collect_hipo_files(input_path, args.input_type, args.recursive)
    original_count = len(files)

    files, skipped_no_run = filter_files(files, args)

    if args.max_files > 0:
        files = files[: args.max_files]

    if not files:
        print("Error: no HIPO files remain after filtering.")
        print(f"Original HIPO files: {original_count}")
        print(f"Skipped, no run tag: {skipped_no_run}")
        sys.exit(1)

    state_dir = Path(resolve_path(args.state_dir, base=work_dir))
    filelist_dir = state_dir / "filelists"
    env_dir = state_dir / "env"
    slurm_log_dir = state_dir / "slurm"

    filelist_dir.mkdir(parents=True, exist_ok=True)
    env_dir.mkdir(parents=True, exist_ok=True)
    slurm_log_dir.mkdir(parents=True, exist_ok=True)

    if as_string(args.output_dir) == "":
        output_dir = work_dir / "generic_slurm_output" / f"{safe_jobtag}_{timestamp}"
    else:
        output_dir = Path(resolve_path(args.output_dir, base=work_dir))

    output_dir.mkdir(parents=True, exist_ok=True)

    filelist_path = filelist_dir / f"{safe_jobtag}_{timestamp}.txt"
    env_path = env_dir / f"{safe_jobtag}_{timestamp}.env"

    write_filelist(filelist_path, files)

    script_path = Path(__file__).resolve()
    worker_path = script_path.with_name("generic_hipo_array_worker.sbatch")

    if not worker_path.is_file():
        print(f"Error: worker script not found: {worker_path}")
        sys.exit(1)

    env_values = {
        "GENERIC_FILELIST": str(filelist_path),
        "GENERIC_WORK_DIR": str(work_dir),
        "GENERIC_OUTPUT_DIR": str(output_dir),
        "GENERIC_COMMAND_TEMPLATE": args.command,
        "GENERIC_LOAD_CLAS12_ENV": "0" if args.no_clas12_env else "1",
        "GENERIC_HIPO_PATH": args.hipo_path,
    }

    write_shell_env(env_path, env_values)

    n_files = len(files)

    if n_files == 1:
        array_spec = "0"
    else:
        array_spec = f"0-{n_files - 1}"

    if args.array_limit > 0 and n_files > 1:
        array_spec = f"{array_spec}%{args.array_limit}"

    job_name = sanitize_tag(f"generic_hipo_{safe_jobtag}")[:80]

    stdout_pattern = str(slurm_log_dir / f"{job_name}_%A_%a.out")
    stderr_pattern = str(slurm_log_dir / f"{job_name}_%A_%a.err")

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
        f"ALL,GENERIC_JOB_ENV={env_path}",
    ]

    if as_string(args.account) != "":
        sbatch_cmd.extend(["--account", args.account])

    sbatch_cmd.append(str(worker_path))

    print("")
    print("Resolved generic HIPO Slurm submission")
    print("------------------------------------------------------------")
    print(f"Work dir:             {work_dir}")
    print(f"Input:                {input_path}")
    print(f"Original HIPO files:  {original_count}")
    print(f"Selected HIPO files:  {n_files}")
    print(f"Skipped, no run tag:  {skipped_no_run}")
    print(f"Job tag:              {safe_jobtag}")
    print(f"Slurm account:        {args.account}")
    print(f"Slurm partition:      {args.partition}")
    print(f"Slurm time:           {args.time}")
    print(f"Slurm mem:            {args.mem}")
    print(f"Array spec:           {array_spec}")
    print(f"Load CLAS12 env:      {'no' if args.no_clas12_env else 'yes'}")
    print(f"HIPO path:            {args.hipo_path}")
    print(f"Output dir:           {output_dir}")
    print(f"State dir:            {state_dir}")
    print(f"File list:            {filelist_path}")
    print(f"Env file:             {env_path}")
    print(f"Slurm stdout:         {stdout_pattern}")
    print(f"Slurm stderr:         {stderr_pattern}")
    print("")
    print("Command template")
    print("------------------------------------------------------------")
    print(args.command)
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
