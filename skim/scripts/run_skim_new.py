#!/usr/bin/env python3

"""
run_skim.py
YAML wrapper for run_all_hipo_new.sh.
Users edit skim.yaml and run:
    ./run_skim.py skim.yaml --dry-run
    ./run_skim.py skim.yaml
This script does not skim events directly. It reads the YAML configuration,
validates it, builds the equivalent run_all_hipo_new.sh command, prints it, and
then optionally executes it.
"""

import argparse
import os
import re
import shlex
import subprocess
import sys


def load_yaml(path):
    try:
        import yaml
    except ImportError:
        print("Error: Python module 'yaml' was not found.")
        print("")
        print("This wrapper needs PyYAML.")
        print("")
        print("Try:")
        print("  python3 -m pip install --user PyYAML")
        sys.exit(1)

    if not os.path.isfile(path):
        print(f"Error: YAML file does not exist: {path}")
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


def as_int_string(value, field_name):
    value_str = as_string(value)

    if value_str == "":
        return ""

    try:
        int(value_str)
    except ValueError:
        print(f"Error: '{field_name}' must be an integer-like value.")
        print(f"Got: {value}")
        sys.exit(1)

    return value_str


def normalize_lower(value):
    return as_string(value).lower()


def require_file_or_dir(path, label):
    if path == "":
        print(f"Error: missing required path for {label}.")
        sys.exit(1)

    if not os.path.exists(path):
        print(f"Error: {label} does not exist:")
        print(f"  {path}")
        sys.exit(1)


def detect_period_from_path(path):
    lower_path = path.lower()

    if "summer22" in lower_path:
        return "summer22"

    if "fall22" in lower_path:
        return "fall22"

    if "spring23" in lower_path:
        return "spring23"

    return ""


def sanitize_folder_tag(raw):
    return re.sub(r"[^A-Za-z0-9._-]", "_", raw)


def infer_dataset_tag(dataset, hipo_path):
    if dataset != "":
        return dataset

    if os.path.isdir(hipo_path):
        return os.path.basename(os.path.normpath(hipo_path))

    return os.path.basename(os.path.dirname(os.path.normpath(hipo_path)))


def build_predicted_output(dataset, hipo_path, period, mode, region, pid, electron_tree, det_pid_cut, kin_cut, targetfilter, jobtag, output_format, output_path):
    dataset_tag = infer_dataset_tag(dataset, hipo_path)

    if mode == "sidis":
        if pid == "211":
            channel_tag = "epip"
        elif pid == "-211":
            channel_tag = "epim"
        elif pid == "321":
            channel_tag = "ekp"
        elif pid == "-321":
            channel_tag = "ekm"
        elif pid == "2212":
            channel_tag = "ep"
        else:
            channel_tag = f"epid{pid}"
    elif mode == "dihadron":
        channel_tag = "epippim"
    else:
        region_for_tag = region
        if region_for_tag == "resonance":
            region_for_tag = "res"
        if region_for_tag in ["quasi", "quasielastic", "quasi-elastic"]:
            region_for_tag = "qe"
        channel_tag = f"inclusive_{region_for_tag}"

    if targetfilter in ["all", "0"]:
        target_filter_tag = "tf0"
    else:
        target_filter_tag = "tf1"

    if mode in ["inclusive", "dihadron"]:
        electron_tree_for_tag = "0"
    else:
        electron_tree_for_tag = electron_tree

    option_tag = f"eT{electron_tree_for_tag}pid{det_pid_cut}kin{kin_cut}{target_filter_tag}"
    subdir_tag = f"skim_{dataset_tag}_{channel_tag}_{option_tag}"

    if period == "none":
        period_tag = ""
    elif period == "auto":
        period_tag = detect_period_from_path(hipo_path)
    else:
        period_tag = period

    root_out_dir = os.path.join("rootfiles", subdir_tag)
    hipo_out_dir = os.path.join("hipofiles", subdir_tag)
    log_out_dir = os.path.join("logs", subdir_tag)

    if period_tag != "":
        root_out_dir = os.path.join(root_out_dir, period_tag)
        hipo_out_dir = os.path.join(hipo_out_dir, period_tag)
        log_out_dir = os.path.join(log_out_dir, period_tag)

    if jobtag != "":
        jobtag_clean = sanitize_folder_tag(jobtag)
        root_out_dir = os.path.join(root_out_dir, jobtag_clean)
        hipo_out_dir = os.path.join(hipo_out_dir, jobtag_clean)
        log_out_dir = os.path.join(log_out_dir, jobtag_clean)
        
    if output_path != "":
        root_out_dir = os.path.join(output_path, root_out_dir)
        hipo_out_dir = os.path.join(output_path, hipo_out_dir, period_tag)
        log_out_dir = os.path.join(output_path, log_out_dir)

    return {
        "period_tag": period_tag if period_tag else "none",
        "root_out_dir": root_out_dir,
        "hipo_out_dir": hipo_out_dir,
        "log_out_dir": log_out_dir,
        "root_file_pattern": f"{dataset_tag}_{channel_tag}_{option_tag}_<run>.root",
        "hipo_file_pattern": f"{dataset_tag}_{channel_tag}_{option_tag}_<run>.hipo",
        "log_file_pattern": f"{dataset_tag}_{channel_tag}_{option_tag}_<run>.log",
        "output_format": output_format,
    }


def validate_choice(value, allowed, field_name):
    if value not in allowed:
        print(f"Error: invalid {field_name}: {value}")
        print("Allowed values:")
        for item in allowed:
            print(f"  {item}")
        sys.exit(1)


def build_command(cfg):
    runner = get_section(cfg, "runner")
    input_cfg = get_section(cfg, "input")
    skim = get_section(cfg, "skim")
    cuts = get_section(cfg, "cuts")
    target = get_section(cfg, "target")
    runsel = get_section(cfg, "run_selection")
    output = get_section(cfg, "output")

    runner_script = as_string(runner.get("script", "./run_all_hipo_new.sh"))
    dataset = as_string(input_cfg.get("dataset", ""))
    hipo_path = as_string(input_cfg.get("path", ""))
    period = normalize_lower(input_cfg.get("period", "none"))

    mode = normalize_lower(skim.get("mode", ""))
    region = normalize_lower(skim.get("region", ""))
    pid = as_int_string(skim.get("pid", ""), "skim.pid")

    electron_tree = as_int_string(cuts.get("electrontree", 0), "cuts.electrontree")
    det_pid_cut = as_int_string(cuts.get("detpidcut", 0), "cuts.detpidcut")
    kin_cut = as_int_string(cuts.get("kincut", 0), "cuts.kincut")

    target_map = as_string(target.get("map", ""))
    polsource = normalize_lower(target.get("polsource", "offline"))
    targetfilter = normalize_lower(target.get("targetfilter", "all"))

    run_mode = normalize_lower(runsel.get("mode", "all"))
    run_single = as_int_string(runsel.get("run", ""), "run_selection.run")
    run_start = as_int_string(runsel.get("start", ""), "run_selection.start")
    run_end = as_int_string(runsel.get("end", ""), "run_selection.end")
    runlist = as_string(runsel.get("runlist", ""))
    runs_value = runsel.get("runs", [])

    output_format = normalize_lower(output.get("format", "root"))
    jobtag = as_string(output.get("jobtag", ""))
    output_path = as_string(output.get("out_path", ""))

    # ------------------------------------------------------------
    # Validate basic inputs
    # ------------------------------------------------------------

    require_file_or_dir(runner_script, "runner.script")
    require_file_or_dir(hipo_path, "input.path")

    validate_choice(mode, ["sidis", "inclusive", "dihadron"], "skim.mode")
    validate_choice(period, ["none", "auto", "summer22", "fall22", "spring23"], "input.period")
    validate_choice(polsource, ["offline", "online"], "target.polsource")
    validate_choice(targetfilter, ["all", "matched", "0", "1"], "target.targetfilter")
    validate_choice(run_mode, ["all", "test_one_run", "single", "range", "list", "file"], "run_selection.mode")
    validate_choice(output_format, ["root", "hipo", "both"], "output.format")

    if electron_tree not in ["0", "1"]:
        print("Error: cuts.electrontree must be 0 or 1.")
        sys.exit(1)

    if det_pid_cut not in ["0", "1"]:
        print("Error: cuts.detpidcut must be 0 or 1.")
        sys.exit(1)

    if kin_cut not in ["0", "1"]:
        print("Error: cuts.kincut must be 0 or 1.")
        sys.exit(1)

    # ------------------------------------------------------------
    # Validate mode specific settings
    # ------------------------------------------------------------

    if mode == "sidis":
        if pid == "":
            print("Error: skim.mode sidis requires skim.pid.")
            print("Example PID values: 211, -211, 321, -321, 2212")
            sys.exit(1)

        if region not in ["", "all"]:
            print("Error: skim.region is only used for inclusive mode.")
            print("For SIDIS, leave region empty.")
            sys.exit(1)

    elif mode == "inclusive":
        if pid != "":
            print("Error: skim.pid should be empty for inclusive mode.")
            sys.exit(1)

        if region == "":
            region = "all"

        valid_regions = [
            "all",
            "dis",
            "res",
            "resonance",
            "qe",
            "quasi",
            "quasielastic",
            "quasi-elastic",
        ]

        validate_choice(region, valid_regions, "skim.region")

        if electron_tree == "1":
            print("Error: cuts.electrontree = 1 is only meaningful for SIDIS.")
            print("Inclusive mode already writes an electron-only tree.")
            sys.exit(1)

    elif mode == "dihadron":
        if pid != "":
            print("Error: skim.pid should be empty for dihadron mode.")
            sys.exit(1)

        if region not in ["", "all"]:
            print("Error: skim.region is only used for inclusive mode.")
            print("For dihadron, leave region empty.")
            sys.exit(1)

    # ------------------------------------------------------------
    # Validate target map settings
    # ------------------------------------------------------------

    if targetfilter in ["matched", "1"] and target_map == "":
        print("Error: target.targetfilter matched requires target.map.")
        sys.exit(1)

    if target_map != "" and not os.path.isfile(target_map):
        print("Error: target.map does not exist:")
        print(f"  {target_map}")
        sys.exit(1)

    # ------------------------------------------------------------
    # Validate run selection
    # ------------------------------------------------------------

    if run_mode == "all":
        pass

    elif run_mode == "test_one_run":
        pass

    elif run_mode == "single":
        if run_single == "":
            print("Error: run_selection.mode single requires run_selection.run.")
            sys.exit(1)

    elif run_mode == "range":
        if run_start == "" or run_end == "":
            print("Error: run_selection.mode range requires run_selection.start and run_selection.end.")
            sys.exit(1)

        if int(run_start) > int(run_end):
            print("Error: run_selection.start is greater than run_selection.end.")
            print(f"  start = {run_start}")
            print(f"  end   = {run_end}")
            sys.exit(1)

    elif run_mode == "list":
        if runs_value is None:
            runs_value = []

        if not isinstance(runs_value, list):
            print("Error: run_selection.runs must be a YAML list.")
            print("Example:")
            print("  runs: [16270, 16271, 16275]")
            sys.exit(1)

        if len(runs_value) == 0:
            print("Error: run_selection.mode list requires run_selection.runs.")
            sys.exit(1)

        normalized_runs = []
        for i, item in enumerate(runs_value):
            item_str = as_int_string(item, f"run_selection.runs[{i}]")
            if item_str == "":
                print(f"Error: empty value in run_selection.runs[{i}].")
                sys.exit(1)
            normalized_runs.append(item_str)

        runs_value = normalized_runs

    elif run_mode == "file":
        if runlist == "":
            print("Error: run_selection.mode file requires run_selection.runlist.")
            sys.exit(1)

        if not os.path.isfile(runlist):
            print("Error: run_selection.runlist does not exist:")
            print(f"  {runlist}")
            sys.exit(1)

    # ------------------------------------------------------------
    # Build command
    # ------------------------------------------------------------

    cmd = [runner_script]

    if dataset != "":
        cmd += ["--dataset", dataset]

    cmd += ["--mode", mode]

    if mode == "sidis":
        cmd += ["--pid", pid]

    if mode == "inclusive":
        cmd += ["--region", region]

    cmd += [hipo_path]

    if mode == "sidis":
        cmd += ["--electrontree", electron_tree]

    cmd += ["--detpidcut", det_pid_cut]
    cmd += ["--kincut", kin_cut]

    if target_map != "":
        cmd += ["--targetmap", target_map]
        cmd += ["--polsource", polsource]

    cmd += ["--targetfilter", targetfilter]
    cmd += ["--period", period]
    cmd += ["--outformat", output_format]

    if run_mode == "test_one_run":
        cmd += ["--test-one-run"]

    elif run_mode == "single":
        cmd += ["--run", run_single]

    elif run_mode == "range":
        cmd += ["--runrange", f"{run_start}:{run_end}"]

    elif run_mode == "list":
        cmd += ["--runs", ",".join(runs_value)]

    elif run_mode == "file":
        cmd += ["--runlist", runlist]

    if jobtag != "":
        cmd += ["--jobtag", jobtag]
        
    if output_path != "":
        cmd += ["--out-path", output_path]

    summary = {
        "runner": runner_script,
        "dataset": dataset if dataset else "(auto)",
        "input": hipo_path,
        "period": period,
        "mode": mode,
        "region": region if mode == "inclusive" else "(not used)",
        "pid": pid if mode == "sidis" else "(not used)",
        "electrontree": electron_tree if mode == "sidis" else "(not used)",
        "detpidcut": det_pid_cut,
        "kincut": kin_cut,
        "output_format": output_format,
        "target_map": target_map if target_map else "(none)",
        "polsource": polsource if target_map else "(not used)",
        "targetfilter": targetfilter,
        "run_selection": run_mode,
        "jobtag": jobtag if jobtag else "(none)",
        "out-path": output_path if output_path else "(none)",
    }

    if run_mode == "test_one_run":
        summary["test_one_run"] = "enabled"
    elif run_mode == "single":
        summary["selected_run"] = run_single
    elif run_mode == "range":
        summary["run_range"] = f"{run_start}:{run_end}"
    elif run_mode == "list":
        summary["selected_runs"] = ",".join(runs_value)
    elif run_mode == "file":
        summary["runlist"] = runlist

    summary["predicted_output"] = build_predicted_output(
        dataset,
        hipo_path,
        period,
        mode,
        region,
        pid,
        electron_tree,
        det_pid_cut,
        kin_cut,
        targetfilter,
        jobtag,
        output_format,
        output_path,
    )

    return cmd, summary


def print_summary(summary, cmd):
    print("")
    print("Resolved skim configuration")
    print("----------------------------------------")
    print(f"Runner:          {summary['runner']}")
    print(f"Dataset:         {summary['dataset']}")
    print(f"Input:           {summary['input']}")
    print(f"Period:          {summary['period']}")
    print(f"Mode:            {summary['mode']}")
    print(f"Region:          {summary['region']}")
    print(f"PID:             {summary['pid']}")
    print(f"Electron tree:   {summary['electrontree']}")
    print(f"Detector/PID cut:{summary['detpidcut']}")
    print(f"Kinematic cut:   {summary['kincut']}")
    print(f"Output format:   {summary['output_format']}")
    print(f"Target map:      {summary['target_map']}")
    print(f"Polarization src:{summary['polsource']}")
    print(f"Target filter:   {summary['targetfilter']}")
    print(f"Run selection:   {summary['run_selection']}")

    if "test_one_run" in summary:
        print("Test one run:    first file passing run/target filters")

    if "selected_run" in summary:
        print(f"Selected run:    {summary['selected_run']}")

    if "run_range" in summary:
        print(f"Run range:       {summary['run_range']}")

    if "selected_runs" in summary:
        print(f"Selected runs:   {summary['selected_runs']}")

    if "runlist" in summary:
        print(f"Run-list file:   {summary['runlist']}")

    print(f"Job tag:         {summary['jobtag']}")
    print(f"Output path:     {summary['out-path']}")
    print("")
    predicted = summary["predicted_output"]
    print("Predicted output")
    print("----------------------------------------")
    print(f"Period tag:       {predicted['period_tag']}")
    if predicted["output_format"] in ["root", "both"]:
        print(f"ROOT output dir:  {predicted['root_out_dir']}")
    if predicted["output_format"] in ["hipo", "both"]:
        print(f"HIPO output dir:  {predicted['hipo_out_dir']}")
    print(f"Log output dir:   {predicted['log_out_dir']}")
    if predicted["output_format"] in ["root", "both"]:
        print(f"ROOT file pattern:{predicted['root_file_pattern']}")
    if predicted["output_format"] in ["hipo", "both"]:
        print(f"HIPO file pattern:{predicted['hipo_file_pattern']}")
    print(f"Log file pattern: {predicted['log_file_pattern']}")
    print("")
    print("Resolved command")
    print("----------------------------------------")
    print(shlex.join(cmd))
    print("")


def main():
    parser = argparse.ArgumentParser(
        description="Run RGC skim using a YAML configuration file."
    )

    parser.add_argument(
        "yaml_file",
        help="Path to skim YAML configuration file, for example skim.yaml",
    )

    parser.add_argument(
        "--dry-run",
        "--dry_run",
        action="store_true",
        help="Print the resolved command but do not run it.",
    )

    parser.add_argument(
        "--yes",
        "-y",
        action="store_true",
        help="Run without asking for confirmation.",
    )

    args = parser.parse_args()

    cfg = load_yaml(args.yaml_file)
    cmd, summary = build_command(cfg)

    print_summary(summary, cmd)

    if args.dry_run:
        print("Dry run only. No skim was executed.")
        return

    if not args.yes:
        answer = input("Proceed? [y/N] ").strip().lower()
        if answer not in ["y", "yes"]:
            print("Cancelled.")
            return

    print("")
    print("Running skim...")
    print("")

    status = subprocess.run(cmd)

    if status.returncode != 0:
        print("")
        print(f"Error: skim command failed with exit code {status.returncode}.")
        sys.exit(status.returncode)

    print("")
    print("Skim command finished successfully.")


if __name__ == "__main__":
    main()
