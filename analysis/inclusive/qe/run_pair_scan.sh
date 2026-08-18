#!/usr/bin/env bash
#******************************************************************
#* run_pair_scan.sh
#*
#* Purpose:
#*   YAML-driven driver for qe_pair_scan.cc plus optional QADB merge
#*   and diagnostic plotting.
#*
#* Cut-profile convention:
#*   kin0 = inclusive QE detector/PID cuts only.
#*   kin1 = inclusive QE detector/PID cuts plus tight QE kinematic cuts.
#*
#* Important:
#*   The YAML may list both detector and kinematic cuts permanently.
#*   qe_pair_scan.cc receives --cut-profile and applies the kinematic
#*   block only for cut_profile=kin1.
#******************************************************************

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
EXEC_DEFAULT="$SCRIPT_DIR/qe_pair_scan"
EXEC="$EXEC_DEFAULT"

CONFIG=""
DO_BUILD=0
DRY_RUN=0

# General controls
dataset="spring23"
input_sample="kin0_skim"
cut_profile="kin0"
indir=""
contains=".root"
runs="all"

# Output controls
auto_output=1
output_dir_base="outputs/pair_scan"
output_prefix=""
output_root=""
run_csv=""
pair_csv=""
qadb_run_csv=""
pair_qadb_csv=""
plot_outdir=""
plot_prefix=""

# Optional post steps
run_qadb=1
run_plots=1
qadb_script="tools/qadb/qadb_pair_postprocess.py"
plot_script="tools/qadb/plot_qadb_pair_diagnostics.py"
qadb_cache_dir="qadb_cache"
qadb_dir=""
misc_table=""
qa_table=""

# Pair scan controls
charge_branch="fcupgated_delta"
charge_rel_err="0.002"
require_target_match="1"
require_is_qe="1"
require_fd_status="0"
min_abs_delta_q="0.0"
min_charge="1.0"
min_final="1"
max_pairs="-1"

# Detector/PID cuts, always active for kin0 and kin1
chi2_max="3.0"
detector_pe_min="0.0"
lv_min="14.0"
lw_min="14.0"
pcal_min="0.07"
sf_max="0.28"
use_sf_cut="1"
use_ecin_pcal_cut="1"
ecin_slope="-0.625"
ecin_intercept="0.15"

# Tight kinematic cuts, applied only for cut_profile=kin1
theta_min="7.8"
theta_max="8.2"
vz_min="-5.758"
vz_max="1.5165"
pe_min="2.0"
q2_min="1.9433"
q2_max="2.0574"
w_min="0.0"
w_max="1.073"

# xB counting/plotting range, active for both profiles
xmin="0.85"
xmax="1.15"

# Plot filters
top_n="50"
include_rejected="0"
pair_class="pos_neg"
min_abs_deltaq_plot="0.08"
max_abs_r2q_plot="5.0"
yield_ratio_min="0.90"
yield_ratio_max="1.10"

print_usage(){
    cat <<USAGE
Usage:
  ./run_pair_scan.sh --config <yaml> [--build] [--dry-run]

Key YAML fields:
  dataset: spring23
  input_sample: kin0_skim
  cut_profile: kin0|kin1
  indir: <directory containing per-run skim ROOT files>

The YAML may list all detector and kinematic cuts. Detector/PID cuts are
active for both profiles; tight kinematic cuts are active only for kin1.

Outputs, when auto_output: 1:
  outputs/pair_scan/<cut_profile>/qe_pair_scan_<dataset>_<cut_profile>.root
  outputs/pair_scan/<cut_profile>/qe_pair_scan_<dataset>_<cut_profile>_run_summary.csv
  outputs/pair_scan/<cut_profile>/qe_pair_scan_<dataset>_<cut_profile>_pair_summary.csv
  outputs/pair_scan/<cut_profile>/qe_pair_scan_<dataset>_<cut_profile>_pair_summary_qadb.csv
  outputs/pair_scan/<cut_profile>/diagnostics/*
USAGE
}

yaml_get(){
    local key="$1"
    local file="$2"
    awk -v k="$key" '
        /^[[:space:]]*#/ {next}
        /^[[:space:]]*$/ {next}
        {
            line=$0
            sub(/^[[:space:]]*/, "", line)
            if(line ~ "^" k "[[:space:]]*:"){
                sub("^" k "[[:space:]]*:[[:space:]]*", "", line)
                sub(/[[:space:]]+#.*$/, "", line)
                gsub(/^[[:space:]"'"'"']+/, "", line)
                gsub(/[[:space:]"'"'"']+$/, "", line)
                print line
                exit
            }
        }
    ' "$file"
}

set_from_yaml_if_present(){
    local var="$1"
    local key="$2"
    local val=""
    val=$(yaml_get "$key" "$CONFIG" || true)
    if [ -n "$val" ]; then
        printf -v "$var" '%s' "$val"
    fi
}

load_config(){
    if [ -z "$CONFIG" ]; then return; fi
    if [ ! -f "$CONFIG" ]; then
        echo "Error: config file not found: $CONFIG"
        exit 1
    fi

    set_from_yaml_if_present dataset dataset
    set_from_yaml_if_present input_sample input_sample
    set_from_yaml_if_present cut_profile cut_profile
    set_from_yaml_if_present indir indir
    set_from_yaml_if_present contains contains
    set_from_yaml_if_present runs runs

    set_from_yaml_if_present auto_output auto_output
    set_from_yaml_if_present output_dir_base output_dir_base
    set_from_yaml_if_present output_prefix output_prefix
    set_from_yaml_if_present output_root output_root
    set_from_yaml_if_present run_csv run_csv
    set_from_yaml_if_present pair_csv pair_csv
    set_from_yaml_if_present qadb_run_csv qadb_run_csv
    set_from_yaml_if_present pair_qadb_csv pair_qadb_csv
    set_from_yaml_if_present plot_outdir plot_outdir
    set_from_yaml_if_present plot_prefix plot_prefix

    set_from_yaml_if_present run_qadb run_qadb
    set_from_yaml_if_present run_plots run_plots
    set_from_yaml_if_present qadb_script qadb_script
    set_from_yaml_if_present plot_script plot_script
    set_from_yaml_if_present qadb_cache_dir qadb_cache_dir
    set_from_yaml_if_present qadb_dir qadb_dir
    set_from_yaml_if_present misc_table misc_table
    set_from_yaml_if_present qa_table qa_table

    set_from_yaml_if_present charge_branch charge_branch
    set_from_yaml_if_present charge_rel_err charge_rel_err
    set_from_yaml_if_present require_target_match require_target_match
    set_from_yaml_if_present require_is_qe require_is_qe
    set_from_yaml_if_present require_fd_status require_fd_status
    set_from_yaml_if_present min_abs_delta_q min_abs_delta_q
    set_from_yaml_if_present min_charge min_charge
    set_from_yaml_if_present min_final min_final
    set_from_yaml_if_present max_pairs max_pairs

    set_from_yaml_if_present chi2_max chi2_max
    set_from_yaml_if_present detector_pe_min detector_pe_min
    set_from_yaml_if_present lv_min lv_min
    set_from_yaml_if_present lw_min lw_min
    set_from_yaml_if_present pcal_min pcal_min
    set_from_yaml_if_present sf_max sf_max
    set_from_yaml_if_present use_sf_cut use_sf_cut
    set_from_yaml_if_present use_ecin_pcal_cut use_ecin_pcal_cut
    set_from_yaml_if_present ecin_slope ecin_slope
    set_from_yaml_if_present ecin_intercept ecin_intercept

    set_from_yaml_if_present theta_min theta_min
    set_from_yaml_if_present theta_max theta_max
    set_from_yaml_if_present vz_min vz_min
    set_from_yaml_if_present vz_max vz_max
    set_from_yaml_if_present pe_min pe_min
    set_from_yaml_if_present q2_min q2_min
    set_from_yaml_if_present q2_max q2_max
    set_from_yaml_if_present w_min w_min
    set_from_yaml_if_present w_max w_max

    set_from_yaml_if_present xmin xmin
    set_from_yaml_if_present xmax xmax

    set_from_yaml_if_present top_n top_n
    set_from_yaml_if_present include_rejected include_rejected
    set_from_yaml_if_present pair_class pair_class
    set_from_yaml_if_present min_abs_deltaq_plot min_abs_deltaq_plot
    set_from_yaml_if_present max_abs_r2q_plot max_abs_r2q_plot
    set_from_yaml_if_present yield_ratio_min yield_ratio_min
    set_from_yaml_if_present yield_ratio_max yield_ratio_max
}

bool01(){
    local v="$1"
    case "$v" in
        1|true|TRUE|yes|YES) return 0 ;;
        0|false|FALSE|no|NO) return 1 ;;
        *) echo "Error: expected 0/1 boolean, got: $v"; exit 1 ;;
    esac
}

validate_cut_profile(){
    case "$cut_profile" in
        kin0|kin1) ;;
        *) echo "Error: cut_profile must be kin0 or kin1, got: $cut_profile"; exit 1 ;;
    esac
}

resolve_paths(){
    if [ -z "$indir" ]; then
        echo "Error: indir is required in YAML or command line."
        exit 1
    fi
    if [ ! -d "$indir" ]; then
        echo "Error: indir is not a directory: $indir"
        exit 1
    fi

    if bool01 "$auto_output"; then
        local base="${output_dir_base%/}/${cut_profile}"
        local prefix="${output_prefix:-qe_pair_scan_${dataset}_${cut_profile}}"
        [ -n "$output_root" ]   || output_root="${base}/${prefix}.root"
        [ -n "$run_csv" ]       || run_csv="${base}/${prefix}_run_summary.csv"
        [ -n "$pair_csv" ]      || pair_csv="${base}/${prefix}_pair_summary.csv"
        [ -n "$qadb_run_csv" ]  || qadb_run_csv="${base}/qadb_${dataset}_run_metadata.csv"
        [ -n "$pair_qadb_csv" ] || pair_qadb_csv="${base}/${prefix}_pair_summary_qadb.csv"
        [ -n "$plot_outdir" ]   || plot_outdir="${base}/diagnostics"
        [ -n "$plot_prefix" ]   || plot_prefix="${dataset}_${cut_profile}"
    fi
}

# First pass.
args=("$@")
i=0
while [ $i -lt ${#args[@]} ]; do
    a="${args[$i]}"
    case "$a" in
        --config) i=$((i+1)); CONFIG="${args[$i]}" ;;
        --build) DO_BUILD=1 ;;
        --exec) i=$((i+1)); EXEC="${args[$i]}" ;;
        --dry-run) DRY_RUN=1 ;;
        --help|-h) print_usage; exit 0 ;;
    esac
    i=$((i+1))
done

load_config

# Second pass overrides.
i=0
while [ $i -lt ${#args[@]} ]; do
    a="${args[$i]}"
    case "$a" in
        --config|--exec) i=$((i+2)); continue ;;
        --build|--dry-run) i=$((i+1)); continue ;;
        --dataset) i=$((i+1)); dataset="${args[$i]}" ;;
        --input-sample) i=$((i+1)); input_sample="${args[$i]}" ;;
        --cut-profile) i=$((i+1)); cut_profile="${args[$i]}" ;;
        --indir) i=$((i+1)); indir="${args[$i]}" ;;
        --contains) i=$((i+1)); contains="${args[$i]}" ;;
        --runs) i=$((i+1)); runs="${args[$i]}" ;;
        *) ;;
    esac
    i=$((i+1))
done

validate_cut_profile
resolve_paths

if [ "$DO_BUILD" -eq 1 ]; then
    echo "Building qe_pair_scan with make..."
    make -C "$SCRIPT_DIR" qe_pair_scan || make -C "$SCRIPT_DIR"
fi

if [ ! -x "$EXEC" ]; then
    echo "Error: executable not found or not executable: $EXEC"
    echo "Run: make qe_pair_scan"
    exit 1
fi

mkdir -p "$(dirname "$output_root")"
mkdir -p "$(dirname "$run_csv")"
mkdir -p "$(dirname "$pair_csv")"

CMD=("$EXEC")
CMD+=(--indir "$indir")
CMD+=(--contains "$contains")
CMD+=(--runs "$runs")
CMD+=(--output-root "$output_root")
CMD+=(--run-csv "$run_csv")
CMD+=(--pair-csv "$pair_csv")
CMD+=(--cut-profile "$cut_profile")
CMD+=(--charge-branch "$charge_branch")
CMD+=(--charge-rel-err "$charge_rel_err")
CMD+=(--require-target-match "$require_target_match")
CMD+=(--require-is-qe "$require_is_qe")
CMD+=(--require-fd-status "$require_fd_status")
CMD+=(--min-abs-delta-q "$min_abs_delta_q")
CMD+=(--min-charge "$min_charge")
CMD+=(--min-final "$min_final")
CMD+=(--max-pairs "$max_pairs")

# Detector/PID cuts, always active in C++.
CMD+=(--chi2-max "$chi2_max")
CMD+=(--detector-pe-min "$detector_pe_min")
CMD+=(--lv-min "$lv_min")
CMD+=(--lw-min "$lw_min")
CMD+=(--pcal-min "$pcal_min")
CMD+=(--sf-max "$sf_max")
CMD+=(--use-sf-cut "$use_sf_cut")
CMD+=(--use-ecin-pcal-cut "$use_ecin_pcal_cut")
CMD+=(--ecin-slope "$ecin_slope")
CMD+=(--ecin-intercept "$ecin_intercept")

# Tight kinematic cuts. C++ applies them only when cut_profile=kin1.
CMD+=(--theta-min "$theta_min")
CMD+=(--theta-max "$theta_max")
CMD+=(--vz-min "$vz_min")
CMD+=(--vz-max "$vz_max")
CMD+=(--pe-min "$pe_min")
CMD+=(--q2-min "$q2_min")
CMD+=(--q2-max "$q2_max")
CMD+=(--w-min "$w_min")
CMD+=(--w-max "$w_max")

# Counting/plotting xB interval.
CMD+=(--xmin "$xmin")
CMD+=(--xmax "$xmax")

echo "================================================================================"
echo "QE many-run pair scanner"
echo "================================================================================"
echo "Working dir:       $(pwd)"
echo "Script dir:        $SCRIPT_DIR"
echo "Config:            ${CONFIG:-none}"
echo "Executable:        $EXEC"
echo "Dataset:           $dataset"
echo "Input sample:      $input_sample"
echo "Cut profile:       $cut_profile"
echo "Input dir:         $indir"
echo "Output ROOT:       $output_root"
echo "Run CSV:           $run_csv"
echo "Pair CSV:          $pair_csv"
echo "QADB pair CSV:     $pair_qadb_csv"
echo "Plot outdir:       $plot_outdir"
echo "================================================================================"
printf 'Command:'
printf ' %q' "${CMD[@]}"
echo

echo "================================================================================"

if [ "$DRY_RUN" -eq 1 ]; then
    echo "Dry run only. Command not executed."
    exit 0
fi

"${CMD[@]}"

if bool01 "$run_qadb"; then
    if [ ! -f "$qadb_script" ]; then
        echo "Error: QADB script not found: $qadb_script"
        exit 1
    fi
    QADB_CMD=(python3 "$qadb_script" --pair-csv "$pair_csv" --out-run-csv "$qadb_run_csv" --out-pair-csv "$pair_qadb_csv" --cache-dir "$qadb_cache_dir")
    if [ -n "$qadb_dir" ]; then QADB_CMD+=(--qadb-dir "$qadb_dir"); fi
    if [ -n "$misc_table" ]; then QADB_CMD+=(--misc-table "$misc_table"); fi
    if [ -n "$qa_table" ]; then QADB_CMD+=(--qa-table "$qa_table"); fi

    echo "================================================================================"
    echo "Running QADB postprocess"
    printf 'Command:'
    printf ' %q' "${QADB_CMD[@]}"
    echo
    echo "================================================================================"
    "${QADB_CMD[@]}"
fi

if bool01 "$run_plots"; then
    if [ ! -f "$plot_script" ]; then
        echo "Error: plot script not found: $plot_script"
        exit 1
    fi
    PLOT_INPUT="$pair_csv"
    if bool01 "$run_qadb"; then
        PLOT_INPUT="$pair_qadb_csv"
    fi
    PLOT_CMD=(python3 "$plot_script" --pair-csv "$PLOT_INPUT" --outdir "$plot_outdir" --prefix "$plot_prefix" --top-n "$top_n" --include-rejected "$include_rejected" --pair-class "$pair_class" --min-abs-deltaq "$min_abs_deltaq_plot" --max-abs-r2q-plot "$max_abs_r2q_plot" --yield-ratio-min "$yield_ratio_min" --yield-ratio-max "$yield_ratio_max")

    echo "================================================================================"
    echo "Running pair diagnostic plots"
    printf 'Command:'
    printf ' %q' "${PLOT_CMD[@]}"
    echo
    echo "================================================================================"
    "${PLOT_CMD[@]}"
fi