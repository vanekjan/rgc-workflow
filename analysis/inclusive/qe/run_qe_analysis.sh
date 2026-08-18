#!/usr/bin/env bash
#******************************************************************
#* run_qe_analysis.sh
#*
#* Purpose:
#*   Driver script for qe_analysis.cc.
#*
#* Intended location:
#*   ../rgc/analysis/qe/
#*
#* Design:
#*   - Builds/runs the standalone QE tensor two-run analysis.
#*   - Reads options either from command line or from a simple flat YAML file.
#*   - Keeps all run-pair choices and cuts outside the C++ source.
#*
#* Notes about the YAML parser:
#*   This script intentionally supports only simple flat key: value YAML.
#*   That keeps it portable on ifarm without requiring PyYAML.
#*
#* Typical use:
#*   make
#*   ./run_qe_analysis.sh --config configs/qe_analysis.yaml
#*
#* Build and run in one command:
#*   ./run_qe_analysis.sh --build --config configs/qe_analysis.yaml
#******************************************************************

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
EXEC_DEFAULT="$SCRIPT_DIR/qe_analysis"
EXEC="$EXEC_DEFAULT"

CONFIG=""
DO_BUILD=0
DRY_RUN=0

#============================================================
# Analysis configuration variables.
# Values can come from YAML and/or command-line options.
# Command-line options override YAML values.
#============================================================

file1=""
file2=""
run1=""
run2=""
output="outputs/qe_R.root"
output_is_user_set=0

# Auto path/output controls. These are driver-script options only; they are
# not forwarded to qe_analysis.cc.
auto_inputs=1
auto_output=1
pair_subdir=1
skim_dir=""
skim_prefix=""
output_dir_base="outputs"
output_stem="qe_Azz_approx"
output_tag=""
cut_profile="kin1"
profile_subdir=1

q1=""
q2=""
qerr1=""
qerr2=""

p1=""
p2=""
perr1=""
perr2=""

xbins=""
xmin=""
xmax=""

theta_min=""
theta_max=""
vz_min=""
vz_max=""
chi2_max=""
detector_pe_min=""
pe_min=""
lv_min=""
lw_min=""
pcal_min=""
sf_max=""
use_sf_cut=""
use_ecin_pcal_cut=""
ecin_slope=""
ecin_intercept=""
q2_min=""
q2_max=""
w_min=""
w_max=""

charge_rel_err=""
charge_branch=""
require_target_match=""
require_is_qe=""
require_fd_status=""
save_plots=""
write_diagnostics=""

print_usage(){
    cat <<USAGE
Usage:
  ./run_qe_analysis.sh [options]

Build/control options:
  --config <yaml>                 Read a simple flat YAML config file.
  --build                         Run make before executing.
  --exec <path>                   Analysis executable. Default: ./qe_analysis.
  --dry-run                       Print command without executing.
  --help, -h                      Show this help.

Required analysis options, unless provided by YAML:
  --run1 <run>                    Run number for run 1.
  --run2 <run>                    Run number for run 2.

Optional direct file options:
  --file1 <root>                  QE skim ROOT file for run 1.
  --file2 <root>                  QE skim ROOT file for run 2.

Auto path/output options:
  --skim-dir <dir>                Directory containing per-run skim ROOT files.
  --skim-prefix <prefix>          File prefix before _0NNNNN.root.
  --auto-inputs 0|1               If 1, derive file1/file2 from run1/run2. Default 1.
  --auto-output 0|1               If 1, derive output when output is not explicitly set. Default 1.
  --pair-subdir 0|1               If 1, put output in outputs/RUN1_RUN2/. Default 1.
  --output-dir-base <dir>         Base output directory. Default outputs.
  --output-stem <name>            Output filename stem. Default qe_Azz_approx.
  --output-tag <tag>              Output filename tag. Default: cut_profile.
  --cut-profile <kin0|kin1>         kin0=QE detector/PID only; kin1=detector/PID + tight kinematics.
  --profile-subdir 0|1              If 1, output under output_dir_base/cut_profile/. Default 1.

Tensor/vector polarization options:
  --q1 <Q1>                       Tensor polarization for run1 as fraction.
  --q2 <Q2>                       Tensor polarization for run2 as fraction.
  --qerr1 <dQ1>                   Absolute tensor-polarization uncertainty.
  --qerr2 <dQ2>                   Absolute tensor-polarization uncertainty.
  --p1 <P1>                       Vector polarization for run1 as fraction.
  --p2 <P2>                       Vector polarization for run2 as fraction.
  --perr1 <dP1>                   Vector-polarization uncertainty for run1.
  --perr2 <dP2>                   Vector-polarization uncertainty for run2.

Binning:
  --xbins <N>
  --xmin <x>
  --xmax <x>

Cuts:
  --theta-min <deg>
  --theta-max <deg>
  --vz-min <cm>
  --vz-max <cm>
  --chi2-max <value>
  --detector-pe-min <GeV>
  --pe-min <GeV>
  --lv-min <value>
  --lw-min <value>
  --pcal-min <GeV>
  --sf-max <value>
  --use-sf-cut 0|1
  --use-ecin-pcal-cut 0|1
  --ecin-slope <value>
  --ecin-intercept <value>
  --q2-min <GeV2>
  --q2-max <GeV2>
  --w-min <GeV>
  --w-max <GeV>

Charge/scaler options:
  --charge-branch <name>          Default in C++: fcupgated_delta.
  --charge-rel-err <value>        Example: 0.002 for 0.2%.

Analysis switches:
  --require-target-match 0|1
  --require-is-qe 0|1
  --require-fd-status 0|1
  --save-plots 0|1
  --write-diagnostics 0|1         Write diagnostic trees/CSVs. C++ default 1.
  --output <root>                 Output ROOT file.

Examples:
  QE check:
    ./run_qe_analysis.sh --build --config configs/qe_analysis.yaml

  Direct command-line run:
    ./run_qe_analysis.sh \\
      --file1 ../../../rgcskim/rootfiles/skim_sidisdvcs_inclusive_qe_eT0pid1kin1tf1/spring23/sidisdvcs_inclusive_qe_eT0pid1kin1tf1_017492.root \\
      --file2 ../../../rgcskim/rootfiles/skim_sidisdvcs_inclusive_qe_eT0pid1kin1tf1/spring23/sidisdvcs_inclusive_qe_eT0pid1kin1tf1_017792.root \\
      --run1 17492 --run2 17792 \\
      --p1 0.5120999 --p2 -0.262899 \\
      --theta-min 7.8 --theta-max 8.2 \\
      --q2-min 1.9433 --q2-max 2.0574 \\
      --w-min 0.0 --w-max 1.073 \\
      --output outputs/qe_R_17492_17792.root
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
    if [ -z "$CONFIG" ]; then
        return
    fi

    if [ ! -f "$CONFIG" ]; then
        echo "Error: config file not found: $CONFIG"
        exit 1
    fi

    set_from_yaml_if_present file1 file1
    set_from_yaml_if_present file2 file2
    set_from_yaml_if_present run1 run1
    set_from_yaml_if_present run2 run2

    local output_yaml=""
    output_yaml=$(yaml_get output "$CONFIG" || true)
    if [ -n "$output_yaml" ]; then
        output="$output_yaml"
        output_is_user_set=1
    fi

    set_from_yaml_if_present skim_dir skim_dir
    set_from_yaml_if_present skim_prefix skim_prefix
    set_from_yaml_if_present auto_inputs auto_inputs
    set_from_yaml_if_present auto_output auto_output
    set_from_yaml_if_present pair_subdir pair_subdir
    set_from_yaml_if_present output_dir_base output_dir_base
    set_from_yaml_if_present output_stem output_stem
    set_from_yaml_if_present output_tag output_tag
    set_from_yaml_if_present cut_profile cut_profile
    set_from_yaml_if_present profile_subdir profile_subdir

    set_from_yaml_if_present q1 q1
    set_from_yaml_if_present q2 q2
    set_from_yaml_if_present qerr1 qerr1
    set_from_yaml_if_present qerr2 qerr2

    set_from_yaml_if_present p1 p1
    set_from_yaml_if_present p2 p2
    set_from_yaml_if_present perr1 perr1
    set_from_yaml_if_present perr2 perr2

    set_from_yaml_if_present xbins xbins
    set_from_yaml_if_present xmin xmin
    set_from_yaml_if_present xmax xmax

    set_from_yaml_if_present theta_min theta_min
    set_from_yaml_if_present theta_max theta_max
    set_from_yaml_if_present vz_min vz_min
    set_from_yaml_if_present vz_max vz_max
    set_from_yaml_if_present chi2_max chi2_max
    set_from_yaml_if_present detector_pe_min detector_pe_min
    set_from_yaml_if_present pe_min pe_min
    set_from_yaml_if_present lv_min lv_min
    set_from_yaml_if_present lw_min lw_min
    set_from_yaml_if_present pcal_min pcal_min
    set_from_yaml_if_present sf_max sf_max
    set_from_yaml_if_present use_sf_cut use_sf_cut
    set_from_yaml_if_present use_ecin_pcal_cut use_ecin_pcal_cut
    set_from_yaml_if_present ecin_slope ecin_slope
    set_from_yaml_if_present ecin_intercept ecin_intercept
    set_from_yaml_if_present q2_min q2_min
    set_from_yaml_if_present q2_max q2_max
    set_from_yaml_if_present w_min w_min
    set_from_yaml_if_present w_max w_max

    set_from_yaml_if_present charge_rel_err charge_rel_err
    set_from_yaml_if_present charge_branch charge_branch
    set_from_yaml_if_present require_target_match require_target_match
    set_from_yaml_if_present require_is_qe require_is_qe
    set_from_yaml_if_present require_fd_status require_fd_status
    set_from_yaml_if_present save_plots save_plots
    set_from_yaml_if_present write_diagnostics write_diagnostics
}

#============================================================
# First pass: find --config, --build, --exec, --dry-run.
#============================================================

args=("$@")
i=0
while [ $i -lt ${#args[@]} ]; do
    a="${args[$i]}"
    case "$a" in
        --config)
            i=$((i+1))
            CONFIG="${args[$i]}"
            ;;
        --build)
            DO_BUILD=1
            ;;
        --exec)
            i=$((i+1))
            EXEC="${args[$i]}"
            ;;
        --dry-run)
            DRY_RUN=1
            ;;
        --help|-h)
            print_usage
            exit 0
            ;;
    esac
    i=$((i+1))
done

load_config

#============================================================
# Second pass: command line overrides YAML.
#============================================================

i=0
while [ $i -lt ${#args[@]} ]; do
    a="${args[$i]}"
    case "$a" in
        --config|--exec)
            i=$((i+2))
            continue
            ;;
        --build|--dry-run)
            i=$((i+1))
            continue
            ;;
        --file1) i=$((i+1)); file1="${args[$i]}" ;;
        --file2) i=$((i+1)); file2="${args[$i]}" ;;
        --run1) i=$((i+1)); run1="${args[$i]}" ;;
        --run2) i=$((i+1)); run2="${args[$i]}" ;;
        --output|-o) i=$((i+1)); output="${args[$i]}"; output_is_user_set=1 ;;

        --skim-dir) i=$((i+1)); skim_dir="${args[$i]}" ;;
        --skim-prefix) i=$((i+1)); skim_prefix="${args[$i]}" ;;
        --auto-inputs) i=$((i+1)); auto_inputs="${args[$i]}" ;;
        --auto-output) i=$((i+1)); auto_output="${args[$i]}" ;;
        --pair-subdir) i=$((i+1)); pair_subdir="${args[$i]}" ;;
        --output-dir-base) i=$((i+1)); output_dir_base="${args[$i]}" ;;
        --output-stem) i=$((i+1)); output_stem="${args[$i]}" ;;
        --output-tag) i=$((i+1)); output_tag="${args[$i]}" ;;
        --cut-profile) i=$((i+1)); cut_profile="${args[$i]}" ;;
        --profile-subdir) i=$((i+1)); profile_subdir="${args[$i]}" ;;

        --q1) i=$((i+1)); q1="${args[$i]}" ;;
        --q2) i=$((i+1)); q2="${args[$i]}" ;;
        --qerr1) i=$((i+1)); qerr1="${args[$i]}" ;;
        --qerr2) i=$((i+1)); qerr2="${args[$i]}" ;;
        --p1) i=$((i+1)); p1="${args[$i]}" ;;
        --p2) i=$((i+1)); p2="${args[$i]}" ;;
        --perr1) i=$((i+1)); perr1="${args[$i]}" ;;
        --perr2) i=$((i+1)); perr2="${args[$i]}" ;;

        --xbins) i=$((i+1)); xbins="${args[$i]}" ;;
        --xmin) i=$((i+1)); xmin="${args[$i]}" ;;
        --xmax) i=$((i+1)); xmax="${args[$i]}" ;;

        --theta-min) i=$((i+1)); theta_min="${args[$i]}" ;;
        --theta-max) i=$((i+1)); theta_max="${args[$i]}" ;;
        --vz-min) i=$((i+1)); vz_min="${args[$i]}" ;;
        --vz-max) i=$((i+1)); vz_max="${args[$i]}" ;;
        --chi2-max) i=$((i+1)); chi2_max="${args[$i]}" ;;
        --detector-pe-min) i=$((i+1)); detector_pe_min="${args[$i]}" ;;
        --pe-min) i=$((i+1)); pe_min="${args[$i]}" ;;
        --lv-min) i=$((i+1)); lv_min="${args[$i]}" ;;
        --lw-min) i=$((i+1)); lw_min="${args[$i]}" ;;
        --pcal-min) i=$((i+1)); pcal_min="${args[$i]}" ;;
        --sf-max) i=$((i+1)); sf_max="${args[$i]}" ;;
        --use-sf-cut) i=$((i+1)); use_sf_cut="${args[$i]}" ;;
        --use-ecin-pcal-cut) i=$((i+1)); use_ecin_pcal_cut="${args[$i]}" ;;
        --ecin-slope) i=$((i+1)); ecin_slope="${args[$i]}" ;;
        --ecin-intercept) i=$((i+1)); ecin_intercept="${args[$i]}" ;;
        --q2-min) i=$((i+1)); q2_min="${args[$i]}" ;;
        --q2-max) i=$((i+1)); q2_max="${args[$i]}" ;;
        --w-min) i=$((i+1)); w_min="${args[$i]}" ;;
        --w-max) i=$((i+1)); w_max="${args[$i]}" ;;

        --charge-rel-err) i=$((i+1)); charge_rel_err="${args[$i]}" ;;
        --charge-branch) i=$((i+1)); charge_branch="${args[$i]}" ;;
        --require-target-match) i=$((i+1)); require_target_match="${args[$i]}" ;;
        --require-is-qe) i=$((i+1)); require_is_qe="${args[$i]}" ;;
        --require-fd-status) i=$((i+1)); require_fd_status="${args[$i]}" ;;
        --save-plots) i=$((i+1)); save_plots="${args[$i]}" ;;
        --write-diagnostics) i=$((i+1)); write_diagnostics="${args[$i]}" ;;

        --help|-h)
            print_usage
            exit 0
            ;;
        *)
            echo "Error: unknown option: $a"
            echo
            print_usage
            exit 1
            ;;
    esac
    i=$((i+1))
done

bool01(){
    local v="$1"
    case "$v" in
        1|true|TRUE|yes|YES) return 0 ;;
        0|false|FALSE|no|NO) return 1 ;;
        *) echo "Error: expected 0/1 for boolean setting, got: $v"; exit 1 ;;
    esac
}

run_padded(){
    printf "%06d" "$1"
}

infer_skim_prefix(){
    if [ -n "$skim_prefix" ]; then
        return
    fi

    if [ -z "$skim_dir" ]; then
        return
    fi

    local d="$skim_dir"
    local cand=""
    local n=0
    while [ $n -lt 4 ]; do
        cand=$(basename "$d")
        if [[ "$cand" == skim_* ]]; then
            skim_prefix="${cand#skim_}"
            return
        fi
        d=$(dirname "$d")
        n=$((n+1))
    done
}

derive_output_tag(){
    if [ -n "$output_tag" ]; then
        printf '%s' "$output_tag"
        return
    fi

    if [ -n "$cut_profile" ]; then
        printf '%s' "$cut_profile"
        return
    fi

    if [[ "$skim_prefix" =~ kin([0-9]+) ]]; then
        printf 'kin%s' "${BASH_REMATCH[1]}"
        return
    fi

    printf 'qe'
}

validate_cut_profile(){
    case "$cut_profile" in
        kin0|kin1) ;;
        *) echo "Error: cut_profile must be kin0 or kin1, got: $cut_profile"; exit 1 ;;
    esac
}

autofill_files_and_output(){
    if [ -z "$run1" ] || [ -z "$run2" ]; then
        echo "Error: --run1 and --run2 are required, either directly or through --config."
        exit 1
    fi

    infer_skim_prefix

    if bool01 "$auto_inputs"; then
        if [ -z "$file1" ] || [ -z "$file2" ]; then
            if [ -z "$skim_dir" ]; then
                echo "Error: file1/file2 are not set and skim_dir is empty."
                echo "       Either set file1/file2 explicitly, or set skim_dir and skim_prefix."
                exit 1
            fi
            if [ -z "$skim_prefix" ]; then
                echo "Error: file1/file2 are not set and skim_prefix could not be inferred."
                echo "       Set skim_prefix, for example: sidisdvcs_inclusive_qe_eT0pid1kin0tf1"
                exit 1
            fi
        fi

        if [ -z "$file1" ]; then
            file1="${skim_dir%/}/${skim_prefix}_$(run_padded "$run1").root"
        fi
        if [ -z "$file2" ]; then
            file2="${skim_dir%/}/${skim_prefix}_$(run_padded "$run2").root"
        fi
    fi

    if bool01 "$auto_output" && [ "$output_is_user_set" -eq 0 ]; then
        local tag=""
        tag=$(derive_output_tag)
        local pair="${run1}_${run2}"
        local outdir="${output_dir_base%/}"
        if bool01 "$profile_subdir"; then
            outdir="${outdir}/${tag}"
        fi
        if bool01 "$pair_subdir"; then
            outdir="${outdir}/${pair}"
        fi
        output="${outdir}/${output_stem}_${pair}_${tag}.root"
    fi
}

validate_cut_profile
autofill_files_and_output

if [ -z "$file1" ] || [ -z "$file2" ] || [ -z "$run1" ] || [ -z "$run2" ]; then
    echo "Error: --file1, --file2, --run1, and --run2 are required after auto resolution."
    exit 1
fi

if [ ! -f "$file1" ]; then
    echo "Error: file1 not found: $file1"
    exit 1
fi

if [ ! -f "$file2" ]; then
    echo "Error: file2 not found: $file2"
    exit 1
fi

if [ "$DO_BUILD" -eq 1 ]; then
    echo "Building qe_analysis with make..."
    make -C "$SCRIPT_DIR"
fi

if [ ! -x "$EXEC" ]; then
    echo "Error: executable not found or not executable: $EXEC"
    echo "Run: make"
    exit 1
fi

mkdir -p "$(dirname "$output")"

CMD=("$EXEC")
CMD+=(--file1 "$file1")
CMD+=(--file2 "$file2")
CMD+=(--run1 "$run1")
CMD+=(--run2 "$run2")
CMD+=(--output "$output")
CMD+=(--cut-profile "$cut_profile")

add_opt(){
    local opt="$1"
    local val="$2"
    if [ -n "$val" ]; then
        CMD+=("$opt" "$val")
    fi
}

add_opt --q1 "$q1"
add_opt --q2 "$q2"
add_opt --qerr1 "$qerr1"
add_opt --qerr2 "$qerr2"
add_opt --p1 "$p1"
add_opt --p2 "$p2"
add_opt --perr1 "$perr1"
add_opt --perr2 "$perr2"

add_opt --xbins "$xbins"
add_opt --xmin "$xmin"
add_opt --xmax "$xmax"

add_opt --theta-min "$theta_min"
add_opt --theta-max "$theta_max"
add_opt --vz-min "$vz_min"
add_opt --vz-max "$vz_max"
add_opt --chi2-max "$chi2_max"
add_opt --detector-pe-min "$detector_pe_min"
add_opt --pe-min "$pe_min"
add_opt --lv-min "$lv_min"
add_opt --lw-min "$lw_min"
add_opt --pcal-min "$pcal_min"
add_opt --sf-max "$sf_max"
add_opt --use-sf-cut "$use_sf_cut"
add_opt --use-ecin-pcal-cut "$use_ecin_pcal_cut"
add_opt --ecin-slope "$ecin_slope"
add_opt --ecin-intercept "$ecin_intercept"
add_opt --q2-min "$q2_min"
add_opt --q2-max "$q2_max"
add_opt --w-min "$w_min"
add_opt --w-max "$w_max"

add_opt --charge-rel-err "$charge_rel_err"
add_opt --charge-branch "$charge_branch"
add_opt --require-target-match "$require_target_match"
add_opt --require-is-qe "$require_is_qe"
add_opt --require-fd-status "$require_fd_status"
add_opt --save-plots "$save_plots"
add_opt --write-diagnostics "$write_diagnostics"

echo "================================================================================"
echo "QE tensor two-run analysis"
echo "================================================================================"
echo "Working dir:   $(pwd)"
echo "Script dir:    $SCRIPT_DIR"
echo "Config:        ${CONFIG:-none}"
echo "Executable:    $EXEC"
echo "Auto inputs:   $auto_inputs"
echo "Auto output:   $auto_output"
echo "Skim dir:      ${skim_dir:-none}"
echo "Skim prefix:   ${skim_prefix:-none}"
echo "Cut profile:   $cut_profile"
echo "Profile subdir: $profile_subdir"
echo "Run 1 file:    $file1"
echo "Run 2 file:    $file2"
echo "Run 1:         $run1"
echo "Run 2:         $run2"
echo "Output:        $output"
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