#!/bin/bash

# RGC skim runner.
# Runs rgcskim over one HIPO file or one folder of HIPO files.
# Supported modes:
#   sidis      : electron + selected hadron
#   inclusive  : electron only, region = all, dis, res, qe
#   dihadron   : electron + pi+ + pi-
# Common examples:
#   ./run_all_hipo.sh --mode inclusive --region qe /path/to/sidisdvcs --detpidcut 1 --kincut 1
#   ./run_all_hipo.sh --mode sidis --pid 211 /path/to/sidisdvcs --detpidcut 1
# Useful options:
#   --detpidcut 0|1 --kincut 0|1
#   --targetmap <csv> --polsource offline|online --targetfilter all|matched
#   --period auto|none|summer22|fall22|spring23
#   --test-one-run
#   --outformat root|hipo|both
#   --run <run> | --runs <r1,r2,r3> | --runrange <start:end> | --runlist <file.txt>
#   --jobtag <tag>

print_usage(){

    echo "Usage:"
    echo "  $0 [--dataset DATASET_TAG] --mode sidis --pid <hadron_pid> <hipo_folder_or_file> [options]"
    echo "  $0 [--dataset DATASET_TAG] --mode inclusive --region all|dis|res|qe <hipo_folder_or_file> [options]"
    echo "  $0 [--dataset DATASET_TAG] --mode dihadron <hipo_folder_or_file> [options]"
    echo ""
    echo "Required:"
    echo "  --mode sidis|inclusive|dihadron"
    echo "  --pid <hadron_pid>              Required only for --mode sidis"
    echo "  --region all|dis|res|qe         Used only for --mode inclusive. Default = all"
    echo "  <hipo_folder_or_file>           Folder containing .hipo files, or one .hipo file"
    echo ""
    echo "Options:"
    echo "  --dataset <tag>                 Override dataset tag. Default = basename of HIPO folder"
    echo "  --electrontree 0|1              SIDIS diagnostic electron tree. Default = 0"
    echo "  --detpidcut 0|1                 Apply detector/PID cuts. Default = 0"
    echo "  --kincut 0|1                    Apply dedicated region kinematic cuts. Default = 0"
    echo "  --targetmap <csv>               Optional run-level target/polarization CSV"
    echo "  --polsource offline|online      Polarization source for target_pol/vector_pol. Default = offline"
    echo "  --targetfilter all|matched      all = process all files, matched = process target-map runs only"
    echo "  --targetfilter 0|1              Alias: 0 = all, 1 = matched"
    echo "  --period auto|none|summer22|fall22|spring23"
    echo "                                  Optional data-period output subfolder. Default = none"
    echo "  --test-one-run                 Process the first HIPO file that passes file-level filters"
    echo "  --outformat root|hipo|both     Output format. Default = root"
    echo "  --run <run>                     Process one run only"
    echo "  --runs <r1,r2,r3>               Process comma-separated run list"
    echo "  --runrange <start:end>          Process inclusive run range. start-end also accepted"
    echo "  --runlist <file.txt>            Process runs from text file, one run per line"
    echo "  --jobtag <tag>                  Optional extra output subfolder"
    echo ""
    echo "Hadron PID:"
    echo "   211   = e pi+"
    echo "  -211   = e pi-"
    echo "   321   = e K+"
    echo "  -321   = e K-"
    echo "   2212  = e proton"
    echo ""
    echo "Examples:"
    echo ""
    echo "  Inclusive QE with target map and automatic period subfolder:"
    echo "    $0 --mode inclusive --region qe /path/to/sidisdvcs --detpidcut 1 --kincut 1 --targetmap target.csv --polsource offline --targetfilter matched --period auto"
    echo ""
    echo "  Inclusive QE test using first target-matched run from a folder:"
    echo "    $0 --mode inclusive --region qe /path/to/sidisdvcs --detpidcut 1 --kincut 1 --targetmap target.csv --polsource offline --targetfilter matched --period auto --test-one-run --jobtag test_one_run"
    echo ""
    echo "  Inclusive QE test writing ROOT and filtered HIPO output:"
    echo "    $0 --mode inclusive --region qe /path/to/sidisdvcs --detpidcut 1 --kincut 1 --targetmap target.csv --polsource offline --targetfilter matched --period auto --test-one-run --outformat both --jobtag test_one_run"
    echo ""
    echo "  Inclusive QE for one run only:"
    echo "    $0 --mode inclusive --region qe /path/to/sidisdvcs --detpidcut 1 --kincut 1 --targetmap target.csv --polsource offline --targetfilter matched --period auto --run 16270 --jobtag test_16270"
    echo ""
    echo "  SIDIS pi-plus:"
    echo "    $0 --mode sidis --pid 211 /path/to/sidisdvcs --detpidcut 1"
    echo ""
    echo "  Dihadron pi+pi-:"
    echo "    $0 --mode dihadron /path/to/sidisdvcs --detpidcut 1"
}

# ------------------------------------------------------------
# Defaults
# ------------------------------------------------------------

EXE=./rgcskim

MODE=""
HADRON_PID=""
HIPO_INPUT=""

# SIDIS only optional diagnostic electron tree.
ELECTRON_TREE=0
ELECTRON_TREE_USER_SET=0

# Detector/PID cut switch passed to rgcskim.
DET_PID_CUT=0

# Kinematic cut switch passed to rgcskim.
KIN_CUT=0

# Inclusive region selector. Used only for --mode inclusive.
REGION="all"

# Optional target/polarization map.
TARGET_MAP=""
POL_SOURCE="offline"
TARGET_FILTER="all"

# Dataset tag. If empty, it is inferred from the HIPO folder name.
DATASET_TAG=""

# Optional data period output subfolder.
PERIOD_MODE="none"
DATA_PERIOD=""

# Optional run selection controls.
RUN_TEST_ONE_RUN=0
RUN_SINGLE=""
RUNS_CSV=""
RUN_RANGE=""
RUNLIST_FILE=""
RUN_SELECTION_MODE="all"
declare -a SELECTED_RUNS
SELECTED_RUNS=()
RUN_RANGE_START=""
RUN_RANGE_END=""

# Optional extra output subfolder.
JOB_TAG=""

# Optional output base path
OUT_PATH=""

# Output format.
OUTPUT_FORMAT="root"

ROOT_BASE_DIR="rootfiles"
HIPO_BASE_DIR="hipofiles"
LOG_BASE_DIR="logs"

# Job timing.
JOB_START_EPOCH=$(date +%s)
JOB_START_TIME=$(date '+%Y-%m-%d %H:%M:%S')

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

normalize_run_number(){

    local raw="$1"

    # Remove whitespace.
    raw=$(echo "$raw" | tr -d '[:space:]')

    if [ -z "$raw" ]; then
        return 1
    fi

    # Only allow digits.
    if ! [[ "$raw" =~ ^[0-9]+$ ]]; then
        return 1
    fi

    # Convert 016270 -> 16270 safely.
    echo $((10#$raw))
    return 0
}

run_passes_selection(){

    local run_num="$1"

    if [ "$RUN_SELECTION_MODE" = "all" ]; then
        return 0
    fi

    if [ "$RUN_SELECTION_MODE" = "test_one_run" ]; then
        return 0
    fi

    if [ "$RUN_SELECTION_MODE" = "range" ]; then
        if [ "$run_num" -ge "$RUN_RANGE_START" ] && [ "$run_num" -le "$RUN_RANGE_END" ]; then
            return 0
        fi
        return 1
    fi

    # single, list, and file modes are stored in SELECTED_RUNS.
    local allowed_run

    for allowed_run in "${SELECTED_RUNS[@]}"
    do
        if [ "$run_num" -eq "$allowed_run" ]; then
            return 0
        fi
    done

    return 1
}

detect_period_from_path(){

    local input_path="$1"
    local lower_input_path

    lower_input_path=$(echo "$input_path" | tr '[:upper:]' '[:lower:]')

    case "$lower_input_path" in
        *summer22*)
            echo "summer22"
            ;;
        *fall22*)
            echo "fall22"
            ;;
        *spring23*)
            echo "spring23"
            ;;
        *)
            echo "unknown"
            ;;
    esac
}

sanitize_folder_tag(){

    local raw="$1"

    raw=$(echo "$raw" | sed 's/[^A-Za-z0-9._-]/_/g')
    echo "$raw"
}

format_elapsed_time(){

    local total_seconds="$1"

    printf "%02d:%02d:%02d" \
        $((total_seconds / 3600)) \
        $(((total_seconds % 3600) / 60)) \
        $((total_seconds % 60))
}

# ------------------------------------------------------------
# Parse command line arguments
# ------------------------------------------------------------

if [ $# -lt 2 ]; then
    print_usage
    exit 1
fi

while [ $# -gt 0 ]; do

    case "$1" in

        --dataset|--sample)
            if [ $# -lt 2 ]; then
                echo "Error: $1 requires a dataset tag."
                exit 1
            fi
            DATASET_TAG="$2"
            shift 2
            ;;

        --period|--dataperiod)
            if [ $# -lt 2 ]; then
                echo "Error: $1 requires auto, none, summer22, fall22, or spring23."
                exit 1
            fi
            PERIOD_MODE="$2"
            shift 2
            ;;

        --mode)
            if [ $# -lt 2 ]; then
                echo "Error: --mode requires sidis, inclusive, or dihadron."
                exit 1
            fi
            MODE="$2"
            shift 2
            ;;

        --pid|--hadronpid)
            if [ $# -lt 2 ]; then
                echo "Error: $1 requires a hadron PID."
                exit 1
            fi
            HADRON_PID="$2"
            shift 2
            ;;

        --region)
            if [ $# -lt 2 ]; then
                echo "Error: --region requires all, dis, res, resonance, qe, quasi, quasielastic, or quasi-elastic."
                exit 1
            fi
            REGION="$2"
            shift 2
            ;;

        --electronTree|--electrontree)
            if [ $# -lt 2 ]; then
                echo "Error: $1 requires 0 or 1."
                exit 1
            fi
            ELECTRON_TREE="$2"
            ELECTRON_TREE_USER_SET=1
            shift 2
            ;;

        --detpidcut)
            if [ $# -lt 2 ]; then
                echo "Error: --detpidcut requires 0 or 1."
                exit 1
            fi
            DET_PID_CUT="$2"
            shift 2
            ;;

        --kincut)
            if [ $# -lt 2 ]; then
                echo "Error: --kincut requires 0 or 1."
                exit 1
            fi
            KIN_CUT="$2"
            shift 2
            ;;

        --targetmap)
            if [ $# -lt 2 ]; then
                echo "Error: --targetmap requires a CSV file path."
                exit 1
            fi
            TARGET_MAP="$2"
            shift 2
            ;;

        --polsource)
            if [ $# -lt 2 ]; then
                echo "Error: --polsource requires offline or online."
                exit 1
            fi
            POL_SOURCE="$2"
            shift 2
            ;;

        --targetfilter)
            if [ $# -lt 2 ]; then
                echo "Error: --targetfilter requires all, matched, 0, or 1."
                exit 1
            fi
            TARGET_FILTER="$2"
            shift 2
            ;;

        --test-one-run|--test_one_run|--testonerun)
            RUN_TEST_ONE_RUN=1
            shift 1
            ;;

        --outformat|--outputformat)
            if [ $# -lt 2 ]; then
                echo "Error: $1 requires root, hipo, or both."
                exit 1
            fi
            OUTPUT_FORMAT="$2"
            shift 2
            ;;

        --run)
            if [ $# -lt 2 ]; then
                echo "Error: --run requires one run number."
                exit 1
            fi
            RUN_SINGLE="$2"
            shift 2
            ;;

        --runs)
            if [ $# -lt 2 ]; then
                echo "Error: --runs requires a comma-separated run list."
                exit 1
            fi
            RUNS_CSV="$2"
            shift 2
            ;;

        --runrange)
            if [ $# -lt 2 ]; then
                echo "Error: --runrange requires start:end."
                exit 1
            fi
            RUN_RANGE="$2"
            shift 2
            ;;

        --runlist)
            if [ $# -lt 2 ]; then
                echo "Error: --runlist requires a text file."
                exit 1
            fi
            RUNLIST_FILE="$2"
            shift 2
            ;;

        --jobtag|--label)
            if [ $# -lt 2 ]; then
                echo "Error: $1 requires a tag name."
                exit 1
            fi
            JOB_TAG="$2"
            shift 2
            ;;
            
        --out-path)
            if [ $# -lt 2 ]; then
                echo "Error: --out-path requires a base path for the output."
                exit 1
            fi
            OUT_PATH="$2"
            shift 2
            ;;

        --help|-h)
            print_usage
            exit 0
            ;;

        --*)
            echo "Error: unknown option: $1"
            echo ""
            print_usage
            exit 1
            ;;

        *)
            if [ -z "$HIPO_INPUT" ]; then
                HIPO_INPUT="$1"
                shift 1
            else
                echo "Error: extra positional argument: $1"
                echo ""
                print_usage
                exit 1
            fi
            ;;
    esac

done

# ------------------------------------------------------------
# Validate mode and required inputs
# ------------------------------------------------------------

if [ -z "$MODE" ]; then
    echo "Error: --mode is required."
    echo ""
    print_usage
    exit 1
fi

if [ "$MODE" != "sidis" ] && [ "$MODE" != "inclusive" ] && [ "$MODE" != "dihadron" ]; then
    echo "Error: --mode must be sidis, inclusive, or dihadron."
    exit 1
fi

if [ "$MODE" = "sidis" ] && [ -z "$HADRON_PID" ]; then
    echo "Error: SIDIS mode requires --pid."
    echo ""
    echo "Example:"
    echo "  $0 --mode sidis --pid 211 /path/to/sidisdvcs --detpidcut 1"
    exit 1
fi

if [ "$MODE" = "dihadron" ] && [ -n "$HADRON_PID" ]; then
    echo "Error: --pid should not be used with --mode dihadron."
    echo "Dihadron mode is fixed to e pi+ pi- for this version."
    exit 1
fi

if [ "$MODE" = "inclusive" ] && [ -n "$HADRON_PID" ]; then
    echo "Error: --pid should not be used with --mode inclusive."
    exit 1
fi

if [ -z "$HIPO_INPUT" ]; then
    echo "Error: no HIPO folder or HIPO file was provided."
    echo ""
    print_usage
    exit 1
fi

if [ ! -x "$EXE" ]; then
    echo "Error: executable $EXE not found or not executable."
    echo "Compile first, for example:"
    echo "  make"
    exit 1
fi

if [ "$ELECTRON_TREE" != "0" ] && [ "$ELECTRON_TREE" != "1" ]; then
    echo "Error: --electrontree must be 0 or 1."
    exit 1
fi

if [ "$DET_PID_CUT" != "0" ] && [ "$DET_PID_CUT" != "1" ]; then
    echo "Error: --detpidcut must be 0 or 1."
    exit 1
fi

if [ "$KIN_CUT" != "0" ] && [ "$KIN_CUT" != "1" ]; then
    echo "Error: --kincut must be 0 or 1."
    exit 1
fi

if [ "$POL_SOURCE" != "offline" ] && [ "$POL_SOURCE" != "online" ]; then
    echo "Error: --polsource must be offline or online."
    exit 1
fi

# ------------------------------------------------------------
# Normalize and validate output format
# ------------------------------------------------------------

OUTPUT_FORMAT=$(echo "$OUTPUT_FORMAT" | tr '[:upper:]' '[:lower:]')

if [ "$OUTPUT_FORMAT" != "root" ] && [ "$OUTPUT_FORMAT" != "hipo" ] && [ "$OUTPUT_FORMAT" != "both" ]; then
    echo "Error: --outformat must be root, hipo, or both."
    exit 1
fi

# ------------------------------------------------------------
# Normalize and validate target filter
# ------------------------------------------------------------

if [ "$TARGET_FILTER" = "0" ]; then
    TARGET_FILTER="all"
fi

if [ "$TARGET_FILTER" = "1" ]; then
    TARGET_FILTER="matched"
fi

if [ "$TARGET_FILTER" != "all" ] && [ "$TARGET_FILTER" != "matched" ]; then
    echo "Error: --targetfilter must be all, matched, 0, or 1."
    exit 1
fi

if [ "$TARGET_FILTER" = "matched" ] && [ -z "$TARGET_MAP" ]; then
    echo "Error: --targetfilter matched requires --targetmap <csv>."
    exit 1
fi

if [ -n "$TARGET_MAP" ] && [ ! -f "$TARGET_MAP" ]; then
    echo "Error: target map file does not exist:"
    echo "  $TARGET_MAP"
    exit 1
fi

# ------------------------------------------------------------
# Normalize and validate inclusive region
# ------------------------------------------------------------

REGION=$(echo "$REGION" | tr '[:upper:]' '[:lower:]')

if [ "$REGION" = "resonance" ]; then
    REGION="res"
fi

if [ "$REGION" = "quasi" ] || [ "$REGION" = "quasielastic" ] || [ "$REGION" = "quasi-elastic" ]; then
    REGION="qe"
fi

if [ "$REGION" != "all" ] && [ "$REGION" != "dis" ] && [ "$REGION" != "res" ] && [ "$REGION" != "qe" ]; then
    echo "Error: --region must be all, dis, res, resonance, qe, quasi, quasielastic, or quasi-elastic."
    exit 1
fi

if [ "$MODE" != "inclusive" ] && [ "$REGION" != "all" ]; then
    echo "Error: --region is only meaningful for --mode inclusive."
    echo "For SIDIS or dihadron mode, do not pass --region."
    exit 1
fi

# Inclusive mode does not use the extra electron diagnostic tree.
if [ "$MODE" = "inclusive" ] && [ "$ELECTRON_TREE_USER_SET" = "1" ] && [ "$ELECTRON_TREE" = "1" ]; then
    echo "Error: --electrontree 1 is only meaningful for --mode sidis."
    echo "Inclusive mode is already an electron-only tree."
    exit 1
fi

# Force inclusive and dihadron output names to eT0.
if [ "$MODE" = "inclusive" ] || [ "$MODE" = "dihadron" ]; then
    ELECTRON_TREE=0
fi

# ------------------------------------------------------------
# Normalize and validate data period mode
# ------------------------------------------------------------

PERIOD_MODE=$(echo "$PERIOD_MODE" | tr '[:upper:]' '[:lower:]')

if [ "$PERIOD_MODE" = "summer" ]; then
    PERIOD_MODE="summer22"
fi

if [ "$PERIOD_MODE" = "fall" ]; then
    PERIOD_MODE="fall22"
fi

if [ "$PERIOD_MODE" = "spring" ]; then
    PERIOD_MODE="spring23"
fi

if [ "$PERIOD_MODE" != "none" ] && \
   [ "$PERIOD_MODE" != "auto" ] && \
   [ "$PERIOD_MODE" != "summer22" ] && \
   [ "$PERIOD_MODE" != "fall22" ] && \
   [ "$PERIOD_MODE" != "spring23" ]; then

    echo "Error: --period must be auto, none, summer22, fall22, or spring23."
    exit 1
fi

if [ "$PERIOD_MODE" = "none" ]; then

    DATA_PERIOD=""

elif [ "$PERIOD_MODE" = "auto" ]; then

    DATA_PERIOD=$(detect_period_from_path "$HIPO_INPUT")

    if [ "$DATA_PERIOD" = "unknown" ]; then
        echo "Warning: --period auto could not infer summer22, fall22, or spring23 from input path."
        echo "Warning: using old output structure without a period subfolder."
        DATA_PERIOD=""
    fi

else

    DATA_PERIOD="$PERIOD_MODE"

fi

# ------------------------------------------------------------
# Validate and build run selection list
# ------------------------------------------------------------

N_RUN_FILTERS=0

if [ "$RUN_TEST_ONE_RUN" = "1" ]; then
    N_RUN_FILTERS=$((N_RUN_FILTERS + 1))
fi

if [ -n "$RUN_SINGLE" ]; then
    N_RUN_FILTERS=$((N_RUN_FILTERS + 1))
fi

if [ -n "$RUNS_CSV" ]; then
    N_RUN_FILTERS=$((N_RUN_FILTERS + 1))
fi

if [ -n "$RUN_RANGE" ]; then
    N_RUN_FILTERS=$((N_RUN_FILTERS + 1))
fi

if [ -n "$RUNLIST_FILE" ]; then
    N_RUN_FILTERS=$((N_RUN_FILTERS + 1))
fi

if [ "$N_RUN_FILTERS" -gt 1 ]; then
    echo "Error: use only one run-selection option at a time."
    echo "Allowed choices:"
    echo "  --test-one-run"
    echo "  --run"
    echo "  --runs"
    echo "  --runrange"
    echo "  --runlist"
    exit 1
fi

if [ "$RUN_TEST_ONE_RUN" = "1" ]; then

    RUN_SELECTION_MODE="test_one_run"

elif [ -n "$RUN_SINGLE" ]; then

    RUN_SELECTION_MODE="single"

    RUN_NORM=$(normalize_run_number "$RUN_SINGLE")
    if [ $? -ne 0 ]; then
        echo "Error: invalid run number for --run: $RUN_SINGLE"
        exit 1
    fi

    SELECTED_RUNS+=("$RUN_NORM")

elif [ -n "$RUNS_CSV" ]; then

    RUN_SELECTION_MODE="list"

    IFS=',' read -ra RAW_RUNS <<< "$RUNS_CSV"

    if [ "${#RAW_RUNS[@]}" -eq 0 ]; then
        echo "Error: --runs list is empty."
        exit 1
    fi

    for raw_run in "${RAW_RUNS[@]}"
    do
        RUN_NORM=$(normalize_run_number "$raw_run")
        if [ $? -ne 0 ]; then
            echo "Error: invalid run number in --runs: $raw_run"
            exit 1
        fi
        SELECTED_RUNS+=("$RUN_NORM")
    done

elif [ -n "$RUN_RANGE" ]; then

    RUN_SELECTION_MODE="range"

    RUN_RANGE_CLEAN=$(echo "$RUN_RANGE" | tr '-' ':')

    IFS=':' read -r raw_start raw_end extra_field <<< "$RUN_RANGE_CLEAN"

    if [ -n "$extra_field" ]; then
        echo "Error: invalid --runrange format: $RUN_RANGE"
        echo "Use start:end, for example 16270:16320."
        exit 1
    fi

    RUN_RANGE_START=$(normalize_run_number "$raw_start")
    if [ $? -ne 0 ]; then
        echo "Error: invalid runrange start: $raw_start"
        exit 1
    fi

    RUN_RANGE_END=$(normalize_run_number "$raw_end")
    if [ $? -ne 0 ]; then
        echo "Error: invalid runrange end: $raw_end"
        exit 1
    fi

    if [ "$RUN_RANGE_START" -gt "$RUN_RANGE_END" ]; then
        echo "Error: --runrange start is greater than end:"
        echo "  $RUN_RANGE_START > $RUN_RANGE_END"
        exit 1
    fi

elif [ -n "$RUNLIST_FILE" ]; then

    RUN_SELECTION_MODE="file"

    if [ ! -f "$RUNLIST_FILE" ]; then
        echo "Error: run-list file does not exist:"
        echo "  $RUNLIST_FILE"
        exit 1
    fi

    while IFS= read -r line
    do
        # Remove inline comments and whitespace.
        line="${line%%#*}"
        line=$(echo "$line" | tr -d '[:space:]')

        if [ -z "$line" ]; then
            continue
        fi

        RUN_NORM=$(normalize_run_number "$line")
        if [ $? -ne 0 ]; then
            echo "Error: invalid run number in run-list file:"
            echo "  file: $RUNLIST_FILE"
            echo "  line: $line"
            exit 1
        fi

        SELECTED_RUNS+=("$RUN_NORM")

    done < "$RUNLIST_FILE"

    if [ "${#SELECTED_RUNS[@]}" -eq 0 ]; then
        echo "Error: run-list file contains no valid runs:"
        echo "  $RUNLIST_FILE"
        exit 1
    fi

else

    RUN_SELECTION_MODE="all"

fi

# ------------------------------------------------------------
# Sanitize optional job tag
# ------------------------------------------------------------

if [ -n "$JOB_TAG" ]; then
    JOB_TAG=$(sanitize_folder_tag "$JOB_TAG")

    if [ -z "$JOB_TAG" ]; then
        echo "Error: --jobtag became empty after sanitizing."
        exit 1
    fi
fi

# ------------------------------------------------------------
# Target filter tag
# ------------------------------------------------------------

if [ "$TARGET_FILTER" = "all" ]; then
    TARGET_FILTER_TAG="tf0"
else
    TARGET_FILTER_TAG="tf1"
fi

# ------------------------------------------------------------
# Optional target map arguments passed to rgcskim
# ------------------------------------------------------------

TARGET_ARGS=()

if [ -n "$TARGET_MAP" ]; then
    TARGET_ARGS+=(--targetmap "$TARGET_MAP")
    TARGET_ARGS+=(--polsource "$POL_SOURCE")
fi

# ------------------------------------------------------------
# Load target map run list for --targetfilter matched
# ------------------------------------------------------------
# The script uses the first CSV column, run, to build a lookup table.
# This filtering happens before rgcskim is launched, using the 6 digit
# run number extracted from the HIPO filename:
#   sidisdvcs_016270.hipo -> RUN_TAG=016270 -> RUN_NUM=16270
# rgcskim itself still uses RUN::config.run internally to fill target metadata.
# ------------------------------------------------------------

declare -A TARGET_RUNS

TARGET_MAP_RUNS=0
TARGET_MAP_SKIPPED_LINES=0

if [ -n "$TARGET_MAP" ]; then

    while IFS=, read -r CSV_RUN CSV_START CSV_STOP CSV_SPECIES CSV_CELL CSV_ONLINE CSV_OFFLINE CSV_OFFLINE_ERR CSV_DOSE CSV_REST
    do
        CSV_RUN=$(echo "$CSV_RUN" | tr -d '[:space:]')

        if [ -z "$CSV_RUN" ]; then
            continue
        fi

        case "$CSV_RUN" in
            \#*)
                continue
                ;;
        esac

        if [[ "$CSV_RUN" =~ ^[0-9]+$ ]]; then
            CSV_RUN_NUM=$((10#$CSV_RUN))
            TARGET_RUNS[$CSV_RUN_NUM]=1
            TARGET_MAP_RUNS=$((TARGET_MAP_RUNS + 1))
        else
            TARGET_MAP_SKIPPED_LINES=$((TARGET_MAP_SKIPPED_LINES + 1))
        fi

    done < "$TARGET_MAP"

    if [ "$TARGET_FILTER" = "matched" ] && [ "$TARGET_MAP_RUNS" -eq 0 ]; then
        echo "Error: --targetfilter matched was requested, but no runs were loaded from:"
        echo "  $TARGET_MAP"
        exit 1
    fi

fi

# ------------------------------------------------------------
# Build HIPO file list
# ------------------------------------------------------------

HIPO_FILES=()

if [ -d "$HIPO_INPUT" ]; then

    HIPO_DIR="$HIPO_INPUT"

    while IFS= read -r HIPO_FILE
    do
        HIPO_FILES+=("$HIPO_FILE")
    done < <(find "$HIPO_DIR" -maxdepth 1 -type f -name "*.hipo" | sort)

elif [ -f "$HIPO_INPUT" ]; then

    HIPO_DIR=$(dirname "$HIPO_INPUT")

    case "$HIPO_INPUT" in
        *.hipo)
            HIPO_FILES+=("$HIPO_INPUT")
            ;;
        *)
            echo "Error: input file is not a .hipo file:"
            echo "  $HIPO_INPUT"
            exit 1
            ;;
    esac

else

    echo "Error: HIPO input does not exist:"
    echo "  $HIPO_INPUT"
    exit 1

fi

NFILES=${#HIPO_FILES[@]}

if [ "$NFILES" -eq 0 ]; then
    echo "No .hipo files found in:"
    echo "  $HIPO_INPUT"
    exit 1
fi

# ------------------------------------------------------------
# Dataset tag
# ------------------------------------------------------------
# Default:
#   folder input: basename of folder
#   file input:   basename of parent folder
# Example:
#   /.../dst/train/sidisdvcs       -> sidisdvcs
#   test_one_hipo/file.hipo        -> test_one_hipo
# For test folders, override with:
#   --dataset sidisdvcs
# ------------------------------------------------------------

if [ -z "$DATASET_TAG" ]; then
    DATASET_TAG=$(basename "$HIPO_DIR")
fi

# ------------------------------------------------------------
# Channel tag
# ------------------------------------------------------------

if [ "$MODE" = "sidis" ]; then

    if [ "$HADRON_PID" = "211" ]; then
        CHANNEL_TAG="epip"
    elif [ "$HADRON_PID" = "-211" ]; then
        CHANNEL_TAG="epim"
    elif [ "$HADRON_PID" = "321" ]; then
        CHANNEL_TAG="ekp"
    elif [ "$HADRON_PID" = "-321" ]; then
        CHANNEL_TAG="ekm"
    elif [ "$HADRON_PID" = "2212" ]; then
        CHANNEL_TAG="ep"
    else
        CHANNEL_TAG="epid${HADRON_PID}"
    fi

elif [ "$MODE" = "dihadron" ]; then

    CHANNEL_TAG="epippim"

else

    CHANNEL_TAG="inclusive_${REGION}"

fi

# ------------------------------------------------------------
# Output folders
# ------------------------------------------------------------

OPTION_TAG="eT${ELECTRON_TREE}pid${DET_PID_CUT}kin${KIN_CUT}${TARGET_FILTER_TAG}"
SUBDIR_TAG="skim_${DATASET_TAG}_${CHANNEL_TAG}_${OPTION_TAG}"

# Output folder can optionally include period and jobtag subfolders.

if [ -n "$DATA_PERIOD" ]; then
    ROOT_OUT_DIR="${ROOT_BASE_DIR}/${SUBDIR_TAG}/${DATA_PERIOD}"
    HIPO_OUT_DIR="${HIPO_BASE_DIR}/${SUBDIR_TAG}/${DATA_PERIOD}"
    LOG_OUT_DIR="${LOG_BASE_DIR}/${SUBDIR_TAG}/${DATA_PERIOD}"
else
    ROOT_OUT_DIR="${ROOT_BASE_DIR}/${SUBDIR_TAG}"
    HIPO_OUT_DIR="${HIPO_BASE_DIR}/${SUBDIR_TAG}"
    LOG_OUT_DIR="${LOG_BASE_DIR}/${SUBDIR_TAG}"
fi

if [ -n "$JOB_TAG" ]; then
    ROOT_OUT_DIR="${ROOT_OUT_DIR}/${JOB_TAG}"
    HIPO_OUT_DIR="${HIPO_OUT_DIR}/${JOB_TAG}"
    LOG_OUT_DIR="${LOG_OUT_DIR}/${JOB_TAG}"
fi

if [ -n "$OUT_PATH" ]; then
    ROOT_OUT_DIR="${OUT_PATH}/${ROOT_OUT_DIR}"
    HIPO_OUT_DIR="${OUT_PATH}/${HIPO_OUT_DIR}"
    LOG_OUT_DIR="${OUT_PATH}/${LOG_OUT_DIR}"
    ROOT_BASE_DIR="${OUT_PATH}/$ROOT_BASE_DIR"
    HIPO_BASE_DIR="${OUT_PATH}/$HIPO_BASE_DIR"
fi

mkdir -p "$ROOT_BASE_DIR"
mkdir -p "$LOG_BASE_DIR"
mkdir -p "$ROOT_OUT_DIR"
mkdir -p "$LOG_OUT_DIR"

if [ "$OUTPUT_FORMAT" = "hipo" ] || [ "$OUTPUT_FORMAT" = "both" ]; then
    mkdir -p "$HIPO_BASE_DIR"
    mkdir -p "$HIPO_OUT_DIR"
fi

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

echo "Executable:        $EXE"
echo "Mode:              $MODE"
echo "Hadron PID:        $HADRON_PID"
echo "Dataset tag:       $DATASET_TAG"
echo "Channel tag:       $CHANNEL_TAG"
echo "Region:            $REGION"
echo "Input HIPO:        $HIPO_INPUT"
echo "Resolved HIPO dir: $HIPO_DIR"
echo "Input HIPO files:  $NFILES"
echo "Electron tree:     $ELECTRON_TREE"
echo "Detector/PID cut:  $DET_PID_CUT"
echo "Kinematic cut:     $KIN_CUT"
echo "Output format:     $OUTPUT_FORMAT"
echo "Target map:        ${TARGET_MAP:-none}"
echo "Polarization src:  $POL_SOURCE"
echo "Target filter:     $TARGET_FILTER"
echo "Target filter tag: $TARGET_FILTER_TAG"
echo "Period mode:       $PERIOD_MODE"
echo "Period tag:        ${DATA_PERIOD:-none}"
echo "Run selection:     $RUN_SELECTION_MODE"

if [ "$RUN_SELECTION_MODE" = "test_one_run" ]; then
    echo "Test one run:      first file passing run/target filters"
fi

if [ "$RUN_SELECTION_MODE" = "single" ] || [ "$RUN_SELECTION_MODE" = "list" ] || [ "$RUN_SELECTION_MODE" = "file" ]; then
    echo "Selected runs:     ${SELECTED_RUNS[*]}"
fi

if [ "$RUN_SELECTION_MODE" = "range" ]; then
    echo "Run range:         ${RUN_RANGE_START}:${RUN_RANGE_END}"
fi

echo "Job tag:           ${JOB_TAG:-none}"
echo "Output path:           ${OUT_PATH:-none}"
echo "Target map runs:   $TARGET_MAP_RUNS"
echo "Target map skipped lines: $TARGET_MAP_SKIPPED_LINES"
echo "Option tag:        $OPTION_TAG"
echo "ROOT output dir:   $ROOT_OUT_DIR"
if [ "$OUTPUT_FORMAT" = "hipo" ] || [ "$OUTPUT_FORMAT" = "both" ]; then
    echo "HIPO output dir:   $HIPO_OUT_DIR"
fi
echo "Log output dir:    $LOG_OUT_DIR"
echo ""

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------

IFILE=0
N_PROCESSED=0
N_SKIPPED_TARGET=0
N_SKIPPED_RUNPARSE=0
N_SKIPPED_RUNSEL=0
N_FAILED=0
N_DONE=0
N_DONE_ROOT=0
N_DONE_HIPO=0

for HIPO_FILE in "${HIPO_FILES[@]}"
do
    IFILE=$((IFILE + 1))

    BASE=$(basename "$HIPO_FILE" .hipo)

    # Extract compact run tag from the HIPO file name.
    # Example:
    #   BASE=sidisdvcs_016270
    #   RUN_TAG=016270
    #   RUN_NUM=16270
    # If no 6 digit run number is found, RUN_TAG remains empty.
    RUN_TAG=$(echo "$BASE" | grep -oE '[0-9]{6}' | tail -1)

    RUN_NUM=""

    if [ -n "$RUN_TAG" ]; then
        RUN_NUM=$((10#$RUN_TAG))
    fi

    # ------------------------------------------------------------
    # Optional user run selection
    # ------------------------------------------------------------
    # This filters the HIPO file list before rgcskim is launched.
    # It does not change the skim code or physics cuts.
    # Filter order:
    #   1. user run selection
    #   2. target map matching, if requested
    #   3. rgcskim execution
    # ------------------------------------------------------------

    if [ "$RUN_SELECTION_MODE" != "all" ]; then

        if [ -z "$RUN_NUM" ]; then
            echo "[$IFILE / $NFILES] Skipping:"
            echo "  input:  $HIPO_FILE"
            echo "  reason: could not extract 6-digit run number from filename"
            echo ""

            N_SKIPPED_RUNPARSE=$((N_SKIPPED_RUNPARSE + 1))
            continue
        fi

        if ! run_passes_selection "$RUN_NUM"; then
            N_SKIPPED_RUNSEL=$((N_SKIPPED_RUNSEL + 1))
            continue
        fi

    fi

    # ------------------------------------------------------------
    # Optional target run filtering
    # ------------------------------------------------------------
    # If --targetfilter matched is selected, skip the HIPO file unless
    # the run number from the filename appears in the target CSV.
    # This is a file level prefilter to save skim time.
    # rgcskim still checks target matching internally using RUN::config.run.
    # ------------------------------------------------------------

    if [ "$TARGET_FILTER" = "matched" ]; then

        if [ -z "$RUN_NUM" ]; then
            echo "[$IFILE / $NFILES] Skipping:"
            echo "  input:  $HIPO_FILE"
            echo "  reason: could not extract 6-digit run number from filename"
            echo ""

            N_SKIPPED_RUNPARSE=$((N_SKIPPED_RUNPARSE + 1))
            continue
        fi

        if [ -z "${TARGET_RUNS[$RUN_NUM]+x}" ]; then
            echo "[$IFILE / $NFILES] Skipping:"
            echo "  input:  $HIPO_FILE"
            echo "  run:    $RUN_NUM"
            echo "  reason: run not found in target map"
            echo ""

            N_SKIPPED_TARGET=$((N_SKIPPED_TARGET + 1))
            continue
        fi
    fi

    if [ -z "$RUN_TAG" ]; then
        RUN_TAG="$BASE"
    fi

    FINAL_BASENAME="${DATASET_TAG}_${CHANNEL_TAG}_${OPTION_TAG}_${RUN_TAG}.root"
    FINAL_ROOT_FILE="${ROOT_OUT_DIR}/${FINAL_BASENAME}"

    FINAL_HIPO_BASENAME="${DATASET_TAG}_${CHANNEL_TAG}_${OPTION_TAG}_${RUN_TAG}.hipo"
    FINAL_HIPO_FILE="${HIPO_OUT_DIR}/${FINAL_HIPO_BASENAME}"

    LOGFILE="${LOG_OUT_DIR}/${DATASET_TAG}_${CHANNEL_TAG}_${OPTION_TAG}_${RUN_TAG}.log"

    # Temporary ROOT file names produced by rgcskim itself.
    # Current C++ output names:
    #   SIDIS:      rootfiles/skim_e<pip/pim/kp/km/p>_<BASE>.root
    #   dihadron:   rootfiles/skim_epippim_<BASE>.root
    #   inclusive:  rootfiles/skim_inclusive_<BASE>.root
    # The script moves/renames these to the final convention.

    if [ "$MODE" = "sidis" ]; then

        if [ "$HADRON_PID" = "211" ]; then
            CPP_HADRON_TAG="pip"
        elif [ "$HADRON_PID" = "-211" ]; then
            CPP_HADRON_TAG="pim"
        elif [ "$HADRON_PID" = "321" ]; then
            CPP_HADRON_TAG="kp"
        elif [ "$HADRON_PID" = "-321" ]; then
            CPP_HADRON_TAG="km"
        elif [ "$HADRON_PID" = "2212" ]; then
            CPP_HADRON_TAG="p"
        else
            CPP_HADRON_TAG="pid${HADRON_PID}"
        fi

        TEMP_ROOT_FILE="${ROOT_BASE_DIR}/skim_e${CPP_HADRON_TAG}_${BASE}.root"
        TEMP_HIPO_FILE="${HIPO_BASE_DIR}/skim_e${CPP_HADRON_TAG}_${BASE}.hipo"

    elif [ "$MODE" = "dihadron" ]; then

        TEMP_ROOT_FILE="${ROOT_BASE_DIR}/skim_epippim_${BASE}.root"
        TEMP_HIPO_FILE="${HIPO_BASE_DIR}/skim_epippim_${BASE}.hipo"

    else

        TEMP_ROOT_FILE="${ROOT_BASE_DIR}/skim_inclusive_${BASE}.root"
        TEMP_HIPO_FILE="${HIPO_BASE_DIR}/skim_inclusive_${BASE}.hipo"

    fi

    echo "[$IFILE / $NFILES] Processing:"
    echo "  input:      $HIPO_FILE"
    echo "  run:        ${RUN_NUM:-unknown}"
    if [ "$OUTPUT_FORMAT" = "root" ] || [ "$OUTPUT_FORMAT" = "both" ]; then
        echo "  temp ROOT:  $TEMP_ROOT_FILE"
        echo "  final ROOT: $FINAL_ROOT_FILE"
    fi
    if [ "$OUTPUT_FORMAT" = "hipo" ] || [ "$OUTPUT_FORMAT" = "both" ]; then
        echo "  temp HIPO:  $TEMP_HIPO_FILE"
        echo "  final HIPO: $FINAL_HIPO_FILE"
    fi
    echo "  log:        $LOGFILE"

    N_PROCESSED=$((N_PROCESSED + 1))

    # Remove a stale temporary output from a previous failed or interrupted run.
    rm -f "$TEMP_ROOT_FILE"
    rm -f "$TEMP_HIPO_FILE"

    if [ "$MODE" = "sidis" ]; then

        "$EXE" --mode sidis --pid "$HADRON_PID" "$HIPO_FILE" \
            --electrontree "$ELECTRON_TREE" \
            --detpidcut "$DET_PID_CUT" \
            --kincut "$KIN_CUT" \
            --outformat "$OUTPUT_FORMAT" \
            --outpath "$OUT_PATH" \
            "${TARGET_ARGS[@]}" \
            > "$LOGFILE" 2>&1

    elif [ "$MODE" = "dihadron" ]; then

        "$EXE" --mode dihadron "$HIPO_FILE" \
            --detpidcut "$DET_PID_CUT" \
            --kincut "$KIN_CUT" \
            --outformat "$OUTPUT_FORMAT" \
            --outpath "$OUT_PATH" \
            "${TARGET_ARGS[@]}" \
            > "$LOGFILE" 2>&1

    else

        "$EXE" --mode inclusive "$HIPO_FILE" \
            --region "$REGION" \
            --detpidcut "$DET_PID_CUT" \
            --kincut "$KIN_CUT" \
            --outformat "$OUTPUT_FORMAT" \
            --outpath "$OUT_PATH" \
            "${TARGET_ARGS[@]}" \
            > "$LOGFILE" 2>&1

    fi

    STATUS=$?

    if [ $STATUS -ne 0 ]; then
        echo "  Error: rgcskim failed with exit code $STATUS"
        echo "  Check log:"
        echo "    $LOGFILE"
        echo ""

        N_FAILED=$((N_FAILED + 1))

        if [ "$RUN_SELECTION_MODE" = "test_one_run" ]; then
            echo "Test-one-run mode: stopping after first file that passed file-level filters."
            echo ""
            break
        fi

        continue
    fi

    OUTPUT_OK=1

    if [ "$OUTPUT_FORMAT" = "root" ] || [ "$OUTPUT_FORMAT" = "both" ]; then
        if [ -f "$TEMP_ROOT_FILE" ]; then
            mv -f "$TEMP_ROOT_FILE" "$FINAL_ROOT_FILE"
            N_DONE_ROOT=$((N_DONE_ROOT + 1))
        else
            echo "  Warning: expected ROOT file was not found:"
            echo "    $TEMP_ROOT_FILE"
            echo "  Check log:"
            echo "    $LOGFILE"
            OUTPUT_OK=0
        fi
    fi

    if [ "$OUTPUT_FORMAT" = "hipo" ] || [ "$OUTPUT_FORMAT" = "both" ]; then
        if [ -f "$TEMP_HIPO_FILE" ]; then
            mv -f "$TEMP_HIPO_FILE" "$FINAL_HIPO_FILE"
            N_DONE_HIPO=$((N_DONE_HIPO + 1))
        else
            echo "  Warning: expected HIPO file was not found:"
            echo "    $TEMP_HIPO_FILE"
            echo "  Check log:"
            echo "    $LOGFILE"
            OUTPUT_OK=0
        fi
    fi

    if [ "$OUTPUT_OK" = "1" ]; then
        echo "  Done."
        N_DONE=$((N_DONE + 1))
    fi

    echo ""

    if [ "$RUN_SELECTION_MODE" = "test_one_run" ]; then
        echo "Test-one-run mode: stopping after first file that passed file-level filters."
        echo ""
        break
    fi

done

echo "All files processed."
echo ""
echo "Summary:"
echo "  Input HIPO files:           $NFILES"
echo "  Files processed:            $N_PROCESSED"
echo "  Files successfully written: $N_DONE"
echo "  ROOT files written:         $N_DONE_ROOT"
echo "  HIPO files written:         $N_DONE_HIPO"
echo "  Files failed:               $N_FAILED"
echo "  Files skipped, no run tag:  $N_SKIPPED_RUNPARSE"
echo "  Files skipped, run filter:  $N_SKIPPED_RUNSEL"
echo "  Files skipped, no target:   $N_SKIPPED_TARGET"
echo ""
if [ "$OUTPUT_FORMAT" = "root" ] || [ "$OUTPUT_FORMAT" = "both" ]; then
    echo "ROOT files are in:"
    echo "  $ROOT_OUT_DIR"
fi
if [ "$OUTPUT_FORMAT" = "hipo" ] || [ "$OUTPUT_FORMAT" = "both" ]; then
    echo "HIPO files are in:"
    echo "  $HIPO_OUT_DIR"
fi
echo "Logs are in:"
echo "  $LOG_OUT_DIR"
echo ""

JOB_END_EPOCH=$(date +%s)
JOB_END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
JOB_ELAPSED_SECONDS=$((JOB_END_EPOCH - JOB_START_EPOCH))
JOB_ELAPSED_TIME=$(format_elapsed_time "$JOB_ELAPSED_SECONDS")

echo "Timing:"
echo "  Start time: $JOB_START_TIME"
echo "  End time:   $JOB_END_TIME"
echo "  Elapsed:    $JOB_ELAPSED_TIME"
