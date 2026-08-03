#!/bin/bash -l
#SBATCH --job-name=generic_hipo
#SBATCH --partition=production
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Generic Slurm array worker for running a user command over HIPO files.
# Normally submitted by:
#   tools/jlab_slurm/submit_generic_hipo_array.py
#
# One Slurm array task processes one HIPO file.

set -u

if [ -z "${GENERIC_JOB_ENV:-}" ]; then
    echo "Error: GENERIC_JOB_ENV is not set."
    exit 2
fi

if [ ! -f "$GENERIC_JOB_ENV" ]; then
    echo "Error: GENERIC_JOB_ENV file does not exist:"
    echo "  $GENERIC_JOB_ENV"
    exit 2
fi

source "$GENERIC_JOB_ENV"

require_env(){
    local name="$1"
    local value="${!name:-}"

    if [ -z "$value" ]; then
        echo "Error: required environment variable is missing: $name"
        exit 2
    fi
}

require_env GENERIC_FILELIST
require_env GENERIC_WORK_DIR
require_env GENERIC_OUTPUT_DIR
require_env GENERIC_COMMAND_TEMPLATE
require_env GENERIC_LOAD_CLAS12_ENV
require_env GENERIC_HIPO_PATH

cd "$GENERIC_WORK_DIR" || exit 2

mkdir -p "$GENERIC_OUTPUT_DIR"

echo "============================================================"
echo "Generic JLab HIPO Slurm array task"
echo "============================================================"
echo "Job ID:        ${SLURM_JOB_ID:-unknown}"
echo "Array job ID:  ${SLURM_ARRAY_JOB_ID:-unknown}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID:-unknown}"
echo "Host:          $(hostname)"
echo "PWD:           $(pwd)"
echo "File list:     $GENERIC_FILELIST"
echo "Output dir:    $GENERIC_OUTPUT_DIR"
echo ""

if [ "$GENERIC_LOAD_CLAS12_ENV" = "1" ]; then
    echo "Loading CLAS12 environment and software for analysis"
    echo ""

    module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
    module purge
    module load clas12

    export HIPO="$GENERIC_HIPO_PATH"

    echo ""
    echo "Environment summary"
    echo "------------------------------------------------------------"
    echo "ROOTSYS: ${ROOTSYS:-}"
    echo "HIPO:    ${HIPO:-}"
    echo ""
else
    echo "Skipping CLAS12 module load because GENERIC_LOAD_CLAS12_ENV=0"
    echo ""
fi

if [ ! -f "$GENERIC_FILELIST" ]; then
    echo "Error: file list does not exist:"
    echo "  $GENERIC_FILELIST"
    exit 2
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
LINE_NO=$((TASK_ID + 1))

HIPO_FILE=$(sed -n "${LINE_NO}p" "$GENERIC_FILELIST")

if [ -z "$HIPO_FILE" ]; then
    echo "Error: no HIPO file found for array task $TASK_ID"
    echo "Line number requested: $LINE_NO"
    exit 2
fi

if [ ! -f "$HIPO_FILE" ]; then
    echo "Error: HIPO file does not exist:"
    echo "  $HIPO_FILE"
    exit 2
fi

HIPO_BASENAME=$(basename "$HIPO_FILE")
HIPO_STEM="${HIPO_BASENAME%.hipo}"

RUN_NUMBER=$(python3 - "$HIPO_FILE" <<'PY'
import re
import sys

matches = re.findall(r"[0-9]{6}", sys.argv[1])
print(str(int(matches[-1])) if matches else "unknown")
PY
)

echo "Input selected for this task"
echo "------------------------------------------------------------"
echo "Array task: $TASK_ID"
echo "Line no:    $LINE_NO"
echo "Run:        $RUN_NUMBER"
echo "HIPO file:  $HIPO_FILE"
echo ""

COMMAND=$(python3 - \
    "$GENERIC_COMMAND_TEMPLATE" \
    "$HIPO_FILE" \
    "$HIPO_BASENAME" \
    "$HIPO_STEM" \
    "$RUN_NUMBER" \
    "$TASK_ID" \
    "${SLURM_JOB_ID:-unknown}" \
    "${SLURM_ARRAY_JOB_ID:-unknown}" \
    "$GENERIC_OUTPUT_DIR" \
    "$GENERIC_WORK_DIR" <<'PY'
import shlex
import sys

template = sys.argv[1]

mapping = {
    "hipo_file": sys.argv[2],
    "hipo_file_q": shlex.quote(sys.argv[2]),
    "hipo_basename": sys.argv[3],
    "hipo_stem": sys.argv[4],
    "run": sys.argv[5],
    "task_id": sys.argv[6],
    "job_id": sys.argv[7],
    "array_job_id": sys.argv[8],
    "output_dir": sys.argv[9],
    "output_dir_q": shlex.quote(sys.argv[9]),
    "work_dir": sys.argv[10],
    "work_dir_q": shlex.quote(sys.argv[10]),
}

try:
    print(template.format(**mapping))
except KeyError as exc:
    missing = exc.args[0]
    print(f"Error: unknown placeholder in command template: {{{missing}}}", file=sys.stderr)
    print("Allowed placeholders:", file=sys.stderr)
    for key in sorted(mapping):
        print(f"  {{{key}}}", file=sys.stderr)
    sys.exit(2)
except ValueError as exc:
    print("Error while formatting command template:", file=sys.stderr)
    print(str(exc), file=sys.stderr)
    print("If your command needs literal curly braces, escape them as {{ and }}.", file=sys.stderr)
    sys.exit(2)
PY
)

STATUS=$?

if [ "$STATUS" -ne 0 ]; then
    echo "Error: failed to build command from template."
    exit "$STATUS"
fi

echo "Command"
echo "------------------------------------------------------------"
echo "$COMMAND"
echo ""

START_TIME=$(date '+%Y-%m-%d %H:%M:%S')
START_EPOCH=$(date +%s)

echo "Start time: $START_TIME"
echo ""

bash -lc "$COMMAND"
STATUS=$?

END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
END_EPOCH=$(date +%s)
ELAPSED=$((END_EPOCH - START_EPOCH))

echo ""
echo "Task summary"
echo "------------------------------------------------------------"
echo "End time:     $END_TIME"
echo "Elapsed sec:  $ELAPSED"
echo "Exit status:  $STATUS"

exit "$STATUS"
