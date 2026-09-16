#!/bin/bash -l
#SBATCH --job-name=rgcskim_array
#SBATCH --partition=production
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Generic Slurm array worker for rgcskim.
# This file is normally submitted by:
#   scripts/batch/submit_slurm_skim.py
#
# One Slurm array task processes one HIPO file.

set -u

require_env(){
    local name="$1"
    local value="${!name:-}"

    if [ -z "$value" ]; then
        echo "Error: required environment variable is missing: $name"
        exit 2
    fi
}

require_env RGC_SKIM_REPO
require_env RGC_SKIM_FILELIST
require_env RGC_SKIM_RUNNER
require_env RGC_SKIM_MODE
require_env RGC_SKIM_DETPIDCUT
require_env RGC_SKIM_KINCUT
require_env RGC_SKIM_PERIOD
require_env RGC_SKIM_OUTFORMAT

cd "$RGC_SKIM_REPO" || exit 2

echo "============================================================"
echo "RGC skim Slurm array task"
echo "============================================================"
echo "Job ID:        ${SLURM_JOB_ID:-unknown}"
echo "Array job ID:  ${SLURM_ARRAY_JOB_ID:-unknown}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID:-unknown}"
echo "Host:          $(hostname)"
echo "PWD:           $(pwd)"
echo "File list:     $RGC_SKIM_FILELIST"
echo ""

echo "Loading CLAS12 environment and software for analysis"
echo ""

module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module purge
module load clas12

# The rgcskim Makefile and executable expect this HIPO layout.
export HIPO=/group/clas12/packages/hipo/dev

echo ""
echo "Environment summary"
echo "------------------------------------------------------------"
echo "ROOTSYS: ${ROOTSYS:-}"
echo "HIPO:    ${HIPO:-}"
echo ""

if [ -z "${ROOTSYS:-}" ]; then
    echo "Error: ROOTSYS is not set after module load."
    exit 2
fi

if [ ! -x "./rgcskim" ]; then
    echo "Error: ./rgcskim is missing or not executable."
    echo "Compile first on ifarm:"
    echo "  source setup_env.sh"
    echo "  make"
    exit 2
fi

if [ ! -f "$RGC_SKIM_FILELIST" ]; then
    echo "Error: file list does not exist:"
    echo "  $RGC_SKIM_FILELIST"
    exit 2
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
LINE_NO=$((TASK_ID + 1))

HIPO_FILE=$(sed -n "${LINE_NO}p" "$RGC_SKIM_FILELIST")

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

echo "Input selected for this task"
echo "------------------------------------------------------------"
echo "Array task: $TASK_ID"
echo "Line no:    $LINE_NO"
echo "HIPO file:  $HIPO_FILE"
echo ""

CMD=(bash "$RGC_SKIM_RUNNER")

if [ -n "${RGC_SKIM_DATASET:-}" ]; then
    CMD+=(--dataset "$RGC_SKIM_DATASET")
fi

CMD+=(--mode "$RGC_SKIM_MODE")

if [ "$RGC_SKIM_MODE" = "sidis" ]; then
    if [ -z "${RGC_SKIM_PID:-}" ]; then
        echo "Error: SIDIS mode requires RGC_SKIM_PID."
        exit 2
    fi
    CMD+=(--pid "$RGC_SKIM_PID")
fi

if [ "$RGC_SKIM_MODE" = "inclusive" ]; then
    CMD+=(--region "${RGC_SKIM_REGION:-all}")
fi

CMD+=("$HIPO_FILE")

if [ "$RGC_SKIM_MODE" = "sidis" ]; then
    CMD+=(--electrontree "${RGC_SKIM_ELECTRONTREE:-0}")
fi

CMD+=(--detpidcut "$RGC_SKIM_DETPIDCUT")
CMD+=(--kincut "$RGC_SKIM_KINCUT")

if [ -n "${RGC_SKIM_TARGETMAP:-}" ]; then
    CMD+=(--targetmap "$RGC_SKIM_TARGETMAP")
    CMD+=(--polsource "${RGC_SKIM_POLSOURCE:-offline}")
    CMD+=(--targetfilter "${RGC_SKIM_TARGETFILTER:-all}")
fi

CMD+=(--period "$RGC_SKIM_PERIOD")
CMD+=(--outformat "$RGC_SKIM_OUTFORMAT")

if [ -n "${RGC_SKIM_JOBTAG:-}" ]; then
    CMD+=(--jobtag "$RGC_SKIM_JOBTAG")
fi

if [ -n "${RGC_SKIM_OUTPATH:-}" ]; then
    CMD+=(--out-path "$RGC_SKIM_OUTPATH")
fi

echo "Command"
echo "------------------------------------------------------------"
printf '%q ' "${CMD[@]}"
echo ""
echo ""

START_TIME=$(date '+%Y-%m-%d %H:%M:%S')
START_EPOCH=$(date +%s)

echo "Start time: $START_TIME"
echo ""

"${CMD[@]}"
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

exit $STATUS
