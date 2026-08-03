# RGC Workflow

This repository contains the analysis workflow for CLAS12 RGC ND3 studies.

The repository is organized into three main parts:

```text
skim/       Production of skimmed ROOT/HIPO files from raw HIPO files
analysis/   Downstream physics analysis using the skimmed files
tools/      General helper tools, including reusable JLab Slurm utilities
```

## Main workflow

```text
raw HIPO files
  -> skim/
  -> skimmed ROOT files
  -> analysis/
  -> plots and extracted observables
```

## Skim workflow

The skim code is located in:

```text
skim/
```

Typical interactive usage:

```bash
cd skim
source setup_env.sh
make
./scripts/run_skim.py configs/skim.yaml --dry-run
./scripts/run_skim.py configs/skim.yaml
```

For JLab Slurm batch submission:

```bash
cd skim
./scripts/batch/submit_slurm_skim.py configs/skim.yaml \
  --account hallc \
  --partition production \
  --array-limit 20 \
  --jobtag spring23_qe_kin1_full
```

More details are in:

```text
skim/README.md
skim/scripts/batch/README.md
```

## Analysis workflow

The downstream physics analysis code will be organized under:

```text
analysis/
```

This part uses the skimmed ROOT files produced by `skim/`.

## Generic Slurm tools

Reusable JLab HIPO Slurm utilities will be organized under:

```text
tools/jlab_slurm/
```

These tools are independent of the repo-specific `skim/` workflow.
