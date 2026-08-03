# RGC Skim Slurm Submission

This folder contains Slurm helpers for submitting `rgcskim` jobs to the JLab farm.

These scripts are specific to the `skim/` workflow in this repository.

## Files

```text
submit_slurm_skim.py   User-facing submission script
rgcskim_array.sbatch   Slurm array worker script
```

Normally, users should run `submit_slurm_skim.py`, not `rgcskim_array.sbatch` directly.

## Basic idea

```text
one Slurm array task = one HIPO file
```

For example, if the YAML config resolves to 100 HIPO files, the submitter creates one Slurm array job with 100 tasks. Each task processes one HIPO file using `scripts/run_all_hipo.sh`.

## Required setup

From the `skim/` directory:

```bash
source setup_env.sh
make
```

Make sure `./rgcskim` exists before submitting jobs.

## Dry run

Always test with `--dry-run` first:

```bash
./scripts/batch/submit_slurm_skim.py configs/skim.yaml \
  --account hallc \
  --partition priority \
  --max-files 1 \
  --jobtag slurm_test_one \
  --dry-run
```

This prints the resolved file list, selected mode, cuts, output paths, and final `sbatch` command without submitting.

## Small test job

Submit one test file:

```bash
./scripts/batch/submit_slurm_skim.py configs/skim.yaml \
  --account hallc \
  --partition priority \
  --max-files 1 \
  --jobtag slurm_test_one
```

Submit three test files:

```bash
./scripts/batch/submit_slurm_skim.py configs/skim.yaml \
  --account hallc \
  --partition priority \
  --max-files 3 \
  --jobtag slurm_test_three
```

Check running jobs:

```bash
squeue -u $USER
```

Check Slurm logs:

```bash
ls -ltr logs/slurm
cat logs/slurm/*slurm_test_one*.out
cat logs/slurm/*slurm_test_one*.err
```

Check output ROOT files:

```bash
find rootfiles -name '*slurm_test_one*.root'
```

## Full production example

For a full production job, use the `production` partition:

```bash
./scripts/batch/submit_slurm_skim.py configs/skim.yaml \
  --account hallc \
  --partition production \
  --array-limit 20 \
  --jobtag spring23_qe_kin1_full
```

The `--array-limit 20` option means at most 20 array tasks run at the same time.

## Important options

```text
--account        Slurm account, for example hallc
--partition      Slurm partition, usually priority for tests and production for full jobs
--max-files      Limit number of HIPO files for testing
--array-limit    Maximum number of simultaneous array tasks
--jobtag         Extra output subfolder name
--dry-run        Print the resolved submission without submitting
```

## Output locations

Slurm stdout/stderr logs:

```text
logs/slurm/
```

Generated file lists:

```text
logs/slurm_filelists/
```

Normal `rgcskim` output logs:

```text
logs/skim_<channel>_<options>/
```

Normal `rgcskim` ROOT files:

```text
rootfiles/skim_<channel>_<options>/
```

## Notes

The Slurm worker uses a bash login shell and loads the CLAS12 module environment directly. It sets:

```bash
export HIPO=/group/clas12/packages/hipo/dev
```

This is intentional because the interactive `setup_env.sh` is csh-style and should not be sourced directly inside the bash Slurm worker.
