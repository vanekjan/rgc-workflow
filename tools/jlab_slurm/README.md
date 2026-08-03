# Generic JLab HIPO Slurm Tools

This folder contains generic Slurm submission helpers for running any command over HIPO files on the JLab farm.

These tools are independent of the repo-specific `skim/` workflow.

## Files

```text
submit_generic_hipo_array.py        User-facing submission script
generic_hipo_array_worker.sbatch    Slurm array worker script
```

Normally, users should run `submit_generic_hipo_array.py`, not the `.sbatch` worker directly.

## Basic idea

```text
one Slurm array task = one HIPO file
```

The submitter creates a HIPO file list, writes a small job environment file, and submits one Slurm array job. Each array task reads one HIPO file from the file list and runs the user-provided command template.

## Minimal example

From the repository root:

```bash
./tools/jlab_slurm/submit_generic_hipo_array.py \
  --input /path/to/hipo/folder \
  --account hallc \
  --partition priority \
  --max-files 1 \
  --jobtag generic_test \
  --command 'echo "task={task_id} run={run} file={hipo_file}"' \
  --dry-run
```

Submit the same test without `--dry-run`:

```bash
./tools/jlab_slurm/submit_generic_hipo_array.py \
  --input /path/to/hipo/folder \
  --account hallc \
  --partition priority \
  --max-files 1 \
  --jobtag generic_test \
  --command 'echo "task={task_id} run={run} file={hipo_file}"'
```

## Command template placeholders

The command supplied with `--command` may use these placeholders:

```text
{hipo_file}        Full HIPO file path
{hipo_file_q}      Shell-quoted HIPO file path
{hipo_basename}    File name only
{hipo_stem}        File name without .hipo
{run}              Run number extracted from the file name, or unknown
{task_id}          Slurm array task ID
{job_id}           Slurm job ID
{array_job_id}     Slurm array job ID
{output_dir}       Output directory
{output_dir_q}     Shell-quoted output directory
{work_dir}         Working directory
{work_dir_q}       Shell-quoted working directory
```

If your command contains literal curly braces, escape them as `{{` and `}}`.

## Example: run a custom executable

```bash
./tools/jlab_slurm/submit_generic_hipo_array.py \
  --input /lustre24/expphy/cache/clas12/rg-c/production/spring23/pass1/ND3/dst/train/sidisdvcs/ \
  --account hallc \
  --partition production \
  --array-limit 20 \
  --jobtag my_analysis \
  --command './my_analysis_code {hipo_file_q} --run {run} --out {output_dir_q}/analysis_{run}.root'
```

## Example: run a ROOT macro

```bash
./tools/jlab_slurm/submit_generic_hipo_array.py \
  --input /path/to/hipo/folder \
  --account hallc \
  --partition production \
  --array-limit 20 \
  --jobtag root_macro_test \
  --command 'root -l -b -q "my_macro.C(\"{hipo_file}\",\"{output_dir}/out_{run}.root\")"'
```

## Run filtering

The submitter can filter by run number before submitting.

Single run:

```bash
--run 17491
```

Run list:

```bash
--runs 17491,17492,17493
```

Run range:

```bash
--runrange 17491:17520
```

Runlist file:

```bash
--runlist configs/runlists/my_runs.txt
```

These filters can be combined. A file must pass all requested filters.

## Useful options

```text
--input          HIPO input file, HIPO directory, or text file list
--command        Command template to run for each HIPO file
--account        Slurm account, for example hallc
--partition      Slurm partition, usually priority for tests and production for full jobs
--max-files      Limit number of files for testing
--array-limit    Maximum number of simultaneous array tasks
--jobtag         Name tag for logs and output folders
--work-dir       Directory where the command should run
--output-dir     Directory available as {output_dir}
--state-dir      Directory for generated file lists, env files, and Slurm logs
--dry-run        Print the resolved submission without submitting
```

## Output locations

By default, generated helper files and Slurm logs go under:

```text
logs/generic_slurm/
```

Default command output directory:

```text
generic_slurm_output/<jobtag>_<timestamp>/
```

## CLAS12 environment

By default, the Slurm worker loads the CLAS12 module environment and sets:

```bash
export HIPO=/group/clas12/packages/hipo/dev
```

Use `--no-clas12-env` only if your command does not need the CLAS12 software environment.
