# RGC skim framework

This repository contains the local RGC HIPO skim framework for inclusive, SIDIS, and dihadron studies.

## Main features

- Inclusive electron skim with region selection: all, dis, res, qe
- Single hadron SIDIS skim
- Dihadron e pi+ pi- skim
- Optional loose detector and PID cuts
- Target map and polarization metadata branches
- Target matched run filtering
- Run selection by single run, run list, run range, or one test run
- ROOT output, filtered HIPO output, or both
- HIPO output preserves selected physics events and HEL::scaler events for FCup and helicity normalization
- Output folders can be separated by run period such as summer22, fall22, and spring23

## Directory layout

rgcskim/
  Makefile
  README.md
  setup_env.sh
  src/
    rgcskim.cc
  scripts/
    run_all_hipo.sh
    run_skim.py
  configs/
    skim.example.yaml
    runlists/
    target_maps/

Produced files are ignored by Git:

rootfiles/
hipofiles/
logs/
build/
rgcskim

## Setup

From the project folder:

source setup_env.sh
make clean
make

## YAML workflow

Start from the example file:

cp configs/skim.example.yaml configs/skim.yaml

Edit:

vi configs/skim.yaml

Always dry run first:

./scripts/run_skim.py configs/skim.yaml --dry-run

Then run:

./scripts/run_skim.py configs/skim.yaml

## Direct script workflow

Example inclusive QE skim:

./scripts/run_all_hipo.sh 
  --dataset sidisdvcs \
  --mode inclusive \
  --region qe \
  /path/to/sidisdvcs/ \
  --detpidcut 1 \
  --targetmap configs/target_maps/your_target_map.csv \
  --polsource offline \
  --targetfilter matched \
  --period auto \
  --outformat root \
  --test-one-run \
  --jobtag test_one_run

## Output formats

root  means write ROOT skim only
hipo  means write filtered HIPO skim only
both  means write ROOT and filtered HIPO skim

## Useful checks

Compile:

make clean
make

Check scripts:

bash -n scripts/run_all_hipo.sh
python3 -m py_compile scripts/run_skim.py

Dry run:

./scripts/run_skim.py configs/skim.yaml --dry-run

## Normal Git workflow

Check changes:

git status
git diff

Stage files:

git add file_name

Commit:

git commit -m "Short description of the change"

Push:

git push
