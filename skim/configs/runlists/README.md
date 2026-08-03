# Run lists

Put optional run list text files here.

Format:

16270
16271
16275

Blank lines and lines starting with # are ignored by scripts/run_all_hipo.sh.

Example use from configs/skim.yaml:

run_selection:
  mode: file
  runlist: configs/runlists/test_runs.txt