# Target maps

This folder contains target and polarization CSV files used by the skim.

The shared baseline target map is tracked in Git:

polarization_vs_run_results_preliminary_v1__1_(polarization_vs_run_results_pre).csv

Other local CSV files are ignored by Git unless they are explicitly added to .gitignore.

Expected CSV columns are:

run,start_time,stop_time,species,cell,charge_avg_online,charge_avg_offline,charge_avg_offline_err,run_dose(Pe/cm2)

Example use from configs/skim.yaml:

target:
  map: configs/target_maps/polarization_vs_run_results_preliminary_v1__1_(polarization_vs_run_results_pre).csv
  polsource: offline
  targetfilter: matched