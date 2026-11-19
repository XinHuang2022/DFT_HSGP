# European station data provenance

The analysis panel was extracted from NOAA's Global Historical Climatology
Network hourly station archive (GHCNh). The retained matrix contains 360
six-hourly wind-speed observations from June through August 2024 at 839
European stations with complete matched records.

Included analysis files:

- `ghcnh_complete_panel_expanded.mat` -- MATLAB analysis structure;
- `ghcnh_complete_station_manifest_expanded.csv` -- ordered station metadata
  and extraction diagnostics;
- `cv_fold_manifest.csv` -- protected deterministic whole-site fold assignment.

Raw station downloads are not duplicated because the extracted panel is the
only input required by the MATLAB analysis. NOAA/GHCNh attribution and current
redistribution terms should be confirmed before public release.
