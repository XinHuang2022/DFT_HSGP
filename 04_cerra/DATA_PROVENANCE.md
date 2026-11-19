# CERRA data provenance

The analysis data were extracted from the Copernicus European Regional
Reanalysis for the North Sea box used in the manuscript. The retained product
contains 10-m eastward wind at six-hour intervals from 1 December 2020 through
28 February 2021 at 3,701 marine native-grid locations.

Included analysis files:

- `u10_full_native_time_by_site.csv` -- 360-by-3,701 wind matrix with a
  leading time index;
- `site_manifest_full_native.csv` -- ordered native-grid coordinates and
  marine-site identifiers;
- `time_manifest_6hourly.csv` -- ordered UTC timestamps;
- `full_native_manifest.json` and `download_manifest.json` -- extraction
  metadata from the original workflow.

The raw reanalysis payloads are not included. Before a public release, the
current CERRA/Copernicus licence and attribution requirements should be
checked and stated in the repository licence documentation.
