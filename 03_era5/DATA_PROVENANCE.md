# ERA5 data provenance

The analysis data were extracted from the ERA5 reanalysis for the North Sea
box used in the manuscript. The retained product contains 10-m eastward wind
at six-hour intervals for 1 December 2020 through 28 February 2021 on a
17-by-17 all-ocean grid.

Included analysis files:

- `u10_time_by_site.csv` -- 360-by-289 wind matrix with a leading time index;
- `site_manifest_17x17.csv` -- ordered coordinates and grid identifiers;
- `time_manifest_6hourly.csv` -- ordered UTC timestamps;
- `download_manifest.json` -- extraction metadata from the original workflow.

The raw NetCDF download is not included because the extracted matrix is the
only data product required by the MATLAB analysis. Before a public release,
the current ERA5/Copernicus licence and attribution requirements should be
checked and stated in the repository licence documentation.
