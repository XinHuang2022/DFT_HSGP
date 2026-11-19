# ERA5 North Sea data

This directory contains the fixed all-ocean 17-by-17 ERA5 North Sea panel:
54--58 degrees N, 2--6 degrees E, 0.25-degree spacing, 2021-12-01 through
2022-02-28 at 00/06/12/18 UTC (360 times).  `u10_time_by_site.csv` is the
primary response matrix used by suite v0.1.0.  The NetCDF and `v10` export are
retained so later experiments can study vector wind without redownloading.

See `download_manifest.json` for the source and retrieval metadata.  The
reproducible downloader is stored at the parent `ERA5` level after delivery.
