# European station code and result provenance

Experiment-specific source and retained results were copied from
`deliverables/ghcnh_trace_hsgp_v0_1_0`. The shared Equation-9 HSGP functions
were copied from `deliverables/era5_northsea_eq9_suite_v0_1_0`, and the
adaptive-rank helpers were copied from
`deliverables/realdata_adaptive_hsgp_v0_1_0`.

The release changes localize data and output paths, standardize the clean
five-fold result layout, and add a five-worker outer-fold driver. They do not
alter the station panel, mean design, covariance model, temporal DFT, Whittle
objective, fold construction, optimizer settings, rank ladders, or stopping
rules used by the retained experiment.
