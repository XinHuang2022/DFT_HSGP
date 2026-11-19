# CERRA code and result provenance

The 31-day source was copied from
`deliverables/cerra_northsea_eq9_suite_v0_2_0`; the 90-day source was copied
from `deliverables/cerra_northsea_eq9_suite_v0_2_1`. Both are byte-identical
to their trusted implementations under `wind_field_implement/Comp/CERRA`
before local path changes.

The preserved tables come from the completed stage-1 fixed-parameter screen,
stage-2 selective HSGP refits, and stage-3 comparison outputs for each record
length. They retain the deterministic 300-site evaluation set and nested
space-filling fitting subsets used by the trusted workflow.

Only path handling and release entry points were changed. The model,
preprocessing, split, mean design, parameter bounds, rank ladder, and
optimization definitions were not changed.
