# Direct-multiplier source provenance

The release files were copied without mathematical changes from:

- matched frequency-adaptive experiment:
  `deliverables/spec_hsgp_adaptive_matched_v0_1_0`;
- direct-multiplier model and numerical-inversion reference:
  `deliverables/spec_hsgp_v1/src`;
- retained matched run:
  `deliverables/spec_hsgp_adaptive_matched_v0_1_0/results/matched_20_20261107_230014`.

Only path construction was changed so the code and outputs are local to this
folder.  The sectioned entry point reports the three methods retained in the
current manuscript.  The underlying validated driver still computes the
initial screened fit because that fit initializes the adaptive workflow.

