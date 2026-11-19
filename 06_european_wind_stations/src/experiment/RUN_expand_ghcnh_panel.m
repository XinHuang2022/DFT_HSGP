clear; clc;
candidates=prepare_ghcnh_expanded_manifest();
fprintf('Expanded manifest: %d stations.\n',height(candidates));
summary=download_ghcnh_expanded_candidates();
disp(groupsummary(summary,'status'));
[coverage,~,~]=audit_ghcnh_expanded_completeness();
fprintf('Complete: %d; at most three gaps: %d.\n', ...
    sum(coverage.complete),sum(coverage.nearComplete));
