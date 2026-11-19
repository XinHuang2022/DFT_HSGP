clear; clc;
C=ghcnh_trace_config();
summary=download_ghcnh_candidates(C.downloadCandidateCount); %#ok<NASGU>
[coverage,Y,targetTimes]=audit_ghcnh_completeness(C.downloadCandidateCount); %#ok<NASGU,ASGLU>
fprintf('Complete stations: %d\n',sum(coverage.complete));
fprintf('Near-complete stations: %d\n',sum(coverage.nearComplete));
