clear; clc;
C=ghcnh_trace_config();
summary=download_ghcnh_candidates(C.pilotDownloadCount); %#ok<NASGU>
[coverage,Y,targetTimes]=audit_ghcnh_completeness(C.pilotDownloadCount); %#ok<NASGU,ASGLU>
disp(coverage(:,{'station','country','name','summerRecordCount','nMatched','nMissing','complete'}));
