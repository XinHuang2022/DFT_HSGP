clear; clc;
suiteRoot=fileparts(mfilename('fullpath')); addpath(suiteRoot);
O=struct();
O.referenceJ=[64 72 80 88 96];
O.candidateJ=[16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92];
O.epsilonReferencePerSite=0.03;
O.epsilonKlPerSite=0.05;
O.auditOnly=false;
O.outputRoot=fullfile(suiteRoot,'Results_extended_relaxed','era5');
run_adaptive_realdata('era5',O);

