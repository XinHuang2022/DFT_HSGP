function result=run_era5_stage2_hsgp_refit_screen(cfg)
%RUN_ERA5_STAGE2_HSGP_REFIT_SCREEN Refit HSGP across a compact J ladder.
    if nargin<1,cfg=era5_northsea_config();end
    referenceFile=fullfile(cfg.paths.stage1Dir,'exact_reference.mat');
    auditFile=fullfile(cfg.paths.stage1Dir,'hsgp_J_audit_summary.csv');
    if ~isfile(referenceFile)||~isfile(auditFile)
        error('Complete Stage 1 before Stage 2.');
    end
    ensure_directory(cfg.paths.stage2Dir);ensure_directory(fullfile(cfg.paths.stage2Dir,'tasks'));
    ensure_directory(fullfile(cfg.paths.stage2Dir,'figures'));
    Jvalues=cfg.stage2.candidateJ(:);tasks=cell(numel(Jvalues),1);
    [parallelUsed,poolWorkers]=era5_open_pool(cfg,numel(Jvalues));
    if parallelUsed
        parfor q=1:numel(Jvalues),tasks{q}=era5_stage2_j_task(cfg,referenceFile,Jvalues(q));end
    else
        for q=1:numel(Jvalues),tasks{q}=era5_stage2_j_task(cfg,referenceFile,Jvalues(q));end
    end
    parts=cellfun(@(x)x.summary,tasks,'UniformOutput',false);summary=vertcat(parts{:});
    audit=readtable(auditFile);
    passMask=era5_logical_column(audit.passesAudit);
    passJ=audit.J(passMask);
    if isempty(passJ),recommendedJ=summary.J(end);else,recommendedJ=min(passJ);end
    writetable(summary,fullfile(cfg.paths.stage2Dir,'hsgp_refit_summary.csv'));
    runSummary=table(string(cfg.version),numel(Jvalues),recommendedJ, ...
        parallelUsed,poolWorkers,'VariableNames',{'version','nJValues', ...
        'recommendedJFromStage1','parallelUsed','poolWorkers'});
    writetable(runSummary,fullfile(cfg.paths.stage2Dir,'run_summary.csv'));
    make_era5_stage2_figures(summary,cfg);
    save(fullfile(cfg.paths.stage2Dir,'complete_result.mat'),'cfg','summary', ...
        'recommendedJ','runSummary','-v7.3');
    result=struct('summary',summary,'recommendedJ',recommendedJ,'runSummary',runSummary);
end
