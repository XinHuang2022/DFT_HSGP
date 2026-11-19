function result = run_pacific_model_only(preset, useParallel)
%RUN_PACIFIC_MODEL_ONLY Run a model preset after an accepted audit.

    if nargin < 1 || isempty(preset)
        preset = 'smoke3';
    end
    if nargin < 2
        result = run_pacific_cv(preset);
    else
        result = run_pacific_cv(preset, useParallel);
    end
end
