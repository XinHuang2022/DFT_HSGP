function mask = era5_logical_column(value)
%ERA5_LOGICAL_COLUMN Convert CSV truth columns to a logical row mask.
%
% MATLAB may infer a column written from logical values as numeric 0/1,
% logical, string, categorical, or cell text depending on release and file
% contents.  Numeric 0/1 values must not be used directly as table indices.

    if islogical(value)
        mask = value;
        return;
    end
    if isnumeric(value)
        if any(~isfinite(value)) || any(~ismember(value, [0, 1]))
            error('Numeric truth column contains values other than 0 and 1.');
        end
        mask = value ~= 0;
        return;
    end
    textValue = lower(strtrim(string(value)));
    trueText = ismember(textValue, ["true", "1", "yes"]);
    falseText = ismember(textValue, ["false", "0", "no"]);
    if any(~trueText & ~falseText)
        error('Text truth column contains an unrecognized value.');
    end
    mask = trueText;
end
