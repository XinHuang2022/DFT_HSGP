function tag = cerra_boundary_tag(boundaryFactor)
%CERRA_BOUNDARY_TAG File-safe deterministic HSGP boundary-factor tag.
    text = regexprep(sprintf('%.3f', boundaryFactor), '0+$', '');
    text = regexprep(text, '\.$', '');
    tag = strrep(text, '.', 'p');
end

