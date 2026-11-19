function label = cerra_hsgp_label(boundaryFactor, J)
%CERRA_HSGP_LABEL Scalar method label carrying c and J.
    label = string(sprintf('hsgp_c%s_J%02d', ...
        cerra_boundary_tag(boundaryFactor), J));
end

