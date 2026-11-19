function r=rcorr(z,R)
%RCORR Interpolate the independently tabulated spatial correlation.
    if any(z(:)<0) || max(z(:))>R.z(end)
        error('Requested z is outside the verified Hankel table.');
    end
    r=interp1(R.z,R.r,z,'pchip');
end

