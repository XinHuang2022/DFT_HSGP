function x=lrsolve(F,b)
%LRSOLVE Apply a cached low-rank or dense covariance inverse.
    if strcmp(F.kind,'woodbury')
        u=F.P'*b;
        x=b/F.vn-F.P*(F.L'\(F.L\u))/(F.vn^2);
    else
        x=F.L'\(F.L\b);
    end
end

