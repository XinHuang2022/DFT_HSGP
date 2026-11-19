function solution = eq9_lowrank_solve(factor, rightHandSide)
%EQ9_LOWRANK_SOLVE Apply the stable QR-based low-rank inverse.

    projected = factor.Q' * rightHandSide;
    orthogonal = rightHandSide - factor.Q * projected;
    reduced = factor.L' \ (factor.L \ projected);
    solution = orthogonal / factor.noiseVariance + factor.Q * reduced;
end
