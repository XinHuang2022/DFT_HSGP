function lambda = eq9_fourier_frequencies(T)
%EQ9_FOURIER_FREQUENCIES DFT-order frequencies in [-pi,pi).

    if T < 2 || fix(T) ~= T
        error('T must be an integer of at least two.');
    end
    lambda = 2 * pi * (0:T-1) / T;
    lambda(lambda >= pi) = lambda(lambda >= pi) - 2 * pi;
end



