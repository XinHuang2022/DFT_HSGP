function lambda=freqs(T)
%FREQS MATLAB-order angular Fourier frequencies.
    k=0:T-1;
    k(k>floor(T/2))=k(k>floor(T/2))-T;
    lambda=2*pi*k/T;
end

