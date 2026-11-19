function S=sim_ref(theta,xy,T,R,seed)
%SIM_REF Exact circular latent field plus independent observation noise.
    rng(seed,'twister'); n=size(xy,1); lam=freqs(T);
    D=complex(zeros(T,n)); jit=zeros(T,1);
    [L,jit(1)]=cholj(cov_ref(theta,xy,xy,lam(1),R,false));
    D(1,:)=(L*randn(n,1))';
    if mod(T,2)==0, last=T/2; ny=T/2+1; else, last=(T+1)/2; ny=[]; end
    for k=2:last
        [L,jit(k)]=cholj(cov_ref(theta,xy,xy,lam(k),R,false));
        d=L*(randn(n,1)+1i*randn(n,1))/sqrt(2);
        partner=T-k+2; D(k,:)=d.'; D(partner,:)=conj(d).';
        jit(partner)=jit(k);
    end
    if ~isempty(ny)
        [L,jit(ny)]=cholj(cov_ref(theta,xy,xy,lam(ny),R,false));
        D(ny,:)=(L*randn(n,1))';
    end
    Z=ifft(D,[],1)*sqrt(T); X=real(Z);
    noise=sqrt(theta(5))*randn(T,n);
    S.X=X; S.noise=noise; S.Y=X+noise; S.D=D;
    S.maxImag=max(abs(imag(Z)),[],'all');
    S.maxJitter=max(jit); S.seed=seed;
end
