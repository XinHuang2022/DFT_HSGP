function R=rtable(C,order)
%RTABLE Universal Hankel correlation table for q(x)=(1+x^3)^-2.
    if nargin<2, order=C.ref.order; end
    edges=0:C.ref.panelWidth:C.ref.xMax;
    if edges(end)<C.ref.xMax, edges=[edges,C.ref.xMax]; end %#ok<AGROW>
    x=[]; w=[];
    for p=1:numel(edges)-1
        [xp,wp]=gl(order,edges(p),edges(p+1));
        x=[x;xp]; w=[w;wp]; %#ok<AGROW>
    end
    z=(0:C.ref.dz:C.ref.zMax)';
    ker=w.*x./(1+x.^3).^2;
    val=zeros(size(z)); block=250;
    constant=9*sqrt(3)/(2*pi);
    for a=1:block:numel(z)
        ix=a:min(a+block-1,numel(z));
        val(ix)=constant*(ker.'*besselj(0,x*z(ix).')).';
    end
    rawZero=val(1); val=val/rawZero;
    R.z=z; R.r=val; R.xMax=C.ref.xMax; R.order=order;
    R.rawZero=rawZero; R.tailBound=constant/(4*C.ref.xMax^4);
end

