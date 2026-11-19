function B=mkbasis(xy,c,M)
%MKBASIS Dirichlet basis, radially truncated by Laplacian eigenvalue.
    lo=min(xy,[],1); hi=max(xy,[],1); ctr=(lo+hi)/2;
    half=(hi-lo)/2; L=c*half;
    if any(L<=0), error('Spatial design must span both coordinates.'); end
    J=ceil(sqrt(M)*1.7)+3;
    [j1,j2]=ndgrid(1:J,1:J);
    ind=[j1(:),j2(:)];
    om=pi*ind./(2*L);
    eigv=sum(om.^2,2);
    [~,ix]=sortrows([eigv,ind(:,1),ind(:,2)],[1,2,3]);
    ix=ix(1:M); ind=ind(ix,:); om=om(ix,:); eigv=eigv(ix);
    B.c=c; B.M=M; B.center=ctr; B.L=L; B.ind=ind;
    B.omega=om; B.radial=sqrt(eigv);
    B.Phi=features(xy,B);
end

