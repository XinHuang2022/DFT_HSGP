function design = make_ghcnh_mean_design(T,xyStandardized)
%MAKE_GHCNH_MEAN_DESIGN Spatial plane plus common diurnal harmonics.

    n=size(xyStandardized,1);
    time=(0:T-1)';
    temporal=[sin(2*pi*time/4),cos(2*pi*time/4),cos(pi*time)];
    p=6; B=zeros(T,n,p);
    B(:,:,1)=1;
    B(:,:,2)=repmat(xyStandardized(:,1)',T,1);
    B(:,:,3)=repmat(xyStandardized(:,2)',T,1);
    for j=1:3, B(:,:,3+j)=repmat(temporal(:,j),1,n); end
    design.model="spatial_plane_plus_diurnal_harmonics";
    design.B=B;
    design.Afrequency=fft(B,[],1)/sqrt(T);
    design.coefficientNames=["intercept";"east";"north"; ...
        "diurnal_sine";"diurnal_cosine";"semidiurnal_cosine"];
    design.nCoefficients=p; design.nStations=n; design.T=T;
end
