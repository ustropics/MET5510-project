function XVx = XV2XVx(XV)
% Compute zonal derivative XVx
% From Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf
% Page 9 (3D meridional wind field)
global jj kk ll m0 Lx cplx

XVx = zeros(ll, 1);

for k = 1:kk+1
    for j = 2:jj
        l = jk2l(j,k);
        XVx(l) = cplx*(2*pi*m0/Lx) * XV(l);
    end
end
end