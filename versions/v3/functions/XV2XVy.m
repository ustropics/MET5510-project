function XVy = XV2XVy(XV)
% Compute meridional derivative XVy
% From Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf
% Page 10 (3D zonal wind field)
global jj kk ll dy

XVy = zeros(ll,1);

for k = 1:kk+1
    for j = 2:jj
        l = jk2l(j,k);
        lnh = jk2l(j+1,k);
        lsh = jk2l(j-1,k);

        if j == 2
            XVsh = 0;
        else
            XVsh = XV(lsh);
        end

        if j == jj
            XVnh = 0;
        else
            XVnh = XV(lnh);
        end

        XVy(l) = (XVnh - XVsh) / (2*dy);
    end
end
end