function field = XV2field(XV,ii,dx)
% Convert from streamfunction vector XV to 3D field
% From Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf
% Page 2 (Eigenmode solutions)
global jj kk ll cplx m0 Lx

field = zeros(ii+1, jj+1, kk+1);

for l = 1:ll
   [j,k] = l2jk(l);
   for i = 1:ii+1
        xlon = (i-1)*dx;
        field(i,j,k) = real(XV(l) * exp(cplx*2*pi*m0*xlon/Lx));
    end
end
end