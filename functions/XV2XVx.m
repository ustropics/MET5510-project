% functions/XV2XVx.m
function XVx = XV2XVx (XV)
global jj kk ll m0 Lx cplx

XVx = zeros(ll, 1);

for k = 1:kk+1

    for j = 2:jj
        l=jk2l(j,kk);
        XVx(l) = cplx*(2*pi*m0/Lx) * XV(l);
    end
end
end