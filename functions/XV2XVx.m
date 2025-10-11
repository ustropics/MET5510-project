%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XV2XVx.m

% Description: Computes the x-derivative (zonal) of the 1D eigenvector (XV) 
% to derive the zonal component of the streamfunction gradient in the 
% quasi-geostrophic model.

% Math/functions: XVx = ∂XV/∂x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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