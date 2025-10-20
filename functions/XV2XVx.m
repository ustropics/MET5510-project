%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XV2XVx.m

% Description: Computes the x-derivative (zonal) of the 1D eigenvector (XV) to 
% derive the zonal component of the streamfunction gradient in the quasi-geostrophic 
% model. This function is used to calculate the meridional wind component or related 
% fields, essential for analyzing wave dynamics and perturbation evolution in the 
% linear QG framework.

% Input:
% - XV: 1D array representing the streamfunction eigenvector (m²/s), with length ll

% Output:
% - XVx: 1D array representing the zonal derivative of the streamfunction (m/s)

% Math/functions: XVx = ∂XV/∂x = i * (2π * m0 / Lx) * XV

% - Variables:
%   - m0: Zonal wavenumber (dimensionless)
%   - Lx: Domain length in the zonal direction (m)
%   - cplx: Complex unit (i = sqrt(-1))
%   - The derivative is computed spectrally using the zonal wavenumber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XVx = XV2XVx (XV)
    global jj kk ll m0 Lx cplx
    
    XVx = zeros(ll, 1);
    
    for k = 1:kk+1
        for j = 2:jj
            l=jk2l(j,k);
            XVx(l) = cplx*(2*pi*m0/Lx) * XV(l);
        end
    end
end