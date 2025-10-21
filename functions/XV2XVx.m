%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: XV2XVx.m

% DESCRIPTION: Computes the x-derivative (zonal) of the 1D eigenvector (XV) to 
% derive the zonal component of the streamfunction gradient in the quasi-geostrophic 
% model. This function is used to calculate the meridional wind component or related 
% fields, essential for analyzing wave dynamics and perturbation evolution in the 
% linear QG framework.

% INPUT:
% - XV: 1D array representing the streamfunction eigenvector (m²/s), with length ll

% OUTPUT:
% - XVx: 1D array representing the zonal derivative of the streamfunction (m/s)

% MATH/FUNCTIONS: XVx = ∂XV/∂x = i * (2π * m0 / Lx) * XV

% - VARIABLES:
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
    
    % Compute zonal derivative using spectral method
    for k = 1:kk+1
        for j = 2:jj
            l=jk2l(j,k);
            XVx(l) = cplx*(2*pi*m0/Lx) * XV(l); % spectral derivative in x
        end
    end

end