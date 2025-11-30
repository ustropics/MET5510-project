%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: w2field.m

% DESCRIPTION: Reconstructs the full 3D vertical velocity field from the 
% interior solution vector w_vec by expanding it into the zonal direction 
% using the Fourier mode assumption in the quasi-geostrophic model.

% INPUT:
% - w: 1D array of interior vertical velocity (m/s), length LW
% - ii: Number of grid points in zonal direction
% - dx: Zonal grid spacing (m)

% OUTPUT:
% - wfield: 3D array of vertical velocity (m/s), size (ii+1, jj+1, kk+1)

% MATH/FUNCTIONS: 
% - wfield(i,j,k) = Re[ w(l) * exp(i k x_i) ]

% VARIABLES:
% - k = 2π m0 / Lx: zonal wavenumber (m⁻¹)
% - x_i = (i-1) * dx: zonal coordinate (m)
% - m0: Zonal wavenumber (dimensionless)
% - Lx: Zonal domain length (m)
% - cplx: Complex unit i = √(-1)
% - lw2jk: Function mapping linear index l → (j,k)
% - Only interior (j=1 to jj+1, k=2 to kk) points are filled; 
%   boundary layers (k=1, k=kk+1) remain zero unless extended elsewhere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wfield = w2field(w,ii,dx)
    global jj kk ll cplx m0 Lx LW

    wfield = zeros(ii+1,jj+1,kk+1);

    %% Expand interior solution into full 3D field
    for l = 1:LW
        [j,k] = lw2jk(l);

        % Fill in zonal expansion
        for i = 1:ii+1
            xlon = (i-1)*dx; % zonal coordinate
            wfield(i,j,k) = real( w(l) * exp(cplx*2*pi*m0*xlon/Lx) ); %
        end
    end
end