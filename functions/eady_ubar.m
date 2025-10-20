%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: init_ubar.m

% Description: Initializes the mean zonal wind field (Ubar) for the Rossby 
% wave model in the quasi-geostrophic framework, based on a linear vertical 
% profile derived from the potential temperature gradient.

% Input:
% - jj: Number of latitude grid points
% - kk: Number of height grid points
% - gg: Gravitational acceleration (m s^-2)
% - f0: Coriolis parameter (s^-1)
% - Theta0: Reference potential temperature (K)
% - dTbar: Potential temperature gradient (K)
% - Ly: Meridional domain length (m)
% - ZZ: Vertical grid coordinates (m)
% - Ubar0: Base zonal wind speed (m/s)

% Output:
% - Ubar: Mean zonal wind field (jj+1 x kk+1 array, m/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ubar = eady_ubar(jj, kk, gg, f0, Theta0, dTbar, Ly, ZZ, Ubar0)
Ubar = zeros(jj + 1, kk + 1);
for k = 1:kk + 1
    Ubar(:, k) = (gg / f0 / Theta0) * (dTbar / Ly) * ZZ(k) + Ubar0;
end
end