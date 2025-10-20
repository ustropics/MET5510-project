%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: hwme_ubar.m

% Description: Initializes the mean zonal wind field (Ubar) for the Hoskins-West
% modified Eady-type model in the quasi-geostrophic framework, using the formula:
% U(y,z) = (g / (f0 * Theta0)) * (H * ΔT / Ly) *
%          [ (z/H) - (μ/2)*(z/H)^2 + sinh(2πLr/Ly * z/H)/sinh(2πLr/Ly) * cos(2π(y - y_s)/Ly) ]

% Input:
% - jj: Number of latitude grid points
% - kk: Number of height grid points
% - gg: Gravitational acceleration (m s^-2)
% - f0: Coriolis parameter (s^-1)
% - Theta0: Reference potential temperature (K)
% - dTbar: Potential temperature gradient (K)
% - HH: Scale height (m)
% - Ly: Meridional domain length (m)
% - ZZ: Vertical grid coordinates (m)
% - YY: Meridional grid coordinates (m)
% - mu: Curvature parameter
% - y_s: Reference latitude (m)
% - Lr: Rossby radius of deformation (m)

% Output:
% - Ubar: Mean zonal wind field (jj+1 x kk+1 array, m/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ubar = hwme_ubar(jj, kk, gg, f0, Theta0, dTbar, HH, Ly, ZZ, YY, mu, y_s, Lr)
Ubar = zeros(jj + 1, kk + 1);
prefac = (gg / (f0 * Theta0)) * (HH * dTbar / Ly);
A = 2 * pi * Lr / Ly;

for j = 1:(jj + 1)
    y = YY(j);
    for k = 1:(kk + 1)
        zH = ZZ(k) / HH;
        sinh_term = sinh(A * zH) / sinh(A);
        cos_term = cos(2 * pi * (y - y_s) / Ly);
        Ubar(j, k) = prefac * (zH - 0.5 * mu * (zH + sinh_term * cos_term));
    end
end
end