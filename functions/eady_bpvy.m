%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: init_bpvy.m

% Description: Initializes the potential vorticity gradient field (BPVy) for 
% the Rossby wave model in the quasi-geostrophic framework, incorporating 
% beta-plane dynamics and second derivatives of the mean zonal wind.

% Input:
% - jj: Number of latitude grid points
% - kk: Number of height grid points
% - beta: Planetary vorticity gradient (m^-1 s^-1)
% - Ubar: Mean zonal wind field (jj+1 x kk+1 array, m/s)
% - dy: Meridional grid spacing (m)
% - f0: Coriolis parameter (s^-1)
% - NN2: Brunt-Vaisala frequency squared (s^-2)
% - dz: Vertical grid spacing (m)

% Output:
% - BPVy: Potential vorticity gradient field (jj+1 x kk+1 array, m^-1 s^-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BPVy = init_bpvy(jj, kk, beta, Ubar, dy, f0, NN2, dz)
BPVy = zeros(jj + 1, kk + 1);
for j = 2:jj
    for k = 2:kk
        BPVy(j, k) = beta - (Ubar(j + 1, k) - 2 * Ubar(j, k) + Ubar(j - 1, k)) / dy^2 ...
            - (f0 * f0 / NN2) * (Ubar(j, k + 1) - 2 * Ubar(j, k) + Ubar(j, k - 1)) / dz^2;
    end
    k = 1;
    BPVy(j, k) = -(Ubar(j, k + 1) - Ubar(j, k)) / dz;
    k = kk + 1;
    BPVy(j, k) = -(Ubar(j, k) - Ubar(j, k - 1)) / dz;
end
end