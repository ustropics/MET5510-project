%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: ubarCalc.m

% Description: Computes the mean zonal wind field for the modified Eady Model 
% in the quasi-geostrophic framework. This function calculates the background 
% zonal wind profile across latitude and height grids, incorporating linear 
% shear and meridional modulation, essential for stability analysis.

% Input:
% - params: Structure containing model parameters (from hwm_config.m)
%   - jj: Number of latitude grid points (integer)
%   - kk: Number of height grid points (integer)
%   - dy: Grid spacing in the meridional direction (m)
%   - dz: Grid spacing in the vertical direction (m)
%   - HH: Scale height of the atmosphere (m)
%   - DeltaT_hor: Horizontal temperature difference (K)
%   - f0: Coriolis parameter (s^-1)
%   - Theta0: Reference potential temperature (K)
%   - Ly: Meridional domain length (m)
%   - y0: Reference latitude for the mean flow (m)
%   - U0: Reference zonal wind speed (m/s)

% Output:
% - Ubar: 2D array of mean zonal wind (m/s) with dimensions (jj+1, kk+1), 
%         representing the background wind field across latitude and height

% Math/functions: Ubar = [(g * DeltaT_hor / (f0 * Theta0 * Ly)) * z + U0] * cos^4(pi * (y - y0) / Ly)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ubar = hwe_ubar(params)

    %% Extract parameters
    jj = params.jj;
    kk = params.kk;
    dy = params.dy;
    dz = params.dz;
    HH = params.HH;
    DeltaT_hor = params.DeltaT_hor;
    f0 = params.f0;
    Theta0 = params.Theta0;
    Ly = params.Ly;
    y0 = params.y0; % Use y0 as the center latitude (adjust if y_s differs)
    U0 = params.U0;
    g = params.gg; % Assuming gg is defined as gravitational acceleration

    %% Initialize output array
    Ubar = zeros(jj+1, kk+1);

    %% Compute Ubar using modified Eady formula
    for j = 1:jj+1
        y = (j-1)*dy; % y coordinate
        for k = 1:kk+1
            z = (k-1)*dz; % z coordinate
            % Shear term from thermal wind: (g * DeltaT_hor / (f0 * Theta0 * Ly)) * z
            shear_term = (g * DeltaT_hor / (f0 * Theta0 * Ly)) * z;
            % Meridional modulation: cos^4(pi * (y - y0) / Ly)
            meridional_mod = cos(pi * (y - y0) / Ly)^4;
            Ubar(j,k) = (shear_term + U0) * meridional_mod;
        end
    end
end