%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: PV2bndgrad.m

% Description: Computes the potential vorticity (PV) gradients at the surface 
% (z=0) and tropopause (z=HH) for the modified Hoskins-West Eady-type Model 
% in the quasi-geostrophic framework. These gradients are critical for boundary 
% condition enforcement in stability analysis and wave dynamics, capturing the 
% meridional variation of PV due to the vertical shear of the mean zonal wind 
% at the boundaries in an f-plane geometry.

% Input:
% - params: Structure containing model parameters (from hwe_config.m)
%   - jj: Number of latitude grid points (integer)
%   - dy: Grid spacing in the meridional direction (m)
%   - HH: Scale height of the atmosphere (m)
%   - y0: Reference latitude for the mean flow (m)
%   - DeltaT_hor: Horizontal temperature difference (K)
%   - f0: Coriolis parameter (s^-1)
%   - NN2: Brunt-Vaisala frequency squared (s^-2)
%   - Theta0: Reference potential temperature (K)
%   - Ly: Meridional domain length (m)
%   - U0: Reference zonal wind speed (m/s)

% Output:
% - dPVdy_surf: 1D array of surface PV gradients (s^-1) with length jj+1, 
%               representing the meridional PV gradient at z=0
% - dPVdy_trop: 1D array of tropopause PV gradients (s^-1) with length jj+1, 
%               representing the meridional PV gradient at z=HH

% - Math/functions: 
%   - dPVdy_surf = -(f₀²/N²) * ∂u/∂z at z=0, where 
%       ∂u/∂z = (g * DeltaT_hor / (f0 * Theta0 * Ly)) * cos^4(pi * (y - y0) / Ly)
%   - dPVdy_trop = (f₀²/N²) * ∂u/∂z at z=HH, where 
%       ∂u/∂z is the same as above (constant with height)

% - Variables:
%   - f₀: Coriolis parameter
%   - N²: Brunt-Vaisala frequency squared
%   - g: Gravitational acceleration
%   - DeltaT_hor: Horizontal temperature difference
%   - Theta0: Reference potential temperature
%   - Ly: Meridional domain length
%   - y: Meridional coordinate
%   - y₀: Reference latitude
%   - U0: Reference zonal wind speed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dPVdy_surf, dPVdy_trop] = PV2bndgrad(params)

    %% Extract parameters
    jj = params.jj;
    dy = params.dy;
    HH = params.HH;
    y0 = params.y0;
    DeltaT_hor = params.DeltaT_hor;
    f0 = params.f0;
    NN2 = params.NN2;
    Theta0 = params.Theta0;
    Ly = params.Ly;
    U0 = params.U0;
    g = params.gg; % Assuming gg is defined as gravitational acceleration
    
    %% Initialize output arrays
    dPVdy_surf = zeros(1, jj+1);
    dPVdy_trop = zeros(1, jj+1);
    
    %% Compute PV gradients at surface (z=0) and tropopause (z=HH)
    for j = 1:jj+1
        y = (j-1) * dy;
        % Vertical shear: ∂u/∂z = (g * DeltaT_hor / (f0 * Theta0 * Ly)) * cos^4(pi * (y - y0) / Ly)
        dU_dz = (g * DeltaT_hor / (f0 * Theta0 * Ly)) * cos(pi * (y - y0) / Ly)^4;
        f02_over_N2 = f0^2 / NN2;
        % Surface PV gradient (negative due to quasi-geostrophic convention)
        dPVdy_surf(j) = -f02_over_N2 * dU_dz;
        % Tropopause PV gradient (positive, opposite sign)
        dPVdy_trop(j) = f02_over_N2 * dU_dz;
    end
end