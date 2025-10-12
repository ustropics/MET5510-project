%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: BPVyCalc.m

% Description: Computes the meridional potential vorticity (PV) gradient for 
% the modified Hoskins-West Eady-type Model in the quasi-geostrophic framework. 
% This function calculates the PV gradient at interior grid points, incorporating 
% the meridional variation of the mean flow, essential for stability analysis 
% and wave propagation studies in an f-plane geometry.

% Input:
% - params: Structure containing model parameters (from hwe_config.m)
%   - jj: Number of latitude grid points (integer)
%   - kk: Number of height grid points (integer)
%   - dy: Grid spacing in the meridional direction (m)
%   - dz: Grid spacing in the vertical direction (m)
%   - HH: Scale height of the atmosphere (m)
%   - y0: Reference latitude for the mean flow (m)
%   - DeltaT_hor: Horizontal temperature difference (K)
%   - f0: Coriolis parameter (s^-1)
%   - NN2: Brunt-Vaisala frequency squared (s^-2)
%   - Theta0: Reference potential temperature (K)
%   - Ly: Meridional domain length (m)
%   - U0: Reference zonal wind speed (m/s)

% Output:
% - BPVy: 2D array of PV gradient values (s^-1) with dimensions (jj+1, kk+1), 
%         representing the meridional PV gradient at interior grid points

% Math/functions: 
% - BPVy = ∂²u/∂y², where
%   - u = [(g * DeltaT_hor / (f0 * Theta0 * Ly)) * z + U0] * cos^4(pi * (y - y0) / Ly)
%   - ∂²u/∂y² = -4 * pi^2 / Ly^2 * [(g * DeltaT_hor / (f0 * Theta0 * Ly)) * z + U0] * cos^4(pi * (y - y0) / Ly) * 
%               [3 - 4 * cos^2(pi * (y - y0) / Ly)]

% - Variables:
%   - g: Gravitational acceleration
%   - DeltaT_hor: Horizontal temperature difference
%   - f0: Coriolis parameter
%   - Theta0: Reference potential temperature
%   - Ly: Meridional domain length
%   - y: Meridional coordinate
%   - y0: Reference latitude
%   - z: Vertical coordinate
%   - U0: Reference zonal wind speed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BPVy = BPVyCalc(params)
    
    %% Extract parameters
    jj = params.jj;
    kk = params.kk;
    dy = params.dy;
    dz = params.dz;
    HH = params.HH;
    y0 = params.y0;
    DeltaT_hor = params.DeltaT_hor;
    f0 = params.f0;
    NN2 = params.NN2;
    Theta0 = params.Theta0;
    Ly = params.Ly;
    U0 = params.U0;
    g = params.gg; % Assuming gg is defined as gravitational acceleration
    
    %% Initialize output array
    BPVy = zeros(jj+1, kk+1);
    
    %% Compute BPVy for interior points
    for j = 2:jj % skip boundaries
        y = (j-1) * dy; % meridional coordinate
        for k = 2:kk
            z = (k-1) * dz; % vertical coordinate
            % Compute the base wind profile term: (g * DeltaT_hor / (f0 * Theta0 * Ly)) * z + U0
            shear_term = (g * DeltaT_hor / (f0 * Theta0 * Ly)) * z + U0;
            % Meridional modulation: cos^4(pi * (y - y0) / Ly)
            cos_term = cos(pi * (y - y0) / Ly);
            cos4_term = cos_term^4;
            % Second derivative of u with respect to y
            d2u_dy2 = -4 * (pi^2 / Ly^2) * shear_term * cos4_term * (3 - 4 * cos_term^2);
            % PV gradient: β = 0, d²u/dz² = 0, so BPVy = d²u/dy²
            BPVy(j,k) = d2u_dy2;
        end
    end
end