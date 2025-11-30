%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: PV2intgrad.m

% DESCRIPTION: Computes the interior potential vorticity (PV) gradient for the 
% Hoskins-West Model in the quasi-geostrophic framework. This function calculates 
% the meridional gradient of PV at interior grid points, incorporating the 
% planetary beta effect, meridional variations of the mean flow, and vertical 
% shear effects. The resulting gradients are essential for analyzing wave 
% propagation, stability, and dynamics within the model atmosphere.

% INPUT:
% - params: Structure containing model parameters (from hwme_config.m)
%   - jj: Number of latitude grid points (integer)
%   - kk: Number of height grid points (integer)
%   - dy: Grid spacing in the meridional direction (m)
%   - dz: Grid spacing in the vertical direction (m)
%   - gamma: Nondimensional parameter for vertical shear variation
%   - HH: Scale height of the atmosphere (m)
%   - m_y: Meridional wavenumber (m^-1)
%   - y0: Reference latitude for the mean flow (m)
%   - Lambda: Amplitude of the mean zonal wind (m/s)
%   - mu: Nondimensional parameter for quadratic vertical shear
%   - f0: Coriolis parameter (s^-1)
%   - NN2: Brunt-Vaisala frequency squared (s^-2)
%   - beta: Planetary vorticity gradient (s^-1 m^-1)

% OUTPUT:
% - dPVdy_interior: 2D array of interior PV gradients (s^-1) with dimensions 
%   (jj+1, kk+1), representing the meridional PV gradient at 
%   interior grid points

% MATH/FUNCTIONS: 
% - dPVdy_interior = β + ∂²U/∂y² - (f₀²/N²) * ∂²U/∂z², where
%   - β is the planetary vorticity gradient
%   - ∂²U/∂y² = -Λ * H * m_y² * sinh(γ * z / H) / sinh(γ) * cos(m_y * (y - y₀))
%   - ∂²U/∂z² = Λ * H * (1/H - μ * 2 * z / H² + γ * cosh(γ * z / H) / sinh(γ) * cos(m_y * (y - y₀)))
%   - For small γ, sinh(γ * z / H) / sinh(γ) ≈ z/H and cosh(γ * z / H) / sinh(γ) ≈ 1

% VARIABLES:
% - Λ: Mean zonal wind amplitude
% - H: Scale height
% - m_y: Meridional wavenumber
% - γ: Shear parameter
% - μ: Quadratic shear parameter
% - f₀: Coriolis parameter
% - N²: Brunt-Vaisala frequency squared
% - β: Planetary vorticity gradient
% - y: Meridional coordinate
% - y₀: Reference latitude
% - z: Vertical coordinate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dPVdy_interior = hwme_PV2intgrad(params)

    %% Extract parameters
    jj = params.jj;
    kk = params.kk;
    dy = params.dy;
    dz = params.dz;
    gamma = params.gamma;
    HH = params.HH;
    m_y = params.m_y;
    y0 = params.y0;
    Lambda = params.Lambda;
    mu = params.mu;
    f0 = params.f0;
    NN2 = params.NN2;
    beta = params.beta;
    
    %% Initialize output array
    dPVdy_interior = zeros(jj+1, kk+1);
    
    %% Compute d(PVbar)/dy for interior points
    for j = 1:jj+1
        y = (j-1) * dy;
        for k = 1:kk+1
            z = (k-1) * dz;
            if abs(gamma) > 1e-10 % avoid division by zero
                sinh_term = sinh(gamma * z / HH) / sinh(gamma);
                cosh_term = cosh(gamma * z / HH) / sinh(gamma);
            else
                sinh_term = z / HH; % approximation for small gamma
                cosh_term = 1;
            end
            d2U_dy2 = -Lambda * HH * m_y^2 * sinh_term * cos(m_y * (y - y0)); % second y derivative
            d2U_dz2 = Lambda * HH * (1/HH - mu * 2 * z / HH^2 + gamma * cosh_term * cos(m_y * (y - y0))); % second z derivative
            dPVdy_interior(j,k) = beta + d2U_dy2 - (f0^2 / NN2) * d2U_dz2; %
        end
    end
end