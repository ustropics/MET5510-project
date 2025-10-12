%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: PV2bndgrad.m

% Description: Computes the potential vorticity (PV) gradients at the surface 
% (z=0) and tropopause (z=HH) for the Hoskins-West Model in the quasi-geostrophic 
% framework. These gradients are critical for boundary condition enforcement in 
% stability analysis and wave dynamics, capturing the meridional variation of PV 
% due to the vertical shear of the mean zonal wind at the boundaries.

% Input:
% - params: Structure containing model parameters (from hoskins_config.m)
%   - jj: Number of latitude grid points (integer)
%   - dy: Grid spacing in the meridional direction (m)
%   - gamma: Nondimensional parameter for vertical shear variation
%   - HH: Scale height of the atmosphere (m)
%   - m_y: Meridional wavenumber (m^-1)
%   - y0: Reference latitude for the mean flow (m)
%   - Lambda: Amplitude of the mean zonal wind (m/s)
%   - mu: Nondimensional parameter for quadratic vertical shear
%   - f0: Coriolis parameter (s^-1)
%   - NN2: Brunt-Vaisala frequency squared (s^-2)

% Output:
% - dPVdy_surf: 1D array of surface PV gradients (s^-1) with length jj+1, 
%               representing the meridional PV gradient at z=0
% - dPVdy_trop: 1D array of tropopause PV gradients (s^-1) with length jj+1, 
%               representing the meridional PV gradient at z=HH

% - Math/functions: 
%   - dPVdy_surf = -(f₀²/N²) * ∂U/∂z at z=0, where 
%       ∂U/∂z = Λ * H * (1 + γ * cosh(0) * cos(m_y * (y - y₀)) / sinh(γ))
%   - dPVdy_trop = (f₀²/N²) * ∂U/∂z at z=HH, where 
%       ∂U/∂z = Λ * H * (1 - μ + γ * cosh(γ) * cos(m_y * (y - y₀)) / sinh(γ))
%   - For small γ, ∂U/∂z is approximated as Λ * H at z=0 and Λ * H * (1 - μ) at z=HH

% - Variables:
%   - f₀: Coriolis parameter
%   - N²: Brunt-Vaisala frequency squared
%   - Λ: Mean zonal wind amplitude
%   - H: Scale height
%   - γ: Shear parameter
%   - μ: Quadratic shear parameter
%   - m_y: Meridional wavenumber
%   - y: Meridional coordinate
%   - y₀: Reference latitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dPVdy_surf, dPVdy_trop] = PV2bndgrad(params)

    %% Extract parameters
    jj = params.jj;
    dy = params.dy;
    gamma = params.gamma;
    HH = params.HH;
    m_y = params.m_y;
    y0 = params.y0;
    Lambda = params.Lambda;
    mu = params.mu;
    f0 = params.f0;
    NN2 = params.NN2;
    
    %% Initialize output arrays
    dPVdy_surf = zeros(1, jj+1);
    dPVdy_trop = zeros(1, jj+1);
    
    %% Compute PV gradients at surface (z=0) and tropopause (z=HH)
    for j = 1:jj+1
        y = (j-1) * dy;
        if abs(gamma) > 1e-10 % avoid division by zero
            sinh_0 = sinh(0) / sinh(gamma); % 0 at surface
            sinh_H = sinh(gamma) / sinh(gamma); % 1 at tropopause
            dU_dz_surf = Lambda * HH * (1 + gamma * cosh(0) * cos(m_y * (y - y0)) / sinh(gamma)); 
            dU_dz_trop = Lambda * HH * (1 - mu + gamma * cosh(gamma) * cos(m_y * (y - y0)) / sinh(gamma)); 
        else
            dU_dz_surf = Lambda * HH * 1; % z=0: 1 + 0/HH - mu * (0/HH)^2
            dU_dz_trop = Lambda * HH * (1 - mu); % z=HH: 1 - mu
        end
        f02_over_N2 = f0^2 / NN2;
        dPVdy_surf(j) = -f02_over_N2 * dU_dz_surf;
        dPVdy_trop(j) = f02_over_N2 * dU_dz_trop;
    end
end