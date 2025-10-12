%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: ubarCalc.m

% Description: Computes the mean zonal wind field for the Hoskins-West Model 
% in the quasi-geostrophic framework. This function calculates the background 
% zonal wind profile across latitude and height grids, incorporating linear, 
% quadratic, and oscillatory components of the wind field, which are essential 
% for stability analysis and wave dynamics in the model atmosphere.

% Input:
% - params: Structure containing model parameters (from hoskins_config.m)
%   - jj: Number of latitude grid points (integer)
%   - kk: Number of height grid points (integer)
%   - dy: Grid spacing in the meridional direction (m)
%   - dz: Grid spacing in the vertical direction (m)
%   - HH: Scale height of the atmosphere (m)
%   - mu: Nondimensional parameter for quadratic vertical shear
%   - gamma: Nondimensional parameter for vertical shear variation
%   - m_y: Meridional wavenumber (m^-1)
%   - y0: Reference latitude for the mean flow (m)
%   - Lambda: Amplitude of the mean zonal wind (m/s)
%   - U0: Reference zonal wind speed (m/s)

% Output:
% - Ubar: 2D array of mean zonal wind (m/s) with dimensions (jj+1, kk+1), 
%         representing the background wind field across latitude and height

% Math/functions: Ubar = Λ * H * (z/H - (μ/2) * (z/H)² + sinh(γ * z / H) / sinh(γ) * cos(m_y * (y - y₀))) - U₀

% - Variables:
%   - Λ: Mean zonal wind amplitude
%   - H: Scale height
%   - μ: Quadratic shear parameter
%   - γ: Shear parameter
%   - m_y: Meridional wavenumber
%   - y: Meridional coordinate
%   - y₀: Reference latitude
%   - z: Vertical coordinate
%   - U₀: Reference zonal wind speed
%   - sinh term is set to 0 for small γ to avoid numerical issues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ubar = hwm_ubar(params)

    %% Extract parameters
    jj = params.jj;
    kk = params.kk;
    dy = params.dy;
    dz = params.dz;
    HH = params.HH;
    mu = params.mu;
    gamma = params.gamma;
    m_y = params.m_y;
    y0 = params.y0;
    Lambda = params.Lambda;
    U0 = params.U0;
    
    %% Initialize output array
    Ubar = zeros(jj+1, kk+1);
    
    %% Compute Ubar using Modified Hoskins-West formula
    for j = 1:jj+1
        y = (j-1)*dy; % y coordinate
        for k = 1:kk+1
            z = (k-1)*dz; % z coordinate
            term_linear = z / HH;
            term_quadratic = -mu / 2 * (z / HH)^2;
            if abs(gamma) > 1e-10  % avoid division by zero
                term_sinh = sinh(gamma * z / HH) / sinh(gamma) * cos(m_y * (y - y0));
            else
                term_sinh = 0;
            end
            Ubar(j,k) = Lambda * HH * (term_linear + term_quadratic + term_sinh) - U0;
        end
    end
end