%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: BPVyCalc.m

% Description: Computes the beta-plane potential vorticity (PV) gradient for 
% the Hoskins-West Model, which is used to represent the meridional gradient 
% of PV in the quasi-geostrophic framework. This function calculates the PV 
% gradient at interior grid points, incorporating the planetary beta effect 
% and the meridional variation of the mean flow, essential for stability 
% analysis and wave propagation studies.

% Input:
% - params: Structure containing model parameters (from hoskins_config.m)
%   - jj: Number of latitude grid points (integer)
%   - kk: Number of height grid points (integer)
%   - dy: Grid spacing in the meridional direction (m)
%   - dz: Grid spacing in the vertical direction (m)
%   - gamma: Nondimensional parameter for vertical shear variation
%   - HH: Scale height of the atmosphere (m)
%   - m_y: Meridional wavenumber (m^-1)
%   - y0: Reference latitude for the mean flow (m)
%   - Lambda: Amplitude of the mean zonal wind (m/s)
%   - beta: Planetary vorticity gradient (s^-1 m^-1)

% Output:
% - BPVy: 2D array of PV gradient values (s^-1) with dimensions (jj+1, kk+1), 
%         representing the meridional PV gradient at interior grid points

% Math/functions: BPVy = β - 2 * Λ * H * m_y^3 * sinh(γ * z / H) / sinh(γ) * sin(m_y * (y - y0))

% - Variables:
%   - β is the planetary vorticity gradient
%   - Λ is the amplitude of the mean zonal wind
%   - H is the scale height
%   - m_y is the meridional wavenumber
%   - γ is the nondimensional shear parameter
%   - y is the meridional coordinate
%   - y0 is the reference latitude
%   - z is the vertical coordinate
%   - sinh_term is approximated as z/H for small γ to avoid numerical issues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BPVy = BPVyCalc(params)
    
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
    beta = params.beta;
    
    %% Initialize output array
    BPVy = zeros(jj+1, kk+1);
    
    %% Compute BPVy for interior points
    for j = 2:jj % skip boundaries
        y = (j-1)*dy; % meridional coordinate
        for k = 2:kk
            z = (k-1)*dz; % vertical coordinate
            
            if abs(gamma) > 1e-10
                sinh_term = sinh(gamma * z / HH) / sinh(gamma); % avoid division by zero
            else
                sinh_term = 0; % for very small gamma, sinh(gamma*z/HH)/sinh(gamma) ~ z/HH
            end

            BPVy(j,k) = beta - 2 * Lambda * HH * m_y^3 * sinh_term * sin(m_y * (y - y0));
        end
    end
end