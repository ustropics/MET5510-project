%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: PV2intgrad.m

% Description: Computes interior PV gradient for the Hoskins-West Model

% Input:
% - params: Structure containing model parameters (from hoskins_config.m)

% Output:
% - dPVdy_interior: 2D array of interior PV gradients (s^-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dPVdy_interior = compute_dPVdy_interior(params)
    % Extract parameters
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
    
    % Initialize output array
    dPVdy_interior = zeros(jj+1, kk+1);
    
    % Compute d(PVbar)/dy for interior points
    for j = 1:jj+1
        y = (j-1) * dy;
        for k = 1:kk+1
            z = (k-1) * dz;
            if abs(gamma) > 1e-10
                sinh_term = sinh(gamma * z / HH) / sinh(gamma);
                cosh_term = cosh(gamma * z / HH) / sinh(gamma);
            else
                sinh_term = z / HH; % Approximation for small gamma
                cosh_term = 1;
            end
            d2U_dy2 = -Lambda * HH * m_y^2 * sinh_term * cos(m_y * (y - y0)); % Second y derivative
            d2U_dz2 = Lambda * HH * (1/HH - mu * 2 * z / HH^2 + gamma * cosh_term * cos(m_y * (y - y0))); % Second z derivative
            dPVdy_interior(j,k) = beta + d2U_dy2 - (f0^2 / NN2) * d2U_dz2;
        end
    end
    end