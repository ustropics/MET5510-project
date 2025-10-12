%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: ubarCalc.m

% Description: Computes mean zonal wind field for the Hoskins-West Model

% Input:
% - params: Structure containing model parameters (from hoskins_config.m)

% Output:
% - Ubar: 2D array of mean zonal wind (m/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ubar = ubarCalc(params)
    % Extract parameters
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
    
    % Initialize output array
    Ubar = zeros(jj+1, kk+1);
    
    % Compute Ubar using Modified Hoskins-West formula
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