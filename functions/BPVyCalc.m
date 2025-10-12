function BPVy = BPVyCalc(params)
% COMPUTE_BPVY Computes beta-plane PV gradient for Hoskins-West Model
% Input: params - Structure containing model parameters (from hoskins_config.m)
% Output: BPVy - 2D array of PV gradient (s^-1)

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
beta = params.beta;

% Initialize output array
BPVy = zeros(jj+1, kk+1);

% Compute BPVy for interior points
for j = 2:jj
    y = (j-1)*dy;
    for k = 2:kk
        z = (k-1)*dz;
        if abs(gamma) > 1e-10
            sinh_term = sinh(gamma * z / HH) / sinh(gamma);
        else
            sinh_term = 0;
        end
        BPVy(j,k) = beta - 2 * Lambda * HH * m_y^3 * sinh_term * sin(m_y * (y - y0));
    end
end
end