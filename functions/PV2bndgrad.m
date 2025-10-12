function [dPVdy_surf, dPVdy_trop] = PV2bndgrad(params)
% COMPUTE_BOUNDARY_PV_GRADIENTS Computes surface and tropopause PV gradients
% Input: params - Structure containing model parameters (from hoskins_config.m)
% Output: dPVdy_surf - 1D array of surface PV gradients (s^-1)
%         dPVdy_trop - 1D array of tropopause PV gradients (s^-1)

% Extract parameters
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

% Initialize output arrays
dPVdy_surf = zeros(1, jj+1);
dPVdy_trop = zeros(1, jj+1);

% Compute PV gradients at surface (z=0) and tropopause (z=HH)
for j = 1:jj+1
    y = (j-1) * dy;
    if abs(gamma) > 1e-10
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