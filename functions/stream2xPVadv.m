%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: stream2xPVadv.m

% DESCRIPTION: Calculates the zonal advection of potential vorticity (PV) based 
% on the potential vorticity field (QV) in the linear quasi-geostrophic (QG) 
% model. This function is used in matrix construction for stability analysis, 
% computing the contribution of the mean zonal wind to the advection of PV 
% perturbations in the zonal direction.

% INPUT:
% - QV: 1D array representing the potential vorticity field (s^-1), with 
%       length ll, derived from the streamfunction

% OUTPUT:
% - xQVadv: 1D array representing the zonal PV advection (s^-2), with length ll

% MATH/FUNCTIONS: 
% - xQVadv = -u * ∂q/∂x

% VARIABLES:
% - u = Ubar(j,k) is the mean zonal wind at grid point (j,k) (m/s)
% - ∂q/∂x = i * (2π * m0 / Lx) * QV is the zonal derivative of the PV field, 
%   computed spectrally using the zonal wavenumber
% - m0: Zonal wavenumber (dimensionless)
% - Lx: Domain length in the zonal direction (m)
% - cplx: Complex unit (i = sqrt(-1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xQVadv = stream2xPVadv(QV)
    global jj kk ll Lx Ubar cplx m0
    
    xQVadv = zeros(ll,1);
    
    %% Compute zonal PV advection over interior grid points
    for k = 1:kk + 1
        for j = 2:jj
            l = jk2l(j,k); % map (j,k) to linear index l
            xQVadv(l) = -Ubar(j,k) * cplx * (2*pi*m0/Lx) * QV(l); % zonal PV advection
        end
    end

end