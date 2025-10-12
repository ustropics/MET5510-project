%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: stream2xPVadv.m

% Description: Calculates the zonal advection of potential vorticity (PV) 
% based on the streamfunction (XV), used in the linear QG model matrix 
% construction.

% Input:
% - QV: 1D array representing the potential vorticity field

% Output:
% - xQVadv: 1D array representing the zonal PV advection

% Math/functions: -u ∂q/∂x, where 
% u = -∂ψ/∂y
% ∂q/∂x is computed via finite differences on the PV field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xQVadv = stream2xPVadv(QV)
    global jj kk ll Lx Ubar cplx m0
    
    xQVadv = zeros(ll,1);
    
    for k = 1:kk + 1
        for j = 2:jj
            l = jk2l(j,k);
            xQVadv(l) = -Ubar(j,k) * cplx * (2*pi*m0/Lx) * QV(l);
        end
    end
    end