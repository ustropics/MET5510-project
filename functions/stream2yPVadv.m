%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: stream2yPVadv.m

% Description: Computes the meridional advection of potential vorticity (PV) 
% including the beta effect, based on the streamfunction (XV), for the 
% linear QG model matrix construction.

% Math/functions: -v ∂q/∂y + β ∂ψ/∂x, where 
% v = ∂ψ/∂x
% β is the meridional PV gradient
% derivatives are calculated via finite differences

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yQVadv = stream2yPVadv(XV)
global jj kk ll BPVy Lx cplx m0

yQVadv = zeros(ll,1);

for k = 1:kk + 1

    for j = 2:jj

        l = jk2l(j,k);

        yQVadv(l) = -BPVy(j,k) * cplx * (2*pi*m0/Lx) * XV(l);
    end
end
end