%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: stream2yPVadv.m

% Description: Computes the meridional advection of potential vorticity (PV) 
% including the beta effect, based on the streamfunction (XV), for the linear 
% quasi-geostrophic (QG) model. This function is used in matrix construction 
% for stability analysis, calculating the contribution of the meridional wind 
% and planetary vorticity gradient to the advection of PV perturbations in the 
% meridional direction.

% Input:
% - XV: 1D array representing the streamfunction eigenvector (m²/s)

% Output:
% - yQVadv: 1D array representing the meridional PV advection (s^-2)

% Math/functions: yQVadv = -β * ∂ψ/∂x
%
% - Variables:
%   - β = BPVy(j,k) is the meridional PV gradient at grid point (j,k) (s^-1 m^-1)
%   - ∂ψ/∂x = i * (2π * m0 / Lx) * XV is the zonal derivative of the streamfunction, 
%     computed spectrally using the zonal wavenumber
%   - m0: Zonal wavenumber (dimensionless)
%   - Lx: Domain length in the zonal direction (m)
%   - cplx: Complex unit (i = sqrt(-1))

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