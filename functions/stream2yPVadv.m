%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: stream2yPVadv.m

% DESCRIPTION: Computes the meridional advection of potential vorticity (PV) 
% including the beta effect, based on the streamfunction (XV), for the linear 
% quasi-geostrophic (QG) model. This function is used in matrix construction 
% for stability analysis, calculating the contribution of the meridional wind 
% and planetary vorticity gradient to the advection of PV perturbations in the 
% meridional direction.

% INPUT:
% - XV: 1D array representing the streamfunction eigenvector (m²/s)

% OUTPUT:
% - yQVadv: 1D array representing the meridional PV advection (s^-2)

% MATH/FUNCTIONS: 
% - yQVadv = -β * ∂ψ/∂x

% VARIABLES:
% - β = BPVy(j,k) is the meridional PV gradient at grid point (j,k) (s^-1 m^-1)
% - ∂ψ/∂x = i * (2π * m0 / Lx) * XV is the zonal derivative of the streamfunction, 
%   computed spectrally using the zonal wavenumber
% - m0: Zonal wavenumber (dimensionless)
% - Lx: Domain length in the zonal direction (m)
% - cplx: Complex unit (i = sqrt(-1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yQVadv = stream2yPVadv(XV)
    global jj kk ll BPVy Lx cplx m0
    
    yQVadv = zeros(ll,1);
    
    %% Compute meridional PV advection over interior grid points
    for k = 1:kk + 1
        for j = 2:jj
            l = jk2l(j,k); % map (j,k) to linear index l
            yQVadv(l) = -BPVy(j,k) * cplx * (2*pi*m0/Lx) * XV(l); % meridional PV advection
        end
    end
    
end