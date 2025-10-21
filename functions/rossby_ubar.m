%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: rossby_ubar.m

% DESCRIPTION: Initializes the mean zonal wind field (Ubar) for the Rossby 
% wave model in the quasi-geostrophic framework, set to a constant value.

% INPUT:
% - jj: Number of latitude grid points
% - kk: Number of height grid points

% OUTPUT:
% - Ubar: Mean zonal wind field (jj+1 x kk+1 array, m/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ubar = rossby_ubar(jj, kk)

    Ubar = zeros(jj + 1, kk + 1);
    Ubar(:,:) = Ubar(:,:) + 10; % set uniform zonal wind speed (constant 10 m/s)

end