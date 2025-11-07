%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: XVz2field.m

% DESCRIPTION: Converts the z-derivative of the eigenvector (XVz) into a 
% 3D vertical velocity or temperature perturbation field over the grid in 
% the quasi-geostrophic model.

% INPUT:
% - XV: Eigenvector representing the streamfunction
% - ii: Number of grid points in the x-direction
% - dx: Grid spacing in the x-direction

% OUTPUT:
% - field: 3D array representing the vertical velocity or temperature perturbation field

% MATH/FUNCTIONS: 
% - w' or T' = (f₀ * H / R) ∂ψ/∂z

% VARIABLES: 
% - ψ is streamfunction
% - f₀ is Coriolis parameter
% - H is scale height
% - R is gas constant
% - ∂/∂z is derived from XVz 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field = XVz2field(XV,ii,dx) 
    global jj kk ll cplx m0 Lx dz

    field = zeros(ii+1,jj+1,kk+1); % initalize 3D field
    
    %% Convert linear index to (j,k)
    for l = 1:ll
        [j,k] = l2jk(l);

        for i = 1:ii
            xlon = (i-1)*dx; % longitude coordinate
            if (k == 1)

                % Difference at bottom boundary
                field(i,j,k) = real((XV(jk2l(j,k+1)) - ... 
                    XV(l))*exp(cplx*2*pi*m0*xlon/Lx))/dz;
           
            elseif (k == kk+1)

                % Backward difference at top boundary
                field(i,j,k) = real((XV(l) - ...
                    XV(jk2l(j,k-1)))*exp(cplx*2*pi*m0*xlon/Lx))/dz;
            else

                % Interior difference 
                field(i,j,k) = real((XV(jk2l(j,k+1)) - ...
                    XV(jk2l(j,k-1)))*exp(cplx*2*pi*m0*xlon/Lx))/(2*dz);
            end
        end
    end

end