%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: XV2XVy.m

% DESCRIPTION: Computes the y-derivative (meridional) of the 1D eigenvector (XV) 
% to derive the meridional component of the streamfunction gradient in the 
% quasi-geostrophic model. This function is used to calculate the zonal wind 
% component or related fields, critical for analyzing wave dynamics and 
% perturbation evolution within the QG framework.

% INPUT:
% - XV: 1D array representing the streamfunction eigenvector (m²/s), with length ll

% OUTPUT:
% - XVy: 1D array representing the meridional derivative of the streamfunction 
%        (m/s), with length ll

% MATH/FUNCTIONS: 
% - XVy = ∂XV/∂y

% VARIABLES:
% - Interior points: ∂XV/∂y = (XV(j+1,k) - XV(j-1,k)) / (2 * dy)
% - Southern boundary (j=2): ∂XV/∂y = (XV(j+1,k) - 0) / (2 * dy)  [since XV=0 at j=1 (ghost)]
% - Northern boundary (j=jj): ∂XV/∂y = (0 - XV(j-1,k)) / (2 * dy)  [since XV=0 at j=jj+1 (ghost)]
% - dy: Meridional grid spacing (m)
% - The derivative is computed using centered finite differences for interior 
%   points, with Dirichlet boundary conditions XV=0 enforced at the southern 
%   (j=1, ghost) and northern (j=jj+1, ghost) walls via zero values at j=2 and j=jj 
%   when evaluating the stencil at the boundary-adjacent interior points.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XVy = XV2XVy (XV)
    global jj kk ll dy
    
    XVy = zeros(ll,1);
    
    %% compute meridional derivative using finite differences
    for k = 1:kk+1
        for j = 2:jj % loop over meridional grid points
            l = jk2l(j,k); % linear index for (j,k)
            lnh = jk2l(j+1,k); % index for (j+1,k)
            lsh = jk2l(j-1,k); % index for (j-1,k)
    
            if(j == 2)
                XVsh=0; % southern boundary condition
            else
                XVsh=XV(lsh); % southern level
            end
    
            if(j == jj)
                XVnh = 0; % northern boundary condition
            else
                XVnh = XV(lnh); % northern level
            end
    
            XVy(l) = ( XVnh - XVsh) /2/dy; % centered finite difference
        end
    end

end