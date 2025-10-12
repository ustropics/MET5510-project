%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XV2XVy.m

% Description: Computes the y-derivative (meridional) of the 1D eigenvector (XV) 
% to derive the meridional component of the streamfunction gradient in the 
% quasi-geostrophic model. This function is used to calculate the zonal wind 
% component or related fields, critical for analyzing wave dynamics and 
% perturbation evolution in the linear QG framework.

% Input:
% - XV: 1D array representing the streamfunction eigenvector (m²/s), with length ll

% Output:
% - XVy: 1D array representing the meridional derivative of the streamfunction 
%        (m/s), with length ll

% Math/functions: XVy = ∂XV/∂y = (XV(j+1,k) - XV(j-1,k)) / (2 * dy), where
% - dy: Meridional grid spacing (m)
% - XV(j+1,k) and XV(j-1,k) are streamfunction values at neighboring meridional 
%   grid points, with boundary conditions setting XV=0 at j=2 (southern) and 
%   j=jj (northern) boundaries
% - The derivative is computed using centered finite differences for interior 
%   points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XVy = XV2XVy (XV)
    global jj kk ll dy
    
    XVy = zeros(ll,1);
    
    for k = 1:kk+1
        for j = 2:jj
            l = jk2l(j,k);
            lnh = jk2l(j+1,k);
            lsh = jk2l(j-1,k);
    
            if(j == 2)
                XVsh=0;
            else
                XVsh=XV(lsh);
            end
    
            if(j == jj)
                XVnh = 0;
            else
                XVnh = XV(lnh);
            end
    
            XVy(l) = ( XVnh - XVsh) /2/dy;
        end
    end
end