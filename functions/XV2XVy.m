%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XV2XVy.m

% Description: Computes the y-derivative (meridional) of the 1D eigenvector 
% (XV) to derive the meridional component of the streamfunction gradient in 
% the quasi-geostrophic model.

% Input:
% - XV: Eigenvector representing the streamfunction

% Output:
% - XVy: 1D array representing the meridional derivative of the streamfunction

% Math/functions: XVy = ∂XV/∂y

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