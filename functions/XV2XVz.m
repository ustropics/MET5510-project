%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XV2XVz.m

% Description: Computes the z-derivative (vertical) of the 1D eigenvector 
% (XV) to derive the vertical component of the streamfunction gradient, 
% used for temperature or vertical velocity in the quasi-geostrophic model.

% Input:
% - XV: Eigenvector representing the streamfunction

% Output:
% - XVz: 1D array representing the vertical derivative of the streamfunction

% Math/functions: XVz = ∂XV/∂z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XVz = XV2XVz(XV)
    global jj kk ll dz
    XVz = zeros(ll,1);
    
    for k = 2:kk
        for j = 2:jj
            l = jk2l(j,k);
            lup = jk2l(j,k+1);
            ldw = jk2l(j, k-1);
            XVz(l) = (XV(lup)-XV(ldw))/2/dz;
        end
    end
    
    k = 1;
    for j = 2:jj
        l = jk2l(j,k);
        lup = jk2l(j, k+1);
        ldw = l;
        XVz(l) = (XV(lup) - XV(ldw))/dz;
    end
    
    k = kk+1;
    for j = 2:jj
        l = jk2l(j,k);
        ldw = jk2l(j, k-1);
        lup=l;
        XVz(l) = (XV(lup) - XV(ldw))/dz;
    end
    end