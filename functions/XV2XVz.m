%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: XV2XVz.m

% DESCRIPTION: Computes the z-derivative (vertical) of the 1D eigenvector (XV) 
% to derive the vertical component of the streamfunction gradient in the 
% quasi-geostrophic model. This function is used to calculate temperature or 
% vertical velocity fields, critical for analyzing thermal structures and vertical 
% motion in wave dynamics and perturbation evolution within the QG framework.

% INPUT:
% - XV: 1D array representing the streamfunction eigenvector (m²/s), with length ll

% OUTPUT:
% - XVz: 1D array representing the vertical derivative of the streamfunction 
%        (m/s), with length ll

% MATH/FUNCTIONS: 
% - XVz = ∂XV/∂z

% VARIABLES:
% - Interior points: ∂XV/∂z = (XV(j,k+1) - XV(j,k-1)) / (2 * dz)
% - Bottom boundary (k=1): ∂XV/∂z = (XV(j,2) - XV(j,1)) / dz
% - For top boundary (k=kk+1): ∂XV/∂z = (XV(j,kk+1) - XV(j,kk)) / dz
% - dz: Vertical grid spacing (m)
% - The derivative is computed using centered finite differences for interior 
%   points and one-sided differences at boundaries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XVz = XV2XVz(XV)

    global jj kk ll dz
    XVz = zeros(ll,1); % initialize output array
    
    %% Interior points
    for k = 2:kk
        for j = 2:jj
            l = jk2l(j,k);
            lup = jk2l(j,k+1);
            ldw = jk2l(j, k-1);
            XVz(l) = (XV(lup)-XV(ldw))/2/dz;
        end
    end
    
    %% Bottom boundary k=1
    k = 1;
    for j = 2:jj
        l = jk2l(j,k);
        lup = jk2l(j, k+1);
        ldw = l;
        XVz(l) = (XV(lup) - XV(ldw))/dz;
    end
    
    %% Top boundary k=kk+1
    k = kk+1;
    for j = 2:jj
        l = jk2l(j,k);
        ldw = jk2l(j, k-1);
        lup=l;
        XVz(l) = (XV(lup) - XV(ldw))/dz;
    end
end