%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: jk2l.m

% Description: Converts 2D grid indices (j, k) representing latitude and height 
% to a single linear index (l) in the quasi-geostrophic model.

% Input:
% - j: Latitude index
% - k: Height index

% Output:
% - l: Linear index

% Math/functions: l = (k-1) * (jj-1) + (j-1)

% - Variables:
%   - jj is the number of latitude points
%   - j is the latitude index (1 to jj-1 for interior)
%   - k is the height index (1 to kk+1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function l = jk2l(j,k)
    global jj  
    
    l = j-1 + (k-1)*(jj-1);
    
    end