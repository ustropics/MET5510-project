%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: jk2l.m

% DESCRIPTION: Converts 2D grid indices (j, k) representing latitude and height 
% to a single linear index (l) in the quasi-geostrophic model.

% INPUT:
% - j: Latitude index
% - k: Height index

% OUTPUT:
% - l: Linear index

% MATH/FUNCTIONS: 
% - l = (k-1) * (jj-1) + (j-1)

% VARIABLES:
% - jj is the number of latitude points
% - j is the latitude index (1 to jj-1 for interior)
% - k is the height index (1 to kk+1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function l = jk2l(j,k)
    
    global jj  
    
    l = j-1 + (k-1)*(jj-1);
    
end