%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: jk2lw.m

% DESCRIPTION: Converts 2D grid indices (j, k) representing latitude and height 
% to a single linear index (l) in the quasi-geostrophic model.

% INPUT:
% - j: Latitude index
% - k: Height index

% OUTPUT:
% - l: Linear index

% MATH/FUNCTIONS: 
% l = j + (k-2)*(jj+1)

% VARIABLES:
% - jj is the number of latitude points
% - j is the latitude index
% - k is the height index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function l = jk2lw(j,k)
   
    global jj  
    
    l = j + (k-2)*(jj+1);
    
end