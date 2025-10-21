%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: lw2jk.m

% DESCRIPTION: Converts a single linear index (l) to 2D grid indices (j, k) 
% representing latitude and height in the quasi-geostrophic model.

% INPUT:
% - l: Linear index

% OUTPUT:
% - j: Latitude index
% - k: Height index

% MATH/FUNCTIONS: [j, k] = [l - (k-2)*(jj+1), floor((l-1)/(jj+1))+2]

% - VARIABLES:
%   - j is the number of latitude points
%   - l is the linear index
%   - jj is the number of latitude points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [j,k] = lw2jk(l)
    
global jj
    
    k = floor((l-1)/(jj+1))+2;
    j = l - (k-2)*(jj+1);
    
end