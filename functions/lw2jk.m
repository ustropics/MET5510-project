%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: lw2jk.m

% Description: Converts a single linear index (l) to 2D grid indices (j, k) 
% representing latitude and height in the quasi-geostrophic model.

% Input:
% - l: Linear index

% Output:
% - j: Latitude index
% - k: Height index

% Math/functions: [j, k] = [l - (k-2)*(jj+1), floor((l-1)/(jj+1))+2], where 
% jj is the number of latitude points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [j,k] = lw2jk(l)
    global jj
    
    k = floor((l-1)/(jj+1))+2;
    j = l - (k-2)*(jj+1);
    end