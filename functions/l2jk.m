%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: l2jk.m

% Description: Converts a single linear index (l) to 2D grid indices (j, k) 
% representing latitude and height for use in the quasi-geostrophic model.

% Input:
% - l: Linear index

% Output:
% - j: Latitude index
% - k: Height index

% Math/functions: [j, k] = [l+1 - (k-1)*(jj-1), floor((l-1)/(jj-1))+1]

% - Variables:
%   - j is the number of latitude points
%   - l is the linear index adjusted for interior points 
%   - (1 to ll, where ll = (jj-1)*(kk+1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [j,k] = l2jk(l)
    global jj
    k = floor((l-1)/(jj-1))+1;
    j = l+1 - (k-1)*(jj-1);
    end