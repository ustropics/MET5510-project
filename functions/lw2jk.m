%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: lw2jk.m

% Description: Converts a single linear index (l) to 2D grid indices (j, k) 
% representing latitude and height


% Math/functions: [j, k] = [mod((l-offset)/w, jj-1)+1,
% floor(((l-offset)/w)/(jj-1))+1], where 
% w is a weighting factor, offset adjusts for grid boundaries
% jj is the number of latitude points 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [j,k] = lw2jk(l)
global jj

k = floor((l-1)/(jj+1))+2;
j = l - (k-2)*(jj+1);
end