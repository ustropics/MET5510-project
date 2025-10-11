function [j,k] = lw2jk(l)
% Convert linear index to (j,k) indices for vertical motion
global jj

k = floor((l-1)/(jj+1))+2;
j = l - (k-2)*(jj+1);
end