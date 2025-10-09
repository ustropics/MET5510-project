function [j,k] = l2jk(l)
global jj
k = floor((l-1)/(jj-1))+1;
j = l+1 - (k-1)*(jj-1);
end