function [j,k] = idx_lw2jk(l)
    global jj
    j = mod(l-1, jj+1) + 1;
    k = floor((l-1)/(jj+1)) + 2;
end