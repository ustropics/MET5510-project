function [j,k] = idx_l2jk(l)
    global jj
    k = floor((l-1)/(jj-1)) + 1;
    j = mod(l-1, (jj-1)) + 1;
end