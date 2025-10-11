function l = jk2l(j,k)
% Convert j,k indices to linear index l
global jj

l = j - 1 + (k-1) * (jj-1);
end