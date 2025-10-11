function l = jk2lw(j,k)
% Convert (j,k) indices to linear index for vertical motion
global jj

l = j + (k-2)*(jj+1);
end