% functions/XV2XVz.m
function XVz = XV2XVz(XV)
global jj kk ll dz
XVz = zeros(ll,1);

for k = 2:kk
    for j = 2:jj

        l = jk2l(j,k);
        lup = jk2l(j,k+1);
        ldw = jk2l(j, k-1);

        XVz(l) = (XV(lup)-XV(ldw))/2/dz;
    end
end

k = 1;
for j = 2:jj
    l = jk2l(j,k);
    lup = jk2l(j, k+1);
    ldw = l;
    XVz(l) = (XV(lup) - XV(ldw))/dz;
end

k = kk+1;

for j = 2:jj
    l = jk2l(j,k);
    ldw = jk2l(j, k-1);
    lup=l;
    XVz(l) = (XV(lup) - XV(ldw))/dz;
end
end