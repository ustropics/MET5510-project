function yQVadv = stream2yPVadv(XV)
% Compute meridional advection of PV: -BPVy * ik * Q, where k is wavenumber
global jj kk ll BPVy Lx cplx m0

yQVadv = zeros(ll,1);

for k = 1:kk + 1
    for j = 2:jj
        l = jk2l(j,k);
        yQVadv(l) = -BPVy(j,k) * cplx * (2*pi*m0/Lx) * XV(l);
    end
end
end