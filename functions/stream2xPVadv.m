function xQVadv = stream2xPVadv(QV)
% Compute zonal advection of PV: -Ubar * ik * Q, where k is wavenumber
global jj kk ll Lx Ubar cplx m0

xQVadv = zeros(ll,1);

for k = 1:kk + 1
    for j = 2:jj
        l = jk2l(j,k);
        xQVadv(l) = -Ubar(j,k) * cplx * (2*pi*m0/Lx) * QV(l);
    end
end
end