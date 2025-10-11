function xtime = XV2streamxtime(XV,ii,dx,omega,ylat,zlev)
% Compute hovmoller diagram (streamfunction over time)
% From Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf
% Page 6 (diagnostic analysis)
global cplx m0 Lx

xtime = zeros(ii+1, 51);

l = jk2l(ylat, zlev);
for day = 1:51
    sec = (day-1)*86400;
    for n = 1:ii+1
        xlon=(n-1)*dx;
        xtime(n,day) = real(XV(l)*exp(cplx*(2*pi*m0*xlon/Lx+omega*sec)));
    end
end
end