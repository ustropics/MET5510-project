%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: stream2pv.m

% Description: Converts streamfunction (XV) to potential vorticity (QV) 
% for the quasi-geostrophic modeapplying boundary conditions and 
% interior PV calculation.

% Math/functions: Q = ∇²ψ + (f₀²/N²) ∂²ψ/∂z² + β, where 
% ψ is streamfunction
% f₀ is Coriolis parameter
% N² is Brunt-Vaisala frequency
% β is meridional PV gradien

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% functions/stream2pv.m
function QV = stream2pv(XV)
global jj kk ll NN2 m0 f0 dy dz Lx

QV = zeros(ll,1);

% bottom boundary (k = 1)
k = 1; 
for j = 2:jj
    l = jk2l(j,k); % from j value we can get eigen value
    lup = jk2l(j, k+1); % get one level from above
    ldn = l; % we've reached bottom boundary

    QV(l) = ( XV(lup) - XV(ldn) )/dz;
end

% top boundary (k=kk+1)
k = kk + 1;
for j = 2:jj
    l = jk2l(j, k);
    ldn = jk2l(j, k-1);
    lup = l; % we've reached top boundary

    QV(l) = ( XV(lup) - XV(ldn) ) / dz;
end

% interior points (quasi-geostrophic PV)
for k = 2:kk

    for j = 2:jj

        l = jk2l(j, k); % current index
        lup = jk2l(j, k+1); % above
        ldn = jk2l(j, k-1); % below
        lnh = jk2l(j+1,k); % north (+j)
        lsh = jk2l(j-1,k); % south (-j)

        % meridional boundary condition
        if (j == 2)
            XVsh=0; % southern boundary
        else
            XVsh = XV(lsh);
        end

        if (j == jj)
            XVnh = 0; % northern boundary

        else
            XVnh = XV(lnh);

        end

        % PV = - (k^2 + laplacian_psi) + (f0^2 / N^2) * d^2 psi / dz^2
        % Where k = 2*pi*m0/Lx  (zonal wavenumber)
        QV(l) = -((2*pi*m0/Lx)^2)*XV(l) ...
            + ( XVsh - 2*XV(l) + XVnh) / dy/dy ... 
            + (f0/dz)^2 * (XV(lup) - 2*XV(l) + XV(ldn)) / NN2;

    end

end
end