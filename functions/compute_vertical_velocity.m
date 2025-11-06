function [w, wfield] = compute_vertical_velocity(XV, Ubar, params)
    % Extract
    jj = params.jj; kk = params.kk; dy = params.dy; dz = params.dz;
    m0 = params.m0; Lx = params.Lx; f0 = params.f0; beta = params.beta;
    NN2 = params.NN2; cplx = params.cplx;
    ii = params.ii; dx = params.dx;

    LW = (jj+1)*(kk-1);
    G  = zeros(LW,LW);
    for l0 = 1:LW
        wtest = zeros(LW,1); wtest(l0) = 1;
        G(:,l0) = w2ellipse(wtest);
    end

    F1 = zeros(LW,1); F2 = F1; F3 = F1;

    % F3: beta
    for j = 2:jj
        for k = 2:kk
            l = idx_jk2lw(j,k);
            F3(l) = (2*pi*m0/Lx)*cplx*beta*...
                    (XV(idx_jk2l(j,k+1)) - XV(idx_jk2l(j,k-1)))/(2*dz);
        end
    end

    % F2: shear
    for j = 2:jj
        for k = 2:kk
            l = idx_jk2lw(j,k);
            dUdy = (Ubar(j+1,k) - 2*Ubar(j,k) + Ubar(j-1,k))/dy^2;
            dXdz = (XV(idx_jk2l(j,k+1)) - XV(idx_jk2l(j,k-1)))/(2*dz);
            F2(l) = -2*(2*pi*m0/Lx)*cplx * dUdy * dXdz;
        end
    end

    % F1: boundary curvature
    for j = [1 2 jj jj+1]
        for k = 2:kk
            l = idx_jk2lw(j,k);
            if j == 1 || j == jj+1
                d2psidy2 = (XV(idx_jk2l(min(j+1,jj),k)) - 2*XV(idx_jk2l(j,k)) + XV(idx_jk2l(max(j-1,2),k))) / dy^2;
            else
                d2psidy2 = (XV(idx_jk2l(j+1,k)) - 2*XV(idx_jk2l(j,k)) + XV(idx_jk2l(j-1,k))) / dy^2;
            end
            lap_psi = -(2*pi*m0/Lx)^2 * XV(idx_jk2l(j,k)) + d2psidy2;
            F1(l) = 2*(2*pi*m0/Lx)*cplx * (Ubar(j,k+1)-Ubar(j,k-1)) * lap_psi / (2*dz);
        end
    end

    w = (f0/NN2) * (G \ (F1+F2+F3));
    wfield = w2wfield(w, ii, dx);
end