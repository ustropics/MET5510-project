%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: bpvy.m

% Description: BPVY computes the background potential vorticity 
% gradient (dPV/dy)for Rossby, Eady, or HWME models.

% Inputs:
% - jj: Number of latitude grid points
% - kk: Number of height grid points
% - beta: Planetary vorticity gradient (m^-1 s^-1)
% - Ubar: Mean zonal wind field (jj+1 x kk+1 array, m/s)
% - dy: Meridional grid spacing (m)
% - f0: Coriolis parameter (s^-1)
% - NN2: Brunt-Vaisala frequency squared (s^-2)
% - dz: Vertical grid spacing (m)

%  Output:
% - BPVy - background PV gradient field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BPVy = bpvy(jj, kk, beta, Ubar, dy, f0, NN2, dz)
    
    BPVy = zeros(jj + 1, kk + 1);
    model = 'HWME';

    % Compute for Rossby model
    if strcmpi(model, 'Rossby')
        for k = 2:kk
            for j = 2:jj
                BPVy(j, k) = beta;
            end
        end

    % Compute for Eady or HWME models
    elseif strcmpi(model, 'Eady') || strcmpi(model, 'HWME')
        for j = 2:jj
            for k = 2:kk
                BPVy(j, k) = beta ...
                    - (Ubar(j + 1, k) - 2 * Ubar(j, k) + Ubar(j - 1, k)) / dy^2 ...
                    - (f0^2 / NN2) * (Ubar(j, k + 1) - 2 * Ubar(j, k) + Ubar(j, k - 1)) / dz^2;
            end

            % Surface and upper boundary conditions
            k = 1;
            BPVy(j, k) = -(Ubar(j, k + 1) - Ubar(j, k)) / dz;
            k = kk + 1;
            BPVy(j, k) = -(Ubar(j, k) - Ubar(j, k - 1)) / dz;
        end

    else
        error('Unknown model type: %s. Use "Rossby", "Eady", or "HWME".', model);
    end

end
