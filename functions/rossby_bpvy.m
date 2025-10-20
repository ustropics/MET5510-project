%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_bpvy.m

% Description: Initializes the potential vorticity gradient field (BPVy) for 
% the Rossby wave model in the quasi-geostrophic framework, set to beta for interior points.

% Input:
% - jj: Number of latitude grid points
% - kk: Number of height grid points
% - beta: Planetary vorticity gradient (m^-1 s^-1)

% Output:
% - BPVy: Potential vorticity gradient field (jj+1 x kk+1 array, m^-1 s^-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BPVy = rossby_bpvy(jj, kk, beta)
BPVy = zeros(jj + 1, kk + 1);
for k = 2 : kk
    for j = 2:jj
        BPVy(j,k) = beta; % assign constant beta to the PV gradient
    end
end
end