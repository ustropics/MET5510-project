%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: w2vec.m

% DESCRIPTION: Solves the elliptic equation G w = (f0/N²) * (F1 + F2 + F3) 
% for the interior vertical velocity vector in the quasi-geostrophic (QG) model. 
% This function computes the vertical motion induced by ageostrophic forcing 
% terms, critical for diagnosing secondary circulation and frontogenetic processes 
% within the QG framework.

% INPUT:
% - G:  Elliptic operator matrix (LW x LW), representing discretized form of 
%       ∇·(∇w) or similar, where LW is the number of interior (j,k) points
% - F1, F2, F3: Forcing vectors (LW x 1), representing components of the 
%       Q-vector divergence or other ageostrophic forcing terms

% OUTPUT:
% - w_vec: 1D array representing the interior vertical velocity (m/s), 
%           with length LW (only interior points, excluding boundaries)

% MATH/FUNCTIONS: 
% - w = (f0/N²) * G⁻¹ * (F1 + F2 + F3)

% VARIABLES:
% - Interior points only: w(j,k) solved via direct matrix inversion
% - f0: Coriolis parameter at reference latitude (s⁻¹)
% - NN2: Reference Brunt-Väisälä frequency squared, N² (s⁻²)
% - G: Sparse elliptic operator matrix (LW x LW)
% - F1, F2, F3: Forcing contributions from geostrophic deformation, 
%   differential vorticity advection, and thermal advection, respectively
% - The solution uses MATLAB's backslash operator (\) for efficient 
%   solution of the linear system G w = rhs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w_vec = w2vec(G,F1,F2,F3)
    global f0 NN2

    %% Solve for interior vertical velocity vector
    w_vec = (f0/NN2) * (G \ (F1 + F2 + F3));
    
end