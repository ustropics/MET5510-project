%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: solve_eigen.m

% Description: Solves the eigenvalue problem for the quasi-geostrophic model,
% computing eigenvectors and eigenvalues of the matrix A = B^(-1) * (C + D),
% and sorts them by descending real part and ascending imaginary part.

% Input:
% - B: PV inversion matrix (ll x ll)
% - C: Zonal PV advection matrix (ll x ll)
% - D: Meridional PV advection matrix (ll x ll)
% - cplx: Imaginary unit

% Output:
% - eigVec: Eigenvectors (unsorted)
% - eigVal: Eigenvalues (unsorted)
% - eigVec2: Eigenvectors sorted by descending real part
% - eigVal2: Eigenvalues sorted by descending real part
% - eigVec3: Eigenvectors sorted by ascending imaginary part
% - eigVal3: Eigenvalues sorted by ascending imaginary part

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eigVec, eigVal, eigVec2, eigVal2, eigVec3, eigVal3] = eigen(B, C, D, cplx)
%% Solve eigenvalue problem
A = B^(-1) * (C + D);
[eigVec, eigValm] = eig(A); % eigVec = eigenvectors, eigValm = eigenvalue matrix
eigVal = diag(eigValm);

% Sort by descending real part
[~, sortdx] = sort(real(eigVal), 'descend');
eigVal2 = eigVal(sortdx);
eigVec2 = eigVec(:, sortdx);

% Sort by ascending real part of imaginary eigenvalue
[~, sortdx] = sort(real(cplx * eigVal), 'ascend');
eigVal3 = eigVal(sortdx);
eigVec3 = eigVec(:, sortdx);
end