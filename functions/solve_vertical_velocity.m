%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: solve_vertical_velocity.m

% DESCRIPTION: Solves the elliptic equation G w = (f0/NN2)*(F1+F2+F3)
% for the interior vertical velocity vector w_vec.

% INPUT:
%   G  – elliptic operator matrix (LW×LW)
%   F1,F2,F3 – forcing vectors (LW×1)

% OUTPUT:
%   w_vec – interior vertical velocity (LW×1)

function w_vec = solve_vertical_velocity(G,F1,F2,F3)
    global f0 NN2
    w_vec = (f0/NN2) * (G \ (F1 + F2 + F3));
end