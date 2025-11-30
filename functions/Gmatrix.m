%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: Gmatrix.m
%
% DESCRIPTION: Constructs the sparse elliptic operator matrix G (LW×LW) for 
% solving the QG omega equation G w = (f₀/N²) (F1+F2+F3). Each column of G is 
% the result of applying the diagnostic elliptic operator (as in w2ellipse) to 
% a unit vector in w-space.
%
% INPUT: None (uses globals)
%
% OUTPUT:
% - G: Dense matrix (LW × LW) representing discretized elliptic operator
%
% MATH/FUNCTIONS: 
% - G w = -k²w + ∂²w/∂y² + (f₀²/N²) ∂²w/∂z²
%
% VARIABLES:
% - LW: Number of interior grid points = (jj+1)(kk-1)
% - w2ellipse: Function applying elliptic operator to input vector
% - jk2lw: Index mapping from (j,k) → linear index in w_vec
% - Matrix built column-by-column via unit vector application
% - Result is dense in memory but can be sparsified externally if needed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = Gmatrix()
    global jj kk LW

    % Initialize full dense matrix (LW x LW)
    G = zeros(LW, LW);

    %% Loop over each column: apply operator to unit vector e_l0
    for l0 = 1:LW

        w = zeros(LW, 1); % create unit vector: w = [0, 0, ..., 1, ..., 0]
        w(l0) = 1;
        
        % Apply elliptic diagnostic operator (same as in w2ellipse)
        EW = w2ellipse(w);
        G(:, l0) = EW(:); 
    end
end