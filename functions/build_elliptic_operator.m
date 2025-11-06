%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: build_elliptic_operator.m

% DESCRIPTION: Constructs the elliptic operator matrix G (LW×LW) used to
% solve for the vertical velocity w.  Applies the operator w2ellipse to
% every unit vector.

% INPUT:  (none – uses globals)
% OUTPUT:
%   G – LW×LW sparse matrix

function G = build_elliptic_operator()
    global jj kk LW

    G = zeros(LW,LW);

    for l0 = 1:LW
        w   = zeros(LW,1);
        w(l0) = 1;
        EW  = w2ellipse(w);                % elliptic operator applied to unit vector
        G(:,l0) = EW(:);
    end
end