%DOCUMENT filename="cmap_PRGn.m">
% PRGN - Custom colormap based on Matplotlib's PRGn
% Usage: colormap(PRGn(256))   or   cmap = PRGn(256);
%
% Interpolates between 10 sampled colors from matplotlib.cm.PRGn
% Colors sourced from Python: cmap(np.linspace(0,1,10))
%
% Hex codes for reference:
% '#40004b', '#793187', '#a17ab2', '#ceb5d7', '#eee3ee',
% '#e6f3e3', '#b7e2b1', '#69b76d', '#217d3b', '#00441b'

function cmap = cmap_PRGn(n)
    if nargin < 1, n = 256; end

    % Node colors in RGB (0-1).  The 4-th column (alpha) is ignored in MATLAB colormaps.
    nodes = [
        0.25098039  0.00000000  0.29411765;   % #40004b
        0.47620146  0.19161861  0.52910419;   % #793187
        0.63152634  0.47996924  0.69826990;   % #a17ab2
        0.80915033  0.70849673  0.84444444;   % #ceb5d7
        0.93294887  0.89058055  0.93517878;   % #eee3ee
        0.90173010  0.95301807  0.88835063;   % #e6f3e3
        0.71764706  0.88627451  0.69411765;   % #b7e2b1
        0.41138024  0.71695502  0.42883506;   % #69b76d
        0.13010381  0.49134948  0.23183391;   % #217d3b
        0.00000000  0.26666667  0.10588235    % #00441b
    ];

    % Node positions evenly spaced from 0 to 1
    pos = linspace(0, 1, size(nodes,1));

    % Interpolate to the requested number of colors
    cmap = interp1(pos, nodes, linspace(0, 1, n), 'pchip');   % change to 'linear' if desired
end