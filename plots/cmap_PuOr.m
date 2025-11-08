%DOCUMENT filename="cmap_PuOr.m">
% PUOR - Custom colormap based on Matplotlib's PuOr
% Usage: colormap(PuOr(256))   or   cmap = PuOr(256);
%
% Interpolates between 10 sampled colors from matplotlib.cm.PuOr
% Colors sourced from Python: cmap(np.linspace(0,1,10))
%
% Hex codes for reference:
% '#7f3b08', '#b75c07', '#e68d23', '#fdc57f', '#fbead2',
% '#e5e7f0', '#bfbbda', '#8a7eb3', '#582e8c', '#2d004b'

function cmap = PuOr(n)
    if nargin < 1, n = 256; end

    % Node colors in RGB (0-1). Alpha channel is ignored in MATLAB colormaps.
    nodes = [
        0.49803922  0.23137255  0.03137255;   % #7f3b08
        0.71926182  0.36124567  0.02891196;   % #b75c07
        0.90073049  0.55132641  0.13917724;   % #e68d23
        0.99346405  0.77385621  0.49673203;   % #fdc57f
        0.98423683  0.91733948  0.82368320;   % #fbead2
        0.89950019  0.90396002  0.94186851;   % #e5e7f0
        0.74771242  0.73202614  0.85620915;   % #bfbbda
        0.54040754  0.49404075  0.70372933;   % #8a7eb3
        0.34632834  0.18216071  0.54717416;   % #582e8c
        0.17647059  0.00000000  0.29411765    % #2d004b
    ];

    % Node positions evenly spaced from 0 to 1
    pos = linspace(0, 1, size(nodes,1));

    % Interpolate to desired number of colors
    cmap = interp1(pos, nodes, linspace(0, 1, n), 'pchip');  % or use 'linear'
end