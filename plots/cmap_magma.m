%DOCUMENT filename="cmap_magma.m">
% MAGMA - Custom colormap based on Matplotlib's magma
% Usage: colormap(magma(256))   or   cmap = magma(256);
%
% Interpolates between 10 sampled colors from matplotlib.cm.magma
% Colors sourced from Python: cmap(np.linspace(0,1,10))
%
% Hex codes for reference:
% '#7f3b08', '#b75c07', '#e68d23', '#fdc57f', '#fbead2',
% '#e5e7f0', '#bfbbda', '#8a7eb3', '#582e8c', '#2d004b'

function cmap = cmap_magma(n)
    if nargin < 1, n = 256; end

    % Node colors in RGB (0-1). Alpha is ignored in MATLAB colormaps.
    nodes = [
        1.46200e-03  4.66000e-04  1.38660e-02;   % #000004
        9.29490e-02  5.99040e-02  2.39164e-01;   % #180f3d
        2.65447e-01  6.02370e-02  4.61840e-01;   % #440f76
        4.45163e-01  1.22724e-01  5.06901e-01;   % #721f81
        6.20005e-01  1.83840e-01  4.97524e-01;   % #9e2f7f
        8.04752e-01  2.49911e-01  4.42102e-01;   % #cd4071
        9.44006e-01  3.77643e-01  3.65136e-01;   % #f1605d
        9.92196e-01  5.87502e-01  4.06299e-01;   % #fd9668
        9.96369e-01  7.91167e-01  5.53499e-01;   % #feca8d
        9.87053e-01  9.91438e-01  7.49504e-01    % #fcfdbf
    ];

    % Node positions evenly spaced from 0 to 1
    pos = linspace(0, 1, size(nodes,1));

    % Interpolate to desired number of colors
    cmap = interp1(pos, nodes, linspace(0, 1, n), 'pchip');   % or 'linear'
end