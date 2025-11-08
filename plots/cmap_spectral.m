    % SPECTRAL_CUSTOM - Custom colormap based on provided spectral colors
    % Usage: colormap(spectral_custom(256))
    %
    % Interpolates between 10 sampled spectral colors
    % Hex codes for reference:
    % '#9e0142', '#d8434e', '#f67a49', '#fdbf6f', '#feeda1',
    % '#f1f9a9', '#bfe5a0', '#74c7a5', '#378ebb', '#5e4fa2'

function cmap = cmap_spectral(n)
    if nargin < 1, n = 256; end

    % Node colors in RGB (0-1)
    nodes = [
        0.61960784, 0.00392157, 0.25882353;  % #9e0142
        0.84721261, 0.26120723, 0.30519031;  % #d8434e
        0.96378316, 0.47743176, 0.28581315;  % #f67a49
        0.99346405, 0.74771242, 0.43529412;  % #fdbf6f
        0.99777009, 0.93087274, 0.63306421;  % #feeda1
        0.94425221, 0.97770088, 0.66205306;  % #f1f9a9
        0.74771242, 0.89803922, 0.62745098;  % #bfe5a0
        0.45305652, 0.78154556, 0.64628989;  % #74c7a5
        0.21607074, 0.55563245, 0.73194925;  % #378ebb
        0.36862745, 0.30980392, 0.63529412   % #5e4fa2
    ];

    % Node positions evenly spaced from 0 to 1
    pos = linspace(0, 1, size(nodes, 1));

    % Interpolate to desired number of colors
    % Using 'pchip' for smooth, shape-preserving interpolation (like matplotlib)
    cmap = interp1(pos, nodes, linspace(0, 1, n), 'pchip');
end