function cmap = RdYlBl(n)
    % RdYlBl - Matplotlib's Red-Yellow-Blue divergent colormap
    % Usage: colormap(RdYlBl(256))
    if nargin < 1, n = 256; end

    % Official Matplotlib RdYlBl node colors (RGB, 0-1)
    % Source: https://github.com/matplotlib/matplotlib/blob/main/lib/matplotlib/_cm.py
    nodes = [
        0.647, 0.000, 0.149;   % Dark red
        0.843, 0.188, 0.153;   % Red
        0.957, 0.427, 0.263;   % Orange-red
        0.996, 0.682, 0.459;   % Light orange
        0.992, 0.878, 0.718;   % Pale yellow
        1.000, 1.000, 0.800;   % Near white (center)
        0.800, 0.922, 0.973;   % Light cyan
        0.529, 0.808, 0.922;   % Sky blue
        0.267, 0.659, 0.871;   % Medium blue
        0.129, 0.443, 0.710;   % Blue
        0.031, 0.188, 0.420    % Dark blue
    ];

    % Node positions from 0 to 1
    pos = linspace(0, 1, size(nodes,1));

    % Interpolate to n colors
    cmap = interp1(pos, nodes, linspace(0, 1, n));
end