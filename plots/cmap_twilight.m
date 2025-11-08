    % TWILIGHT_SHIFTED - Custom colormap based on Matplotlib's twilight_shifted
    % Usage: colormap(twilight_shifted(256))
    %
    % Interpolates between 10 sampled colors from matplotlib.cm.twilight_shifted
    % Colors sourced from Python: cmap(np.linspace(0,1,10))
    %
    % Hex codes for reference:
    % '#301437', '#572385', '#6066b6', '#7ca2c2', '#c8d0d6',
    % '#dacac3', '#c6896c', '#a84750', '#6b1b4d', '#2f1436'

function cmap = cmap_twilight(n)
    if nargin < 1, n = 256; end

    % Node colors in RGB (0-1), alpha ignored in MATLAB colormaps
    nodes = [
        0.18739228, 0.0771021 , 0.21618875;  % #301437
        0.33954169, 0.13704802, 0.52328798;  % #572385
        0.37471686, 0.39955665, 0.71257369;  % #6066b6
        0.48673596, 0.63416646, 0.7625578 ;  % #7ca2c2
        0.78399101, 0.81755426, 0.83936747;  % #c8d0d6
        0.8542922 , 0.79097197, 0.76583801;  % #dacac3
        0.77590791, 0.53554218, 0.42413368;  % #c6896c
        0.65832808, 0.27803211, 0.31346434;  % #a84750
        0.41796105, 0.10420645, 0.30278652;  % #6b1b4d
        0.18488036, 0.07942573, 0.21307652   % #2f1436
    ];

    % Node positions evenly spaced from 0 to 1
    pos = linspace(0, 1, size(nodes, 1));

    % Interpolate to desired number of colors
    cmap = interp1(pos, nodes, linspace(0, 1, n), 'pchip'); % or 'linear'
end