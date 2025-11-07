%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_hovmoller.m

% DESCRIPTION: Plots the Hovmoller diagram and saves it as an image,
% including hlevel in title and filename.

% INPUT:
% - xx: Longitude coordinates (degrees)
% - time: Time coordinates (days)
% - gpt_h_hovmoler: Hovmoller data
% - model: Model name for filename
% - m0: Wavenumber for filename
% - n_mode: Mode number for filename
% - fig_path: Directory path for saving figure
% - hlevel: Height level (e.g., 1 or 51) for title and filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_gph_hovmoller(xx, time, gpt_h_hovmoler, m0, n_mode, fig_path, hlevel)

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, time, gpt_h_hovmoler', 'linestyle', 'none');
    xlabel('Longitude')
    ylabel('Time (days)')

    title_str = ['Hovmoller Diagram (hlevel = ', num2str(hlevel), ...
    ', zonal wave # = ', num2str(m0), ...
    ', eMode # = ', num2str(n_mode), ')'];

    title(title_str);
    colorbar

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
    
    %% Save figure
    outFile = fullfile(fig_path, ['gph-hovmoller_','hlevel-', num2str(hlevel), ... 
        '_eMode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);
    saveas(gcf, outFile);
    close(gcf);

end