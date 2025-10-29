%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_hovmoller.m

% DESCRIPTION: Plots the Hovmoller diagram and saves it as an image, including hlevel in title and filename.

% INPUT:
% - xx: Longitude coordinates (degrees)
% - time: Time coordinates (days)
% - gpt_h_hovmoler: Hovmoller data
% - model: Model name for filename
% - m0: Wavenumber for filename
% - n_mode: Mode number for filename
% - fig_path: Directory path for saving figure
% - hlevel: Height level (e.g., 1 or 51) for title and filename

% OUTPUT:
% - Saves plot to 'fig_path/model_hovmoller_nmode-n_mode_m0-m0_hlevel-hlevel.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_hovmoller(xx, time, gpt_h_hovmoler, model, m0, n_mode, fig_path, hlevel)

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, time, gpt_h_hovmoler', 'linestyle', 'none');
    xlabel('Longitude')
    ylabel('Time (days)')
    title(['Hovmoller Diagram (hlevel = ', num2str(hlevel), ')']);
    colorbar

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
    
    %% Save figure
    outFile = fullfile(fig_path, [model, '_hovmoller_', '_nmode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '_hlevel-', num2str(hlevel), '.png']);
    saveas(gcf, outFile);
    close(gcf);

end