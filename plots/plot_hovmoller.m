%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_hovmoller.m

% DESCRIPTION: Plots the Hovmoller diagram and saves it as an image.

% INPUT:
% - xx: Longitude coordinates (degrees)
% - time: Time coordinates (days)
% - gpt_h_hovmoler: Hovmoller data

% OUTPUT:
% - Saves plot to 'output/figures/hovmoller.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_hovmoller(xx, time, gpt_h_hovmoler, model, m0)
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, time, gpt_h_hovmoler', 'linestyle', 'none');
    xlabel('Longitude')
    ylabel('Time (days)')
    title('Hovmoller Diagram');
    colorbar

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);
    
    % Save plot
    saveas(gcf, ['output', filesep, 'figures', filesep, 'hovmoller_', model, '_m0_', num2str(m0), '.png']);
    close(gcf);
end