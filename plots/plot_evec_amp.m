%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_evec_amp.m

% DESCRIPTION: Plots the eigenvector amplitude contour and saves it as an image.

% INPUT:
% - yy: Latitude coordinates (degrees)
% - zz: Height coordinates (km)
% - eVec_amp: Eigenvector amplitude matrix
% - m0: Zonal wave number
% - n_mode: Mode number
% - growth_rate: Growth rate
% - omega: Imaginary eigenvalue

% OUTPUT:
% - Saves plot to 'output/figures/evec_amp.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_evec_amp(yy, zz, eVec_amp, m0, n_mode, growth_rate, omega, model)

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(yy, zz, eVec_amp', 'linestyle', 'none');
    xlabel('Latitude')
    ylabel('Height (km)')
    set(gca, 'fontsize', 22, 'color', 'w')
    title_str = ['Zonal wave # = ', num2str(m0), ', eMode # =', num2str(n_mode), ...
        ', real eVal = ', num2str(growth_rate), ', imag eVal = ', num2str(omega)];
    title(title_str);
    colorbar
    
    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);

    %% Save plot
    saveas(gcf, ['output', filesep, 'figures', filesep, 'evec_amp_', model, '_m0_', num2str(m0), '.png']);
    close(gcf);

end