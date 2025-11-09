%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: combined_ubar_with_evec_amp.m

% DESCRIPTION: Plots background zonal wind Ubar (filled contour) in the
%              latitude-height plane with eigenvector amplitude overlaid
%              as black contour lines. The figure is saved as a PNG
%              including wavenumber, mode number, growth rate and frequency
%              in title and filename.

% INPUT:
%   yy          – latitude vector (degrees) , size (jj+1)
%   zz          – height vector   (km)      , size (kk+1)
%   Ubar        – background zonal wind (jj+1 x kk+1)  [m/s]
%   eVec_amp    – eigenvector amplitude (jj+1 x kk+1)
%   m0          – zonal wavenumber
%   n_mode      – eigenmode index
%   growth_rate – real part of eigenvalue
%   omega       – imaginary part of eigenvalue
%   fig_path    – directory where the figure is saved

% OUTPUT:
%   Saves:  fig_path/ubar_with_evec_amp_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combined_ubar_with_evec_amp(yy, zz, Ubar, eVec_amp, m0, n_mode, ...
                                     growth_rate, omega, fig_path)
    
    %% --------------------------------------------------------------------
    %% 1. Create figure
    %% --------------------------------------------------------------------
    fig = figure('units','inch','position',[4 2 16 12],'Visible','off');
    hold on;

    %% --------------------------------------------------------------------
    %% 2. Filled contour of Ubar (background)
    %% --------------------------------------------------------------------
    contourf(yy, zz, Ubar', 'LineStyle','none');
    cb1 = colorbar;
    colormap(flipud(cmap_spectral(256)));
    cb1.Label.String = 'Ubar (m s^{-1})';
    cb1.Label.FontSize = 20;

    %% --------------------------------------------------------------------
    %% 3. Overlay eigenvector amplitude as black contour lines
    %% --------------------------------------------------------------------
    ampRange = max(abs(eVec_amp(:)));
    levels   = linspace(-ampRange, ampRange, 13);   % symmetric, 13 lines
    [C,h] = contour(yy, zz, eVec_amp', levels, ...
                    'LineColor','k','LineWidth',4);
    clabel(C, h, 'FontSize',14, 'Color','k', 'LabelSpacing',300);

    %% --------------------------------------------------------------------
    %% 4. Axes, title, fonts
    %% --------------------------------------------------------------------
    xlabel('Latitude (degrees)','FontSize',20);
    ylabel('Height (km)','FontSize',20);

    title_str = sprintf(['Ubar (filled) & Eigenvector Amplitude (contours)\n' ...
                         'zonal wave # = %d, eMode # = %d, ' ...
                         'real eVal = %.3g, imag eVal = %.3g'], ...
                         m0, n_mode, growth_rate, omega);
    title(title_str,'FontSize',20);

    set(findall(fig,'-property','FontSize'),'FontSize',20);

    %% --------------------------------------------------------------------
    %% 5. Save & close
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ...
        sprintf('ubar_with_evec_amp_eMode-%d_m0-%d.png', n_mode, m0));
    fprintf('\nSaving Ubar + evec_amp plot to: %s\n', outFile);
    saveas(fig, outFile);
    close(fig);
end