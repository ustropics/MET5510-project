clc;

%% Compute values

sec_to_day      = 86400;                                    % day                                    
growth_rate     = real(eigVal2(n_mode));                    % s^-1
efolding        = 1/growth_rate/sec_to_day;                 % days                     
phase_speed     = -imag(eigVal2(n_mode)) / (2*pi*m0/Lx);    % m s^-1
phase_speed_day = Lx / phase_speed / sec_to_day;            % days

period          = 2*pi/abs(omega)/sec_to_day;               % days
k               = 2*pi*m0/Lx;                               % rad m^-1
group_speed     = phase_speed + (omega/(2*pi*m0/Lx))*0;     % for non-dispersive wave = c_p
rossby_radius   = params.Lr;                                % km
rossby_number   = max(abs(ug(:))) / (f0*rossby_radius);     % dimensionless
max_T           = max(abs(temp(:)));                        % K
max_w           = max(abs(wfield(:)));                      % cm s^-1

max_bottom = max(eVec_amp(:,1));
max_top    = max(eVec_amp(:,51));
fprintf('bottom = %.4e , top = %.4e\n', max_bottom, max_top);

%% Print values and units
fprintf('eMode # = %d (zonal wave # = %d): phase speed = %.3f m/s\n', ...
        n_mode, m0, phase_speed);

fprintf('Phase speed in days : %.3f m/days\n', phase_speed_day)

% Calculate the Rossby number and print it
fprintf('Rossby number: %.2e \n', rossby_number);

% Calculate the e-folding time and print it
fprintf('E-folding time: %.3f days\n', efolding);

% Calculate the group speed and print it
fprintf('Group speed: %.3f m/s\n', group_speed);

% Calculate and print the maximum values of the fields
fprintf('Maximum temperature: %.3f K\n', max_T);
fprintf('Maximum wind speed: %.3e cm/s\n', max_w);

% Additional calculations or analyses can be performed here if needed
% For example, calculating the energy associated with the wave
energy = 0.5 * (max_gph^2 + max_T^2 + max_w^2);  % Example energy calculation
% Print the calculated energy associated with the wave
fprintf('Energy associated with the wave: %.3f J\n', energy);

% Calculate the total energy density associated with the wave
total_energy_density = energy / (Lx * rossby_radius);  % J/km^2
fprintf('Total energy density: %.3e J/km^2\n', total_energy_density);

% Calculate the potential energy associated with the wave
potential_energy = 0.5 * max_gph^2;  % Example potential energy calculation
fprintf('Potential energy associated with the wave: %.0f J\n', potential_energy);

% Calculate the kinetic energy associated with the wind field
kinetic_energy = 0.5 * max_w^2;  % Example kinetic energy calculation
fprintf('Kinetic energy associated with the wind field: %.3e J\n', kinetic_energy);

% Calculate the ratio of kinetic to potential energy
kinetic_to_potential_ratio = kinetic_energy / potential_energy;  
fprintf('Kinetic to potential energy ratio: %.3e \n', kinetic_to_potential_ratio);

% Calculate the wave energy flux
wave_energy_flux = group_speed * energy;  % J/s
fprintf('Wave energy flux: %.3f J/s\n', wave_energy_flux);

% Calculate the wave number based on the group speed and phase speed
wave_number = group_speed / phase_speed;  % rad/m
fprintf('Wave number: % rad/m\n', wave_number);

% Calculate the wave energy density
wave_energy_density = wave_energy_flux / (Lx * sec_to_day);  % J/km^2
fprintf('Wave energy density: %.4e J/km^2\n', wave_energy_density);

% Calculate the dispersion relation for the wave
dispersion_relation = (omega.^2 - (k.^2 * phase_speed^2)); 
fprintf('Dispersion relation: %.4e \n', dispersion_relation);

