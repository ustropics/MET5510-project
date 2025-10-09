% process_variables.m
% Process variables for Rossby wave analysis
global jj kk ll BPVy NN2 f0 dy dz m0 Lx Ubar f0 beta cplx HH gg ii dx eigVec3 eigVal3 B

function [eVec_amp, gpt_h, omega, phase_speed, growth_rate, eFolding, Amp, QV, XVy, XVx, XVz, temp, ug, vg, pvfield, gpt_h_hovmoler, pv, time] = process_variables
    % Load required data (assumes Rossby_wave_2.mat is already loaded)
    global eigVec3 eigVal3 B
    
    XV = zeros(ll,1); % initialize streamfunction vector
    n_mode = 7; % this goes west fast, 2 is unstable
    XV(:) = eigVec3(:,n_mode);
    omega = imag(eigVal3(n_mode));
    phase_speed = -omega/(2*pi*m0/Lx);

    % from Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf
    % slide 6
    growth_rate = real(eigVal3(n_mode)); % real part of eigenvalue
    eFolding = (1/growth_rate)/86400; % in units of days

    % Amplitude of eigenvector
    eVec_amp = zeros(jj+1, kk+1); % slide 7-8
    for l = 1:ll
        [j,k] = l2jk(l);
        eVec_amp(j,k) = XV(l).*conj(XV(l));
    end

    % set geopotential field
    gpt_h = XV2field(XV,ii,dx) * f0/gg;

    latmax = (jj/2 + 1)/2;
    levelmax = 1;

    % finds the maximum value index converting 3d to 1d
    [valuemax, indexmax] = max(gpt_h(:));

    % normalize the value
    XV = (10/valuemax) * XV;

    Amp = sqrt(sum(XV.*XV)/ll);
    QV = B * XV;
    XVy = XV2XVy(XV); % d/dy at N and S boundary not calculated
    XVx = XV2XVx(XV);
    XVz = XV2XVz(XV);

    % from Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf
    gpt_h = XV2field(XV,ii,dx)*f0/gg; % 3d geopotential height from slide 8
    temp = (f0*HH/287) * XVz2field(XVz, ii, dx); % from slide 11
    ug = -XV2field(XVy, ii, dx);
    vg = XV2field(XVx, ii, dx);
    pvfield = XV2field(QV, ii, dx);

    % run series from day 0 to day 50, x from 0 to 360 degrees
    hlat = floor(jj/4 + 1); % lat for hovmoller diagram
    hlevel = 1; % vertical level for hovmoller diagram

    gpt_h_hovmoler = XV2streamxtime(XV, ii, dx, omega, hlat, hlevel) * f0/gg;
    pv = XV2field(QV,ii,dx);

    time = 0:1:50;
end