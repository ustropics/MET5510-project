%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: matricesBCD.m

% Description: Constructs matrices B, C, D for PV inversion and advection

% Input: 
% - ll: Total number of linear indices
% - stream2pv: Function handle to compute PV from streamfunction
% - stream2xPVadv: Function handle for zonal PV advection
% - stream2yPVadv: Function handle for meridional PV advection

% Output:
% - B: Matrix for PV inversion
% - C: Matrix for zonal PV advection
% - D: Matrix for meridional PV advection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B, C, D] = matricesBCD(ll, stream2pv, stream2xPVadv, stream2yPVadv)

%% Initialize matrices
B = zeros(ll, ll);
C = zeros(ll, ll);
D = zeros(ll, ll);

%% Construct matrices
for l0 = 1:ll
    XV = zeros(ll, 1); % reset streamfunction vector
    XV(l0) = 1; % set to 1 at index l0
    
    %% B MATRIX: PV from streamfunction
    QV = stream2pv(XV);
    B(:, l0) = QV(:);
    
    %% C MATRIX: Zonal PV advection
    xQVadv = stream2xPVadv(QV);
    C(:, l0) = xQVadv(:);
    
    %% D MATRIX: Meridional PV advection
    yQVadv = stream2yPVadv(XV);
    D(:, l0) = yQVadv(:);
end
end