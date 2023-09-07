% Copyright (c) 2023 Luis Romero-Ben <luis.romero.ben@upc.edu>
% 
% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

clear;clc

%% Path configuration %%

path = pwd;
ppath = erase(path,'\code');
data_path = [ppath '\data'];
hyd_data_path = [data_path '\hydraulic\'];
struc_data_path = [data_path '\structural\'];

%% Load data %%

% Hydralic data %

% CONFIGURATION: Change the leaktionary name parameters to perform interpolation over other datasets.
network_name = 'Modena';
min_leak_size = 4.5; min_leak_size = erase(num2str(min_leak_size),'.');
max_leak_size = 6.5; max_leak_size = erase(num2str(max_leak_size),'.');
uncertainty = 0;   uncertainty = erase(num2str(uncertainty),'.');
data_purpose = 'testing';

try 
load([hyd_data_path 'leaktionary_' network_name ...
      '_' min_leak_size 'to' max_leak_size ...
      '_' uncertainty 'p100_' data_purpose '.mat'])
catch
    error(['The configured leaktionary name does not exist. Please, check the available names in: ' hyd_data_path])
end

% Structural data %

load([struc_data_path 'ModenaInfo.mat']); % Structural info of the network
load([struc_data_path 'ModenaPipeDistance.mat']); % Matrix of pipe distances
load([struc_data_path 'ModenaIncidenceMatrix.mat']); % Incidence matrix
load([struc_data_path 'best_individual_20_withReservoir.mat']); % Set of actual sensors

d = graphInfo.d; A = graphInfo.A; WA = graphInfo.WA;
node_coordenades = graphInfo.nc; Elevation = graphInfo.Elevation;
G = graphInfo.G; Lengths = graphInfo.Lengths; N = double(graphInfo.N);
reservoirsID = graphInfo.reservoirsID;

%% Generate interpolation related matrices (Laplacians, sensor matrix) %%

% WLap  -> smoothing-Laplacian        | WLap is the actual Laplacian of the
%                                             network graph, considering the network
%                                       lengths as the graph weights
% WLap2 -> analytic-weights Laplacian | WLap2 is the actual Laplacian of the
%                                       network graph, considering the analytic
%                                       weights as the graph weights
% WLap3 -> original GSI Laplacian     | WLap3 is the matrix of the
%                                       quadratic term of the original GSI problem.

WDeg = diag(sum(WA));
WLap = WDeg - WA;
WDeg2 = WDeg.^-2; WDeg2(isinf(WDeg2))=0;
WLap3 = WLap*WDeg2*WLap;

% Derive the sensor matrix S for the interpolation problem
s = zeros(N,1); s(sensors) = 1;
S = eye(N).*(s*s');
    
%% Interpolation %%

time = 1:size(Leaktionary{N}.Time,1); 
tau = 1000;

%%%%%%%%%%%% ARGUMENTS TO LOOP %%%%%%%%%%%%%%
tic
leaks = [N 1:length(leaky_node_v)]; % includes the nominal case and the leaks

for l = leaks       
for t = time     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

duration_leak = 0; % only to settle the starting time instant
precision = 100;   % sensor precision - 10 = 10 cm; 100 = 1 cm

fprintf('\n####### LEAK %d | TIME INSTANT %d ####### \n',l,t);

% Some actions are different for the nominal and leak cases

if ~isempty(intersect(reservoirsID,l))

    [~,x_nom_orig,x_leak_orig] = load_leak(...
        reservoirsID(end),Leaktionary,reservoirsID(end));
    pressureRaw = x_nom_orig(t+duration_leak,:)';
    pressureRaw = fix(pressureRaw*precision)/precision;
    hN = pressureRaw(sensors);
    hnom(:,t) = hN;

    % ###################### smoothing GSI ######################

    disp('sGSI')

    [f,~] = GSI(hN,WLap,S,I,...
        reservoirsID,...
        sensors,tau);
    fnom(:,t) = f;
    x_leak_sGSI_M(l,:,t) = fnom(:,t);
    res_leak_sGSI_M(l,:,t) = f - fnom(:,t);
    x_leak_sGSIABS_M(l,:,t) = fnom(:,t);
    res_leak_sGSIABS_M(l,:,t) = f - fnom(:,t);

    % ###################### GSI ######################

    disp('GSI')

    [f,diagnosis] = GSI(hN,WLap3,S,I,...
        reservoirsID,...
        sensors,tau);
    fnom_gsi(:,t) = f;
    x_leak_GSI_M(l,:,t) = fnom_gsi(:,t);
    res_leak_GSI_M(l,:,t) = f - fnom_gsi(:,t);
    x_leak_GSIABS_M(l,:,t) = fnom_gsi(:,t);
    res_leak_GSIABS_M(l,:,t) = f - fnom_gsi(:,t);

    % ###################### AW-GSI ######################

    disp('AW-GSI');

    [G,B,I2] = compute_GandB(graphInfo,fnom(:,t));
    WLap2 = compute_LapfromGandB(fnom(:,t),G,B,2,graphInfo);

    [x_AW,optinfo_AW] = GSI(hN,WLap2,S,-I2,...
        graphInfo.reservoirsID,...
        sensors,tau);
    fnomAW(:,t) = x_AW;
    x_leak_AW_M(l,:,t) = fnomAW(:,t);
    res_leak_AW_M(l,:,t) = x_AW - fnomAW(:,t);
    x_leak_AWABS_M(l,:,t) = fnomAW(:,t);
    res_leak_AWABS_M(l,:,t) = x_AW - fnomAW(:,t);

else

    [~,x_nom_orig,x_leak_orig] = load_leak(...
        leaky_node_v(l),Leaktionary,reservoirsID(end));
    pressureRaw = x_leak_orig(t+duration_leak,:)';
    pressureRaw = fix(pressureRaw*precision)/precision;
    hN = pressureRaw(sensors);
    hN_h = hN;
    hN = hN - hnom(:,t);

    % ###################### smoothing GSI ######################

    disp('sGSI')

    [f,diagnosis] = GSI(hN_h,WLap,S,I,...
        reservoirsID,...
        sensors,tau);
    x_leak_sGSIABS_M(l,:,t) = f;
    res_leak_sGSIABS_M(l,:,t) = f - fnom(:,t);

    error_sGSIABS(t,l) = norm(f-pressureRaw);
    optinfo_leak_sGSIABS(l,t)= diagnosis.exitflag;

    [f,diagnosis] = GSI_res(hN,WLap,S,sensors);
    x_leak_sGSI_M(l,:,t) = f + fnom(:,t);
    res_leak_sGSI_M(l,:,t) = f;

    error_sGSI(t,l) = norm(x_leak_sGSI_M(l,:,t)'-pressureRaw);
    optinfo_leak_sGSI(l,t)= diagnosis.exitflag;

    % ###################### GSI ######################

    disp('GSI')

    [f,diagnosis] = GSI(hN_h,WLap3,S,I,...
        reservoirsID,...
        sensors,tau);
    x_leak_GSIABS_M(l,:,t) = f;
    res_leak_GSIABS_M(l,:,t) = f - fnom_gsi(:,t);

    error_GSIABS(t,l) = norm(f-pressureRaw);
    optinfo_leak_GSIABS(l,t)= diagnosis.exitflag;

    [f,diagnosis] = GSI_res(hN,WLap3,S,sensors);
    x_leak_GSI_M(l,:,t) = f + fnom_gsi(:,t);
    res_leak_GSI_M(l,:,t) = f;

    error_GSI(t,l) = norm(x_leak_GSI_M(l,:,t)'-pressureRaw);
    optinfo_leak_GSI(l,t)= diagnosis.exitflag;

    % ###################### AW-GSI ######################

    disp('AW-GSI');

    [G,B,I2] = compute_GandB(graphInfo,fnom(:,t));
    WLap2 = compute_LapfromGandB(fnom(:,t),G,B,2,graphInfo);

    [x_AW,diagnosis] = GSI(hN_h,WLap2,S,-I2,...
        reservoirsID,...
        sensors,tau);
    x_leak_AWABS_M(l,:,t) = x_AW;
    res_leak_AWABS_M(l,:,t) = x_AW - fnomAW(:,t);
    error_AWABS(t,l) = norm(x_AW-pressureRaw);
    optinfo_leak_AWABS(l,t)= diagnosis.exitflag;

    [res_AW,optinfo_AW] = GSI_res(hN,WLap2,S,sensors);
    x_leak_AW_M(l,:,t) = res_AW + fnomAW(:,t);
    res_leak_AW_M(l,:,t) = res_AW;
    error_AW(t,l) = norm(x_leak_AW_M(l,:,t)'-pressureRaw);
    optinfo_leak_AW(l,t)= optinfo_AW.exitflag;

end

x_leak_REAL_M(l,:,t) = pressureRaw;

end
end

%% Save interpolated data %%

int_data_path = [data_path '/interpolated_data/'];

if ~exist(int_data_path, 'dir')
    mkdir(int_data_path)
end

save([int_data_path 'interpolation_Modena_with_RESID_' 'tau_' num2str(tau) '_fnom_' uncertainty 'noise_' min_leak_size 'to' max_leak_size '_' data_purpose '.mat'], ...
                                                                 'x_leak_sGSI_M','optinfo_leak_sGSI','error_sGSI','res_leak_sGSI_M',...
                                                                 'x_leak_GSI_M','optinfo_leak_GSI','error_GSI','res_leak_GSI_M',...
                                                                 'x_leak_AW_M','optinfo_leak_AW','error_AW','res_leak_AW_M', ...
                                                                 'x_leak_sGSIABS_M','optinfo_leak_sGSIABS','error_sGSIABS','res_leak_sGSIABS_M',...
                                                                 'x_leak_GSIABS_M','optinfo_leak_GSIABS','error_GSIABS','res_leak_GSIABS_M',...
                                                                 'x_leak_AWABS_M','optinfo_leak_AWABS','error_AWABS','res_leak_AWABS_M',...
                                                                 'x_leak_REAL_M','res_leak_REAL_M','fnom','fnomAW');
