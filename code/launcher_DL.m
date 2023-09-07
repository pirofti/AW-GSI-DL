% Copyright (c) 2023 Luis Romero-Ben <luis.romero.ben@upc.edu>
% Copyright (c) 2018, 2019, 2023 Paul Irofti <paul@irofti.net>
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

rng('default')

clc; close all; clear all

%% Path configuration %%

path = pwd;
ppath = erase(path,'\code');
data_path = [ppath '\data'];
hyd_data_path = [data_path '\hydraulic\'];
struc_data_path = [data_path '\structural\'];
vs_data_path = [data_path '\vsensors\'];

%% Load general data %%

% Virtual sensors %

nvs = 10;

try 
load([vs_data_path 'best_individual_' num2str(nvs + 20) ... % modify this 20 by the number of real sensors
      '_withReservoir_VIRTUAL' num2str(nvs) 'SENSORS.mat'])
virtual_sensors = sensors;
catch
    error(['The configured VS file name does not exist. Please, check the available names in: ' data_path])
end

% Hydralic data %
% CONFIGURATION: Change the leaktionary name parameters to perform interpolation over other datasets.

network_name = 'Modena';
uncertainty = 0; uncertainty = erase(num2str(uncertainty),'.');

% Training %

min_leak_size_tr = 4; min_leak_size_tr = erase(num2str(min_leak_size_tr),'.');
max_leak_size_tr = 7; max_leak_size_tr = erase(num2str(max_leak_size_tr),'.');
data_purpose = 'training';

try 
load([hyd_data_path 'leaktionary_' network_name ...
      '_' min_leak_size_tr 'to' max_leak_size_tr ...
      '_' uncertainty 'p100_' data_purpose '.mat'])
Leaktionary_train = Leaktionary;
catch
    error(['The configured leaktionary name does not exist. Please, check the available names in: ' data_path])
end

% Testing %

min_leak_size_te = 4.5; min_leak_size_te = erase(num2str(min_leak_size_te),'.');
max_leak_size_te = 6.5; max_leak_size_te = erase(num2str(max_leak_size_te),'.');
data_purpose = 'testing';

try 
load([hyd_data_path 'leaktionary_' network_name ...
      '_' min_leak_size_te 'to' max_leak_size_te ...
      '_' uncertainty 'p100_' data_purpose '.mat'])
Leaktionary_test = Leaktionary;
catch
    error(['The configured leaktionary name does not exist. Please, check the available names in: ' data_path])
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
  
%% Load interpolated data %%

tau = 1000;
precision = 100;
int_type = 2; % 1 = sGSI, 2 = AW-GSI, 3 = GSI, 4 = real
leaky_node_v = double(leaky_node_v);

% Training %

int_type_str = ['with_RESID_tau_' num2str(tau) '_fnom_' uncertainty...
                'noise_' min_leak_size_tr 'to' max_leak_size_tr '_training'];

if exist([data_path '/interpolated_data/interpolation_' network_name '_' int_type_str '.mat'],'file')
    load([data_path '/interpolated_data/interpolation_' network_name '_' int_type_str '.mat']);
    disp(['Interpolation data loaded - Training: ' int_type_str]);
    res_leak_sGSIABS_M_train = res_leak_sGSIABS_M;
    res_leak_AW_M_train = res_leak_AW_M;
    res_leak_GSIABS_M_train = res_leak_GSIABS_M;
    res_leak_REAL_M_train = res_leak_REAL_M;
else
    error('Interpolation data not exists. Run launcher_interpolation.m to generate the interpolation data.')
end

% Testing %

int_type_str = ['with_RESID_tau_' num2str(tau) '_fnom_' uncertainty...
    'noise_' min_leak_size_te 'to' max_leak_size_te '_testing'];

if exist([data_path '/interpolated_data/interpolation_' network_name '_' int_type_str '.mat'],'file')
    load([data_path '/interpolated_data/interpolation_' network_name '_' int_type_str '.mat']);
    disp(['Interpolation data loaded - Testing: ' int_type_str]);
    res_leak_sGSIABS_M_test = res_leak_sGSIABS_M;
    res_leak_AW_M_test = res_leak_AW_M;
    res_leak_GSIABS_M_test = res_leak_GSIABS_M;
    res_leak_REAL_M_test = res_leak_REAL_M;
else
    error('Interpolation data not exists. Run launcher_interpolation.m to generate the interpolation data.')
end

switch int_type
    case 1
        res_leak_M_train = res_leak_sGSIABS_M_train;
        res_leak_M_test = res_leak_sGSIABS_M_test;
        disp('Smoothing GSI - ABS')
    case 2
        res_leak_M_train = res_leak_AW_M_train;
        res_leak_M_test = res_leak_AW_M_test;
        disp('AW-GSI')
    case 3
        res_leak_M_train = res_leak_GSIABS_M_train;
        res_leak_M_test = res_leak_GSIABS_M_test;
        disp('GSI - ABS')
    case 4
        res_leak_M_train = res_leak_REAL_M_train;
        res_leak_M_test = res_leak_REAL_M_test;
        disp('REAL')
    otherwise
        error('Not permited interpolation type!')
end

%% Dictionary Learning %%

total_sensors = [virtual_sensors]; % substitute virtual_sensors with sensors to use 0 VS

[~,~,ntimes_train] = size(res_leak_M_train);
[~,~,ntimes_test] = size(res_leak_M_test);

% Stack train data %

R_train=[]; labels_train=[];
for i=1:length(leaky_node_v)
R_train=[R_train reshape(res_leak_M_train(i,total_sensors,:),[length(total_sensors) ntimes_train])];%x_leak_M{i}(total_sensors,interval)];%tionary(:,interval)];%
labels_train=[labels_train i*ones(1,ntimes_train)];
end
R_train = -R_train;
 
% Stack test data %

R_test=[]; labels_test=[];
for i=1:length(leaky_node_v)
R_test=[R_test reshape(res_leak_M_test(i,total_sensors,:),[length(total_sensors) ntimes_test])];%x_leak_M{i}(total_sensors,interval)];%tionary(:,interval)];%
labels_test=[labels_test i*ones(1,ntimes_test)];
end
R_test = -R_test;

% DL settings %

cumacc1 = 0; cumacc2 = 0;

alpha = 4;                         % LC-KSVD alpha parameter
beta = 16;                         % LC-KSVD beta parameter
init_method = 2;                   % initialization for LC-KSVD, ie. trained 
                                   % atoms for each class and shared atoms                  
oc = 8;

for instance = 1:1 

disp(['Instance: ' num2str(instance) ' | Interp. type: ' num2str(int_type)])

% Random permutation of data matrices %

R_train_sel = []; labels_train_sel = []; indices_train = [];
for i=1:length(leaky_node_v)
    indices = randperm(ntimes_train)+(i-1)*(ntimes_train);
    indices_train = [indices_train indices];
end
R_train_sel = R_train(:,indices_train);
labels_train_sel = labels_train(:,indices_train);
indices_train_orig = indices_train;

R_test_sel = []; labels_test_sel = []; indices_test = [];
for i=1:length(leaky_node_v)
    indices = randperm(ntimes_test)+(i-1)*(ntimes_test);
    indices_test = [indices_test indices];
end
R_test_sel = R_test(:,indices_test);
labels_test_sel = labels_test(:,indices_test);

% Learning %

labels_train = labels_train_sel;
labels_test = labels_test_sel;
Y_train=normc(double(R_train_sel));
Y_test=normc(double(R_test_sel));
Y_train_orig = Y_train;
Y_test_orig = Y_test;

[m,N] = size(Y_train);
labels_list = unique(labels_train);  
c = length(labels_list);           % number of classes

n = oc * (c+1);                    % number of dictionary atoms
s = sqrt(n);                       % sparsity


% Form label matrix
H_train = zeros(c, N);
for cls = 1: c                  
    H_train(cls, labels_train == labels_list(cls)) = 1;
end

ctest = length(labels_list); 
labels_list_test = unique(labels_test);
H_test = zeros(ctest, size(Y_test,2));
for cls = 1: ctest                 
    H_test(cls, labels_test == labels_list_test(cls)) = 1;
end

% Form label consistency matrix
nc = floor(n/(c+1));                % evenly divide atoms per classes
nr = nc*(c+1);                      % total number of atoms (nr <= n)

Q = zeros(nr, N);
jj = 0;
for i = 1 : c                       % allocate atoms for each signal
    jc = find(H_train(i,:)==1);           % indices of signals from class i
    Q(jj+1:jj+nc,jc) = 1;
    jj = jj + nc;
end
Q(jj+1:jj+nc,:) = 1;                % shared dictionary

% Perform Label Consistent Dictionary Learning
tic
profile on -history
[W, D, A_DL] = clas_labelcon_dl(Y_train, H_train, Q, n, s, alpha, beta, init_method);
p_train = profile('info');
profile clear
time_lapse_training = toc;
%%
tic
profile on -history
[accuracy1,~,estimate,truth] = classification(Y_test, H_test, D, W, s);%
p_test = profile('info');
profile clear
time_lapse_testing = toc;

success0(instance)=0;     % classification is successful only if the estimation is equal with the node under fault or is inside the cummunity containing the node
success1(instance)=0;     % classification is successful if the estimation is equal to the node or any of its neighbors
success2(instance)=0;     % classification is successful if the estimation is equal to the node or any of its neighbors or any of the neighbirs' neighbors

labels_success = [labels_test];

for i=1:length(estimate)
    if labels_success(i)==estimate(i)
        success0(instance)=success0(instance)+1;
    end
    if any(unique([labels_success(i); neighbors(G,labels_success(i))])==estimate(i))
        success1(instance)=success1(instance)+1;
    end
    if any(unique(cell2mat(arrayfun(@(x) neighbors(G,x), [labels_success(i); neighbors(G, labels_success(i))],'UniformOutput',false)))==estimate(i))
        success2(instance)=success2(instance)+1;
    end
end
% transform into percentages
success0(instance)=success0(instance)*100/length(estimate)
success1(instance)=success1(instance)*100/length(estimate)
success2(instance)=success2(instance)*100/length(estimate)

end
