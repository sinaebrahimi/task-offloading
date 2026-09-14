clc;
clear;
close all;
cvx_quiet(true);

%% ========== SYSTEM PARAMETERS ==========

% MaxLatency = 5 : 5 : 50; % Max acceptable latency [ms]

MaxLatency = 20;

params.K = 30;                 % Number of single-antenna users
params.N_ant = 32;            % Number of antennas
params.N = 6;                 % Number of NFV-enabled nodes
params.L_R = 4;               % Number of RRHs

params.L = ones(params.K, 1);                    % Task size [Mega Cycles]
DataSizeRange = 5 * [1, 2, 3, 4];

% Channel between users and RRHs

params.H = randn(params.N_ant, params.L_R, params.K)/sqrt(2) + 1i * randn(params.N_ant, params.L_R, params.K)/sqrt(2);


%% ========== NETWORK TOPOLOGY & Parameters ==========
% Adjacency matrix for the network graph
params.A = [1 1 1 1 0 0;
     1 1 1 0 1 0;
     1 1 1 1 1 1;
     1 0 1 1 0 1;
     0 1 1 0 1 1;
     0 0 1 1 1 1];

params.delta_link = 20 * ones(params.N, params.N) + 0.00001 * randn(params.N, params.N);  % Link propagation delay
params.C_node = ones(params.N, 1);                % Computing capacity [GHz]
params.C_link = 20 * ones(params.N, params.N);           % Link capacity
params.C_front = 30 * ones(params.L_R, 1);        % Front-haul capacity


% Make delay and cost matrices symmetric with zero diagonal
params.delta_link = triu(params.delta_link, 1) + triu(params.delta_link, 1)';
params.delta_link = params.delta_link - diag(diag(params.delta_link));
params.delta_link = params.delta_link .* params.A;

params.C_link = params.C_link .* params.A;
params.C_link = params.C_link .* (10000 * eye(params.N)) + params.C_link;  % Large diagonal values

%% JTO Function Call

num_sizes = numel(DataSizeRange);
AR_JTO = zeros(num_sizes, 1);

CostPower_JTO = zeros(num_sizes, 1);
CostCompute_JTO = zeros(num_sizes, 1);



mean_Tau_tx_JTO = zeros(num_sizes, 1);
mean_Tau_exe_JTO = zeros(num_sizes, 1);
mean_Tau_prop_JTO = zeros(num_sizes, 1);


for n = 1 : num_sizes
    params.D = DataSizeRange(n) * ones(params.K, 1);                 % Data size [20 Kbits/Normalized to 20 MHz Bandwidth];
    fprintf('Task Data Size : %d \n', DataSizeRange(n))
    [AR_JTO(n), Tau_tx_JTO, Tau_exe_JTO, Tau_prop_JTO, CostPower_JTO(n), CostCompute_JTO(n)] =...
        JTO_function(MaxLatency, params);
    
    mean_Tau_tx_JTO(n) = mean(Tau_tx_JTO);
    mean_Tau_exe_JTO(n) = mean(Tau_exe_JTO);
    mean_Tau_prop_JTO(n) = mean(Tau_prop_JTO);
end

data = [mean_Tau_tx_JTO(:)'; mean_Tau_exe_JTO(:)'];

figure;

b = bar(data', 'stacked');

xticks(1:length(DataSizeRange));
xticklabels(string(DataSizeRange/50));

legend({'Average Tx. Delay', 'Average Exe. Delay'}, 'Orientation', 'horizontal',...
    'Location', 'northoutside', 'Interpreter', 'latex', 'fontsize', 12);

xlabel('$D$ (Mbits)', 'Interpreter', 'latex', 'fontsize', 12)
ylabel('Average Delay (ms)', 'Interpreter', 'latex', 'fontsize', 12)
grid on