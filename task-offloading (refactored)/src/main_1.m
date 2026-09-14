clc;
clear;
close all;
cvx_quiet(true);

%% ========== SYSTEM PARAMETERS ==========

MaxLatency = 5 : 5 : 50; % Max acceptable latency [ms]

% MaxLatency = 20;

params.K = 30;                 % Number of single-antenna users
params.N_ant = 32;            % Number of antennas
params.N = 6;                 % Number of NFV-enabled nodes
params.L_R = 4;               % Number of RRHs

params.L = ones(params.K, 1);                     % Task size [Mega Cycles]
params.D = 5 * ones(params.K, 1);                 % Data size [20 Kbits/Normalized to 20 MHz Bandwidth]

% Channel between users and RRHs
params.H = randn(params.N_ant, params.L_R, params.K)/sqrt(2) + 1i * randn(params.N_ant, params.L_R, params.K)/sqrt(2);

% T^RAN = T_Ratio * MaxLatency
T_Ratio = 0.5;

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

num_latencies = numel(MaxLatency);
AR_JTO = zeros(num_latencies, 1);
AR_Radio_DTO = zeros(num_latencies, 1);
AR_NonRadio_DTO = zeros(num_latencies, 1);

CostPower_JTO = zeros(num_latencies, 1);
CostCompute_JTO = zeros(num_latencies, 1);

CostPower_DTO = zeros(num_latencies, 1);
CostCompute_DTO = zeros(num_latencies, 1);


mean_Tau_tx_JTO = zeros(num_latencies, 1);
mean_Tau_exe_JTO = zeros(num_latencies, 1);
mean_Tau_prop_JTO = zeros(num_latencies, 1);

mean_Tau_tx_DTO = zeros(num_latencies, 1);
mean_Tau_exe_DTO = zeros(num_latencies, 1);
mean_Tau_prop_DTO = zeros(num_latencies, 1);


parpool('local', 10)

parfor n = 1 : num_latencies
    max_delay = MaxLatency(n);
    fprintf('Max. Delay : %d \n', max_delay)
    [AR_JTO(n), Tau_tx_JTO, Tau_exe_JTO, Tau_prop_JTO, CostPower_JTO(n), CostCompute_JTO(n)] =...
        JTO_function(max_delay, params);
    
    mean_Tau_tx_JTO(n) = mean(Tau_tx_JTO);
    mean_Tau_exe_JTO(n) = mean(Tau_exe_JTO);
    mean_Tau_prop_JTO(n) = mean(Tau_prop_JTO);
    
    
    [AR_Radio_DTO(n), AR_NonRadio_DTO(n), Tau_tx_DTO, Tau_exe_DTO, Tau_prop_DTO, CostPower_DTO(n), CostCompute_DTO(n)] =...
        DTO_function(max_delay, T_Ratio, params);

    mean_Tau_tx_DTO(n) = mean(Tau_tx_DTO);
    mean_Tau_exe_DTO(n) = mean(Tau_exe_DTO);
    mean_Tau_prop_DTO(n) = mean(Tau_prop_DTO);
end
delete(gcp('nocreate'))

figure;
hold on
grid on
plot(MaxLatency, AR_JTO, '-o', 'LineWidth', 2, 'MarkerSize', 8)
plot(MaxLatency, AR_NonRadio_DTO, '-square', 'LineWidth', 2, 'MarkerSize', 8)
legend('JTO', 'DTO', 'Interpreter', 'latex', 'fontsize', 12, 'location', 'best')
xlabel('MaxDelay $T$ (ms)', 'Interpreter', 'latex', 'fontsize', 12)
ylabel('Acceptance Ratio', 'Interpreter', 'latex', 'fontsize', 12)