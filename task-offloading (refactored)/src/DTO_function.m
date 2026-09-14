function [AR_Radio, AR_Non_Radio, Tau_tx, Tau_exe, Tau_prop, Radio_Cost, Non_Radio_Cost] = ...
    DTO_function(max_delay, T_Ratio, params)
%% DTO.m - Disjoint Task Offloading and Resource Allocation
% This script implements a disjoint optimization framework for task offloading,
% where radio and computational resource allocation are performed separately.
%
% The algorithm:
% 1. Radio part feasibility analysis via CCP
% 2. Power allocation optimization for radio part
% 3. Non-radio part (computing) feasibility analysis
% 4. Heuristic task placement refinement
% 5. Computational resource allocation optimization

cvx_quiet(true);

%% ========== SYSTEM PARAMETERS ==========

H = params.H;
A = params.A;

[N, ~] = size(A);
[N_ant, L_R, K] = size(H);


% Number of possible paths between node 1 and each node j
path = ones(N, 1);
for j = 2:N
    path(j) = length(pathbetweennodes(A, 1, j));
end

%% ========== NETWORK PARAMETERS ==========
delta_link = params.delta_link;
C_node = params.C_node;
C_link = params.C_link;
C_front = params.C_front;

%% ========== LINK-TO-PATH INDICATOR ==========
% I_l2p(m, mm, b, n) = 1 if link (m, mm) is on path b to node n
I_l2p = zeros(N, N, max(path), N);
for m = 1:N
    for mm = 1:N
        for n = 1:N
            if n ~= 1
                paths2n = pathbetweennodes(A, 1, n);
                for b = 1:path(n)
                    bth_path = cell2mat(paths2n(b));
                    for i = 1:length(bth_path) - 1
                        if bth_path(i) == m && bth_path(i+1) == mm
                            I_l2p(m, mm, b, n) = 1;
                            I_l2p(mm, m, b, n) = 1;
                        end
                    end
                end
            else  % n == 1 (source node)
                I_l2p(m, mm, :, 1) = 0;
                I_l2p(1, 1, 1, 1) = 1;
            end
        end
    end
end

%% ========== PROPAGATION DELAY CALCULATION ==========
% prop(n, b) = propagation delay along path b from node 1 to node n
prop = nan(N, max(path));
for n = 1:N
    for b = 1:path(n)
        if n == 1
            prop(n, b) = 0;
        else
            prop(n, b) = 0;
            paths2n = pathbetweennodes(A, 1, n);
            bth_path = cell2mat(paths2n(b));
            for i = 1:length(bth_path) - 1
                prop(n, b) = prop(n, b) + delta_link(bth_path(i), bth_path(i+1));
            end
        end
    end
end

sigma = 1;  % Noise power

%% ========== TASK INITIALIZATION ==========
Tau = max_delay * ones(K, 1);
Tau_ran = T_Ratio * max_delay * ones(K, 1);
L = params.L;
D = params.D;

% load('channel.mat');
% H = H_initial(:, :, 1:K);

K_set = 1:K;
v = zeros(N, length(K_set));           % Task placement indicator
e = zeros(N, max(path), length(K_set)); % Path assignment

% Initially place all tasks at source node (node 1)
for k = K_set
    e(1, 1, k) = 1;
end

%% ========== INITIAL PROPAGATION DELAY ==========
Tau_prop = zeros(length(K_set), 1);
for k = K_set
    [n, b] = find(e(:, :, k) == 1);
    if n == 1
        Tau_prop(k) = 0;
    else
        paths2n = pathbetweennodes(A, 1, n);
        bth_path = cell2mat(paths2n(b));
        for i = 1:length(bth_path) - 1
            Tau_prop(k) = Tau_prop(k) + delta_link(bth_path(i), bth_path(i+1));
        end
    end
end

%% ========== INITIAL RRH ASSIGNMENT ==========
p_var = 1e-4 * ones(length(K_set), 1);
RRH_assign = zeros(length(K_set), 1);
for k = K_set
    channel_gain = zeros(L_R, 1);
    for l = 1:L_R
        channel_gain(l) = norm(H(:, l, k));
    end
    [~, RRH_assign(k)] = max(channel_gain);
end

% Subset of users per RRH
subset = cell(0);
for l = 1:L_R
    subset = [subset; {find(RRH_assign == l)}];
end

%% ========== TASK-TO-RRH INDICATOR ==========
I_t2r = zeros(L_R, length(K_set));
for l = 1:L_R
    for j = K_set
        if any(j == cell2mat(subset(l))')
            I_t2r(l, j) = 1;
        end
    end
end

%% ========== INITIAL SINR, RATE, AND TX DELAY ==========
H_user = zeros(N_ant, length(K_set));
for k = K_set
    H_user(:, k) = sum((H(:, :, k) * diag(I_t2r(:, k))).');
end

SINR = zeros(length(K_set), 1);
r = zeros(length(K_set), 1);
Tau_tx = zeros(length(K_set), 1);

for k = K_set
    interference = 0;
    for l = 1:L_R
        for j = cell2mat(subset(l))'
            if j ~= k
                interference = interference + (p_var(j) * abs((H(:, l, k)' * H(:, l, j))^2)) / (norm(H(:, l, j))^2);
            end
        end
    end
    SINR(k) = (norm(H_user(:, k))^2 * p_var(k)) / (interference + sigma^2);
    r(k) = log2(1 + SINR(k));
    Tau_tx(k) = D(k) / r(k);
end

% Initialize slack variables
s = Tau_tx + Tau_prop + 1e5 * max(Tau) * ones(length(K_set), 1);

%% ========== FEASIBILITY ANALYSIS FOR RADIO PART ==========
I_max = 100;

while true  % Radio feasibility loop
    count = 0;
    elastic_stack = zeros(length(K_set), I_max);

    while true  % ASM for radio part
        convergence_check = sum(s);
        count = count + 1;
        % disp(['ASM (Radio) Iteration: ', num2str(count)]);

        %% ===== RRH RE-ASSIGNMENT =====
        RRH_assign = zeros(length(K_set), 1);
        for k = K_set
            channel_gain = zeros(L_R, 1);
            for l = 1:L_R
                channel_gain(l) = norm(H(:, l, k));
            end
            [~, RRH_assign(k)] = max(channel_gain);
        end

        subset = cell(0);
        for l = 1:L_R
            subset = [subset; {find(RRH_assign == l)}];
        end

        I_t2r = zeros(L_R, length(K_set));
        for l = 1:L_R
            for j = K_set
                if any(j == cell2mat(subset(l))')
                    I_t2r(l, j) = 1;
                end
            end
        end

        %% ===== COMPUTE CCP TERMS =====
        % nabla_h
        nabla_h = zeros(length(K_set), length(K_set));
        for k = K_set
            for i = K_set
                interference = 0;
                numerator = 0;
                for l = 1:L_R
                    for j = cell2mat(subset(l))'
                        interference = interference + (p_var(j) * abs((H(:, l, k)' * H(:, l, j))^2)) / (norm(H(:, l, j))^2);
                    end
                    numerator = numerator + (I_t2r(l, i) * abs(H(:, l, k)' * H(:, l, i))^2) / norm(H(:, l, i))^2;
                end
                nabla_h(k, i) = numerator / (log(2) * (interference + sigma^2));
            end
        end

        % h and g
        h = zeros(length(K_set), 1);
        g = zeros(length(K_set), 1);
        for k = K_set
            int_h = 0;
            int_g = 0;
            for l = 1:L_R
                for j = cell2mat(subset(l))'
                    int_h = int_h + (p_var(j) * abs(H(:, l, k)' * H(:, l, j))^2) / (norm(H(:, l, j))^2);
                    if j ~= k
                        int_g = int_g + (p_var(j) * abs(H(:, l, k)' * H(:, l, j))^2) / (norm(H(:, l, j))^2);
                    end
                end
            end
            h(k) = log2(int_h + sigma^2);
            g(k) = log2(int_g + sigma^2);
        end

        % nabla_g
        nabla_g = zeros(length(K_set), length(K_set));
        for k = K_set
            for i = K_set
                interference = 0;
                numerator = 0;
                for l = 1:L_R
                    for j = cell2mat(subset(l))'
                        if j ~= k
                            interference = interference + (p_var(j) * abs((H(:, l, k)' * H(:, l, j))^2)) / (norm(H(:, l, j))^2);
                        end
                    end
                    numerator = numerator + (I_t2r(l, i) * abs(H(:, l, k)' * H(:, l, i))^2) / norm(H(:, l, i))^2;
                end
                if i == k
                    nabla_g(k, i) = 0;
                else
                    nabla_g(k, i) = numerator / (log(2) * (interference + sigma^2));
                end
            end
        end

        %% ===== h_tilde & Ie =====
        h_tilde = zeros(length(K_set), L_R, length(K_set));
        for k = K_set
            for l = 1:L_R
                for j = K_set
                    h_tilde(k, l, j) = (I_t2r(l, j) * abs(H(:, l, k)' * H(:, l, j))^2) / abs(H(:, l, j)' * H(:, l, j));
                end
            end
        end

        X = zeros(length(K_set), length(K_set));
        Y = zeros(length(K_set), 1);
        for k = K_set
            for j = K_set
                X(k, j) = sum(h_tilde(k, :, j));
            end
            Y(k, 1) = sum(h_tilde(k, :, k));
        end

        Ie = zeros(N, N, length(K_set));
        for m = 1:N
            for mm = 1:N
                for k = K_set
                    temp1 = reshape(I_l2p(m, mm, :, :), [max(path), N]);
                    temp2 = reshape(e(:, :, k), [N, max(path)]);
                    Ie(m, mm, k) = sum(sum(temp1 .* temp2'));
                end
            end
        end

        %% ===== POWER ALLOCATION VIA CCP =====
        % disp('Power Allocation (Radio Feasibility)');
        count_p = 0;
        while true
            count_p = count_p + 1;
            % disp(['CCP Iteration: ', num2str(count_p)]);

            cvx_begin
            variables p(length(K_set)) alpha_p
            minimize 1e2 * alpha_p
            subject to
            % Rate constraints
            for k = K_set
                temp = reshape(h_tilde(k, :, :), [L_R, length(K_set)]);
                log(sum(temp * p) + sigma^2) / log(2) - g(k) - nabla_g(k, :) * (p - p_var) >= ...
                    D(k) * inv_pos(Tau_ran(k) + s(k) + alpha_p);
                p(k) >= 0;
                p(k) <= 0.5;
            end
            alpha_p >= 0;

            % Front-haul capacity constraints
            for l = 1:L_R
                I_t2r(l, :) * (h + nabla_h * (p - p_var) - ...
                    log(X * p - diag(Y) * p + sigma^2 * ones(length(K_set), 1)) / log(2)) <= C_front(l, 1);
            end
            cvx_end

            % Compute concave rates
            r_concave = zeros(length(K_set), 1);
            for k = K_set
                temp = reshape(h_tilde(k, :, :), [L_R, length(K_set)]);
                r_concave(k) = log(sum(temp * p) + sigma^2) / log(2) - g(k) - nabla_g(k, :) * (p - p_var);
            end

            %% ===== UPDATE CCP TERMS =====
            % Update g
            for k = K_set
                int_g = 0;
                for l = 1:L_R
                    for i = cell2mat(subset(l))'
                        if i ~= k
                            int_g = int_g + (p(i) * abs((H(:, l, k)' * H(:, l, i))^2)) / (norm(H(:, l, i))^2);
                        end
                    end
                end
                g(k) = log2(int_g + sigma^2);
            end

            r_convex = h + nabla_h * (p - p_var) - g;

            % Update h
            h = zeros(length(K_set), 1);
            r_original = zeros(length(K_set), 1);
            for k = K_set
                int_h = 0;
                for l = 1:L_R
                    for i = cell2mat(subset(l))'
                        int_h = int_h + (p(i) * abs((H(:, l, k)' * H(:, l, i))^2)) / (norm(H(:, l, i))^2);
                    end
                end
                h(k) = log2(int_h + sigma^2);
                r_original(k) = h(k) - g(k);
            end

            % Update nabla_h
            nabla_h = zeros(length(K_set), length(K_set));
            for k = K_set
                for i = K_set
                    interference = 0;
                    numerator = 0;
                    for l = 1:L_R
                        for j = cell2mat(subset(l))'
                            interference = interference + (p(j) * abs((H(:, l, k)' * H(:, l, j))^2)) / (norm(H(:, l, j))^2);
                        end
                        numerator = numerator + (I_t2r(l, i) * abs(H(:, l, k)' * H(:, l, i))^2) / norm(H(:, l, i))^2;
                    end
                    nabla_h(k, i) = numerator / (log(2) * (interference + sigma^2));
                end
            end

            % Update nabla_g
            nabla_g = zeros(length(K_set), length(K_set));
            for k = K_set
                for i = K_set
                    interference = 0;
                    numerator = 0;
                    for l = 1:L_R
                        for j = cell2mat(subset(l))'
                            if j ~= k
                                interference = interference + (p(j) * abs((H(:, l, k)' * H(:, l, j))^2)) / (norm(H(:, l, j))^2);
                            end
                        end
                        numerator = numerator + (I_t2r(l, i) * abs(H(:, l, k)' * H(:, l, i))^2) / norm(H(:, l, i))^2;
                    end
                    if i == k
                        nabla_g(k, i) = 0;
                    else
                        nabla_g(k, i) = numerator / (log(2) * (interference + sigma^2));
                    end
                end
            end

            if norm(p_var - p) <= 1e-2 || count_p == 1
                p_var = p;
                break;
            end
            p_var = p;
        end

        Tau_tx = D ./ r_concave;

        %% ===== COMPUTE SLACK VARIABLES =====
        % disp('Computing Slack Variables (Radio)');
        cvx_begin
        variable s_new(length(K_set))
        minimize(sum(s_new))
        subject to
        for k = K_set
            s_new(k) >= Tau_tx(k) - Tau_ran(k);
            s_new(k) >= 0;
        end
        cvx_end
        s = s_new;

        % Check convergence
        elastic_stack(:, count) = s;
        if (convergence_check - sum(s)) / convergence_check < 1e-1 || count >= I_max || ...
                convergence_check - sum(s) < 1e-2 || sum(s) < 1e-5
            break;
        end
    end

    % Remove infeasible tasks if necessary
    if max(s) >= 1e-6
        if isscalar(K_set)
            disp('THE PROBLEM IS INFEASIBLE');
            AR_Radio = 0;
            AR_Non_Radio = 0;
            Tau_tx = nan;
            Tau_exe = nan;
            Tau_prop = nan;
            Radio_Cost = nan;
            Non_Radio_Cost = nan;
            return;
        end
        [~, k_reject] = max(s);
        K_set = 1:length(K_set) - 1;
        H(:, :, k_reject) = [];
        p(k_reject) = [];
        p_var(k_reject) = [];
        e(:, :, k_reject) = [];
        v(:, k_reject) = [];
        s(k_reject) = [];
        L(k_reject) = [];
        D(k_reject) = [];
        Tau(k_reject) = [];
        Tau_prop(k_reject) = [];
        Tau_tx(k_reject) = [];
        Tau_ran(k_reject) = [];
    else
        break;
    end
end

AR_Radio = length(K_set) / K;
Transmit_Power = sum(p);

%% ========== OPTIMIZATION FOR RADIO PART ==========
% disp('Power Allocation (Radio Optimization)');
count_p = 0;
while true
    count_p = count_p + 1;
    % disp(['CCP Iteration: ', num2str(count_p)]);

    cvx_begin
    variables p(length(K_set)) alpha_p
    minimize(sum(p) + 1e2 * alpha_p)
    subject to
    for k = K_set
        temp = reshape(h_tilde(k, :, :), [L_R, length(K_set)]);
        log(sum(temp * p) + sigma^2) / log(2) - g(k) - nabla_g(k, :) * (p - p_var) >= ...
            D(k) * inv_pos(Tau_ran(k) + alpha_p);
        p(k) >= 0;
        p(k) <= 0.5;
    end
    alpha_p >= 0;
    for l = 1:L_R
        I_t2r(l, :) * (h + nabla_h * (p - p_var) - ...
            log(X * p - diag(Y) * p + sigma^2 * ones(length(K_set), 1)) / log(2)) <= C_front(l, 1);
    end
    cvx_end

    cost_power = cvx_optval;

    % Compute concave rates
    r_concave = zeros(length(K_set), 1);
    for k = K_set
        temp = reshape(h_tilde(k, :, :), [L_R, length(K_set)]);
        r_concave(k) = log(sum(temp * p) + sigma^2) / log(2) - g(k) - nabla_g(k, :) * (p - p_var);
    end

    %% ===== UPDATE CCP TERMS =====
    % Update g
    for k = K_set
        int_g = 0;
        for l = 1:L_R
            for i = cell2mat(subset(l))'
                if i ~= k
                    int_g = int_g + (p(i) * abs((H(:, l, k)' * H(:, l, i))^2)) / (norm(H(:, l, i))^2);
                end
            end
        end
        g(k) = log2(int_g + sigma^2);
    end

    r_convex = h + nabla_h * (p - p_var) - g;

    % Update h
    h = zeros(length(K_set), 1);
    r_original = zeros(length(K_set), 1);
    for k = K_set
        int_h = 0;
        for l = 1:L_R
            for i = cell2mat(subset(l))'
                int_h = int_h + (p(i) * abs((H(:, l, k)' * H(:, l, i))^2)) / (norm(H(:, l, i))^2);
            end
        end
        h(k) = log2(int_h + sigma^2);
        r_original(k) = h(k) - g(k);
    end

    % Update nabla_h
    nabla_h = zeros(length(K_set), length(K_set));
    for k = K_set
        for i = K_set
            interference = 0;
            numerator = 0;
            for l = 1:L_R
                for j = cell2mat(subset(l))'
                    interference = interference + (p(j) * abs((H(:, l, k)' * H(:, l, j))^2)) / (norm(H(:, l, j))^2);
                end
                numerator = numerator + (I_t2r(l, i) * abs(H(:, l, k)' * H(:, l, i))^2) / norm(H(:, l, i))^2;
            end
            nabla_h(k, i) = numerator / (log(2) * (interference + sigma^2));
        end
    end

    % Update nabla_g
    nabla_g = zeros(length(K_set), length(K_set));
    for k = K_set
        for i = K_set
            interference = 0;
            numerator = 0;
            for l = 1:L_R
                for j = cell2mat(subset(l))'
                    if j ~= k
                        interference = interference + (p(j) * abs((H(:, l, k)' * H(:, l, j))^2)) / (norm(H(:, l, j))^2);
                    end
                end
                numerator = numerator + (I_t2r(l, i) * abs(H(:, l, k)' * H(:, l, i))^2) / norm(H(:, l, i))^2;
            end
            if i == k
                nabla_g(k, i) = 0;
            else
                nabla_g(k, i) = numerator / (log(2) * (interference + sigma^2));
            end
        end
    end

    if norm(p_var - p) <= 1e-2 || count_p == 10
        p_var = p;
        break;
    end
    p_var = p;
end

Tau_tx = D ./ r_concave;
Radio_Cost = sum(p);

%% ========== FEASIBILITY ANALYSIS FOR NON-RADIO PART ==========
s = Tau_tx + Tau_prop + 1e3 * max(Tau) * ones(length(K_set), 1);
c = zeros(length(K_set), 1);
Tau_exe = zeros(length(K_set), 1);

while true  % Non-radio feasibility loop
    if isempty(s)
        break;
    end
    count = 0;
    elastic_stack = zeros(length(K_set), I_max);

    while true  % ASM for non-radio part
        convergence_check = sum(s);
        count = count + 1;
        % disp(['ASM (Non-Radio) Iteration: ', num2str(count)]);

        %% ===== COMPUTATIONAL RESOURCE ALLOCATION =====
        % disp('Computational Resource Allocation (Feasibility)');
        cvx_begin
        variable c(length(K_set))
        minimize (sum(L .* inv_pos(c)))
        subject to
        for k = K_set
            c(k) >= L(k) / (Tau(k) + s(k) - Tau_prop(k) - Tau_ran(k));
            c(k) >= 0;
        end
        for n = 1:N
            temp = squeeze(e(n, :, :));
            sum(temp * c) <= C_node(n);
        end
        cvx_end
        Tau_exe = L ./ c;

        %% ===== COMPUTE SLACK VARIABLES =====
        % disp('Computing Slack Variables (Non-Radio)');
        cvx_begin
        variable s_new(length(K_set))
        minimize(sum(s_new))
        subject to
        for k = K_set
            s_new(k) >= Tau_prop(k) + Tau_exe(k) + Tau_ran(k) - Tau(k);
            s_new(k) >= 0;
        end
        cvx_end
        s = s_new;

        %% ===== HEURISTIC TASK PLACEMENT =====
        r_aug = zeros(N, max(path), length(K_set));
        for n = 1:N
            for b = 1:path(n)
                r_aug(n, b, :) = r_convex;
            end
        end

        % disp('Heuristic Task Placement');
        s = s + 1e-9 * rand(length(s), 1);  % For sorting stability

        Sorted_Vars = sort(s(s > 0));
        Sorted_Tasks = zeros(length(Sorted_Vars), 1);
        for j = 1:length(Sorted_Vars)
            Sorted_Tasks(j) = find(Sorted_Vars(j) == s);
        end

        C_tilde_node = zeros(N, 1);
        for k = Sorted_Tasks'
            % Residual computing capacity
            for n = 1:N
                temp = squeeze(e(n, :, :));
                C_tilde_node(n) = C_node(n) - sum(temp * c) + sum(c(k) * e(n, :, k)) - 1e-4;
            end

            % Residual link capacity
            C_tilde_link = zeros(N, N);
            for m = 1:N
                for mm = 1:N
                    temp = zeros(length(K_set), 1);
                    temp1 = reshape(I_l2p(m, mm, :, :), [max(path), N]);
                    for i = K_set
                        if i ~= k
                            temp2 = reshape(r_aug(:, :, i), [N, max(path)])';
                            temp(i) = sum(sum(temp1 .* reshape(e(:, :, i), [N, max(path)])' .* temp2));
                        end
                    end
                    C_tilde_link(m, mm) = C_link(m, mm) - sum(temp);
                end
            end

            % Find feasible nodes and paths
            Nodes_e = [];
            feasiblepaths = zeros(N, max(path));
            for n = 1:N
                feasiblepaths_b = [];
                for b = 1:path(n)
                    flag = [];
                    for m = 1:N
                        for mm = 1:N
                            if I_l2p(m, mm, b, n) == 1
                                flag = [flag; (r_convex(k) <= C_tilde_link(m, mm))];
                            end
                        end
                    end
                    if all(flag)
                        if ~ismember(n, Nodes_e)
                            Nodes_e = [Nodes_e, n];
                        end
                        feasiblepaths_b = [feasiblepaths_b, b];
                    end
                end
                feasiblepaths(n, 1:length(feasiblepaths_b)) = feasiblepaths_b;
            end

            Node_feasible = [];
            for n = Nodes_e
                if C_tilde_node(n) >= c(k) - 1e-3
                    Node_feasible = [Node_feasible, n];
                end
            end

            % Select best node-path pair
            Tau_exe_prop = nan(N, max(path));
            for n = Node_feasible
                for b = setdiff(feasiblepaths(n, :), 0)
                    Tau_exe_prop(n, b) = L(k) / C_tilde_node(n) + prop(n, b);
                end
            end
            [n_star, b_star] = find(Tau_exe_prop == min(min(Tau_exe_prop)));

            % Update placement and resource allocation
            Tau_prop(k) = prop(n_star, b_star);
            v(:, k) = zeros(N, 1);
            v(n_star, k) = 1;
            e(:, :, k) = zeros(N, max(path));
            e(n_star, b_star, k) = 1;

            s_tilde = Tau_ran(k) + Tau_prop(k) + L(k) / C_tilde_node(n_star) - Tau(k);
            if s_tilde < 0
                c(k) = L(k) / (Tau(k) - Tau_ran(k) - Tau_prop(k));
                s(k) = 0;
            else
                s(k) = s_tilde;
                c(k) = C_tilde_node(n_star);
            end
            Tau_exe(k) = L(k) / c(k);
        end

        % Check convergence
        elastic_stack(:, count) = s;
        if (convergence_check - sum(s)) / convergence_check < 1e-1 || count >= I_max || ...
                convergence_check - sum(s) < 1e-1
            break;
        end
    end

    % Remove infeasible tasks if necessary
    if max(s) >= 1e-6
        if isscalar(K_set)
            AR_Non_Radio = 0;
            Tau_tx = nan;
            Tau_exe = nan;
            Tau_prop = nan;
            Radio_Cost = nan;
            Non_Radio_Cost = nan;
            return;
        end
        [~, k_reject] = max(s);
        K_set = 1:length(K_set) - 1;
        H(:, :, k_reject) = [];
        c(k_reject) = [];
        p(k_reject) = [];
        p_var(k_reject) = [];
        e(:, :, k_reject) = [];
        v(:, k_reject) = [];
        s(k_reject) = [];
        L(k_reject) = [];
        D(k_reject) = [];
        Tau(k_reject) = [];
        Tau_prop(k_reject) = [];
        Tau_tx(k_reject) = [];
        Tau_exe(k_reject) = [];
        Tau_ran(k_reject) = [];
        r_convex(k_reject) = [];
    else
        break;
    end
end

%% ========== COMPUTATIONAL RESOURCE ALLOCATION OPTIMIZATION ==========
% disp('Computational Resource Allocation (Optimization)');
cvx_begin
variables c(length(K_set)) alpha_c
minimize(sum(pow_p(c, 3)) + 1e3 * max(C_node) * alpha_c)
subject to
for k = K_set
    c(k) >= L(k) * inv_pos(Tau(k) - Tau_prop(k) - Tau_tx(k) + alpha_c);
    c(k) >= 0;
end
alpha_c >= 0;
for n = 1:N
    temp = squeeze(e(n, :, :));
    sum(temp * c) <= C_node(n);
end
cvx_end

Tau_exe = L ./ c;

% Compute costs
c_n = zeros(N, 1);
for n = 1:N
    c_n(n) = v(n, :) * c;
end
Non_Radio_Cost = 1e-1 * sum(c_n.^3);
AR_Non_Radio = length(K_set) / K;
end