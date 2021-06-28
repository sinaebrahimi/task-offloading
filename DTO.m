clc;
clear;
close all;
cvx_quiet true

%% Inputs:
K = input('Number of Single-Antenna Users: ');
Tau_max = input('Maximum Acceptable Latency of Each Task: ' );
T_scaler = input('Ratio of Ran Latency to the Maximum Accepteable latency: ' );
L_max = input('Computational Load of Each Task: ');
D_max = input('Data size of Each Task: ');
%% Initialization
%K=30;%input('number of single-antenna users: '); %number of single-antenna users
N_ant=32;%input('number of antennas: '); %number of antennas
N=6;%input('number of nfv-enabled nodes: '); %number of nodes
L_R=4;%input('number of RRHs: '); %number of RRHs
%%Definition of matrixes
%H=(1/sqrt(2))*randn(N_ant,L_R,K)+1i*(1/sqrt(2))*randn(N_ant,L_R,K); %complex gaussian   (H in SINR (12) formula)
while true
    flag_connectedgraph=zeros(N-1,1);
%     A=randi([0 1],N); %adj matrix of nfv_enabled nodes of Graph G
%     A=triu(A,1)+triu(A)';%making A symmetric (it is necessary for an adj matrix for an undirected graph)
%     A=A-diag(diag(A))+diag(ones(N,1)); %diag entries=1
%     A=[1     1     0     1     1     0
%         1     1     1     0     1     1
%         0     1     1     1     1     1
%         1     0     1     1     1     1
%         1     1     1     1     1     0
%         0     1     1     1     0     1];
           A=[1 1 1 1 0 0
              1 1 1 0 1 0
              1 1 1 1 1 1
              1 0 1 1 0 1
              0 1 1 0 1 1
              0 0 1 1 1 1];
    
    parfor n=2:N
        paths2n=pathbetweennodes(A,1,n);
        flag_connectedgraph(n-1,1)=isempty(paths2n);%if the node is not connected it returns an empty cell, which this answer is 1
    end
    
    if sum(flag_connectedgraph)==0
        break;
    end
end
path=ones(N,1);%number of possible paths between 1 and j (for each j)
parfor j=2:N
    path(j)=length(pathbetweennodes(A,1,j));%returns number of available paths between 1 and j
end
%% Graph's links and nodes setup
delta_link=20*ones(N,N)+ .00001*randn(N,N);
%C_node=10*rand(N,1); %C_n %Computational Capacity of node n %
C_node=5*ones(N,1);
%C_link=10*rand(N,N).*A;
C_link=20*ones(N,N);
%C_front=10*rand(L_R,1);%Fronthaul links capacity vector
C_front=30*ones(L_R,1);
cost_link=10*ones(N,N)+ .001*randn(N,N);
delta_link=triu(delta_link,1)+triu(delta_link)';%symmetric
delta_link=delta_link-diag(diag(delta_link));%make diag entries 0
delta_link=delta_link.*A;
cost_link=triu(cost_link,1)+triu(cost_link)';%symmetric
cost_link=cost_link-diag(diag(cost_link));%make diag entries 0
cost_link=cost_link.*A;
C_link=C_link.*A;
C_link=C_link.*(10000*eye(N))+C_link; %Diag should be big enough
%% Link to path indicator
I_l2p=zeros(N,N,max(path),N);  %(3) link to path indicator
for m=1:N
    for mm=1:N
        for n=1:N
            if n~=1
                paths2n=pathbetweennodes(A,1,n);
                for b=1:path(n)
                    bth_path=cell2mat(paths2n(b));
                    for i=1:length(bth_path)-1
                        if bth_path(i)==m && bth_path(i+1)==mm
                            I_l2p(m,mm,b,n)=1;
                            I_l2p(mm,m,b,n)=1;
                        end
                    end
                end
            else %n==1
                I_l2p(m,mm,:,1)=0;
                I_l2p(1,1,1,1)=1;
            end
        end
    end
end
%% Calculation of Delay for All paths in the Graph
prop=nan(N,max(path));
for n=1:N
    for b=1:path(n)
        if(n==1) %Caused bugs
            prop(n,b)=0;
        else
            prop(n,b)=0;
            paths2n=pathbetweennodes(A,1,n);%Cell Array containing all paths from 1 to n
            bth_path=cell2mat(paths2n(b));%bth path for task k
            %%%bth_path contians error:Index exceeds array bounds.
            for i=1:length(bth_path)-1
                prop(n,b)=prop(n,b)+delta_link(bth_path(i),bth_path(i+1)); %(8)
            end
        end
    end
end
%% Channel & Noise Power & Path Loss Generation
H_tilde=(1/sqrt(2))*randn(N_ant,L_R,K)+1i*(1/sqrt(2))*randn(N_ant,L_R,K); %complex gaussian   (H in SINR (12) formula)
radius=100*rand(K,1);
azimuth=2*pi*rand(K,1);
Location=[radius.*cos(azimuth), radius.*sin(azimuth)];
Location_site=[50,50;-50,50;-50,-50;50,-50];
Distance=zeros(K,L_R);
for k=1:K
    for l=1:L_R
        Distance(k,l)=10e-3*norm(Location(k,:)-Location_site(l,:));%Distance in Km
        %   Distance(k,l)=norm(Location(k,:)-Location_site(l,:));%Distance in meter
    end
end
Path_loss_dB=128.1+37.6*log10(Distance);
Path_loss_linear=10.^(Path_loss_dB/10);
% alpha=3; %Path loss exponent
% Path_loss_linear=Distance.^alpha;
H=zeros(N_ant,L_R,K);
for k=1:K
    for l=1:L_R
        %H(:,l,k)=(1/sqrt(Path_loss_linear(k,l))).*H_tilde(:,l,k);
        H(:,l,k)=sqrt(1/Path_loss_linear(k,l)).*H_tilde(:,l,k);
    end
end
sigma=2.0e7*10^(-180/10); %variance in SINR formula
sigma=sqrt(sigma);
H=H/sigma;
H_initial=H;
load('channel.mat')
sigma=1;
%% Simulation Setup
% Tau_scaler=30;
%Tau_ran_scaler=Tau_scaler-Tau_scaler/20;
%count_tau=Tau_scaler/20:Tau_scaler/20:Tau_scaler-Tau_scaler/20;
% count_tau=1:1:29;
count_tau=1;
Acceptance_ratio=zeros(length(count_tau),1);
Tx_power=zeros(length(count_tau),1);
Exe_cost=zeros(length(count_tau),1);
Acceptance_ratio_power=zeros(length(count_tau),1);
Task_offloading_cost=zeros(length(count_tau),1);
count_tau_integer=0;
for Tau_ran_scaler=count_tau
    count_tau_integer=count_tau_integer+1;
    H=H_initial;
    %Tau=rand(K,1); %max tolerable latency for kth task
    Tau=Tau_max*ones(K,1);
    % Tau_ran=Tau_ran_scaler*ones(K,1);
    %Tau_ran=Tau_ran_scaler*ones(K,1);
    Tau_ran=T_scaler*Tau_max*ones(K,1);
    %L=rand(K,1);%Load vector for kth task
    % L=10*ones(K,1);
    L = L_max*ones(K,1);
    %D=rand(K,1);%Data Size vector for kth task
    % D=5*ones(K,1);
    D = D_max*ones(K,1);
    %Tau=sort(Tau,'descend');
    %delta_link=.01*rand(N,N);%prop delay of EACH link %
    
    %% we change K to K_set here
    K_set=1:K;%Set of all users
    v=zeros(N,length(K_set)); %dec var for placing task k in node n
    e=zeros(N,max(path),length(K_set));
    
    % for k=K_set
    %     [n,b]=find(prop==max(max(prop)));
    %     v(n,k)=1; %determines offloaded node (only one) for task k
    %     e(n,b,k)=1;
    %  end
    % for k=K_set
    %     x=N*rand+1;%random variable
    %     x=floor(x);
    %     v(x,k)=1; %determines offloaded node (only one) for task k
    %     y=path(x)*rand+1;
    %     b=floor(y);
    %     e(x,b,k)=1; %just offloaded k in a single path
    % end
    for k=K_set
        v(1,k)=1; %determines offloaded node (only one) for task k
        e(1,1,k)=1; %just offloaded k in a single path
    end
    %% Tau_prop (8)
    
    Tau_prop=zeros(length(K_set),1);%Tau Prop for each task k
    
    parfor k=K_set
        [n,b]=find(e(:,:,k)==1);%n=offloaded node and bth offloaded path number
        if(n==1) %Caused bugs
            Tau_prop(k,1)=0;
        else
            paths2n=pathbetweennodes(A,1,n);%Cell Array containing all paths from 1 to n
            bth_path=cell2mat(paths2n(b));%bth path for task k
            %%%bth_path contians error:Index exceeds array bounds.
            for i=1:length(bth_path)-1
                Tau_prop(k,1)=Tau_prop(k,1)+delta_link(bth_path(i),bth_path(i+1)); %(8)
            end
            %paths2n={};%clear paths2n
        end
    end
    %% (12) calculating SINR for task k (12)
    %p=rand(length(K_set),1);%p (p_k)
    p_var=10^-4*ones(length(K_set),1);%Initialization of p for DC
    %p_var=10^-4* rand(length(K_set),1);
    %sigma=1; %variance in SINR formula
    
    RRH_assign=zeros(length(K_set),1);%each user is served by which RRH
    for k=K_set
        for l=1:L_R
            channel_gain(l)=norm(H(:,l,k));
        end
        [~,RRH_assign(k)]=max(channel_gain);%~:max channel gain  , which l has the most channel gain
    end %now we know which users use which RRH (l)
    %now we want to know which RRH serves UE k:
    subset=cell(0);%Each row will specifies which users are assigned to which RRH
    for l=1:L_R
        subset=[subset;{find(RRH_assign==l)}];%WARNING: some subset elements maybe empty (cell2mat will cause problem in case of being empty)
    end
    
    SINR=zeros(length(K_set),1);%Initialization SINR (12)
    r=zeros(length(K_set),1);%Initialization %r_k %achievable data rate for k (13)
    Tau_tx=zeros(length(K_set),1);%Initialization of transmission delay
    %% Calculation of I_l^k
    I_t2r=zeros(L_R,length(K_set)); % task k to RRH L_R indicator
    
    for l=1:L_R
        for j=K_set
            if sum(j==cell2mat(subset(l))')
                I_t2r(l,j)=1;
            end
        end
    end
    %% Calculation of h_k
    H_user=zeros(N_ant,length(K_set));
    for k=K_set
        H_user(:,k)=sum((H(:,:,k)*diag(I_t2r(:,k))).');
    end
    %% Calculation of (12) and (13) and Tau_tx
    parfor k=K_set
        int=0;%accumulator for interference of k (in denominator of SINR formula (12))
        for l=1:L_R
            for j=cell2mat(subset(l))'
                if (j~=k)
                    int=int+(p_var(j)*abs((H(:,l,k)'*H(:,l,j))^2))/(norm(H(:,l,j))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                end
            end
            %numerator=numerator+ismember(k,cell2mat(subset(l)))*H(:,l,k);%numerator of SINR formula  %dont know why it gets 0 most of the time
        end
        numerator=norm(H_user(:,k))^2;
        SINR(k)=(numerator*p_var(k))/(int+sigma^2);%SINR Formula
        r(k)=log2(1+SINR(k));%formula (13)
        Tau_tx(k)=D(k)/r(k);%Calculating transmission delay
    end %some values of SINR is 0, but it shouldn't be
    s=Tau_tx+Tau_prop+ones(length(K_set),1)*1.0e20*max(Tau); %s_k (1000)
    %% begin while
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ADMISSION CONTROL FOR DISJOINT METHOD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_max=100; %max allowed iteration of feasibilty problem algorithm
    K_reject=[];
    count_convergence=1;
    while true
        count=0;%counter of while loop for convergence of feasibility problem (16)
        elastic_stack=zeros(length(K_set),I_max);
        %rate_lower_bound=zeros(length(K_set),I_max);
        
        while true
            convergence_check=sum(s);
            count=count+1;
            disp('count');
            disp(count);
            %% ADMISSION CONTROL ALGORITHM FOR DISJOINT METHOD - RADIO PART
            channel_gain=zeros(L_R,1);
            RRH_assign=zeros(length(K_set),1);%each user is served by which RRH
            for k=K_set
                for l=1:L_R
                    channel_gain(l)=norm(H(:,l,k));
                end
                [~,RRH_assign(k)]=max(channel_gain);%~:max channel gain  , which l has the most channel gain
            end %now we know which users use which RRH (l)
            %now we want to know which RRH serves UE k:
            subset=cell(0);%Each row will specifies which users are assigned to which RRH
            for l=1:L_R
                subset=[subset;{find(RRH_assign==l)}];%WARNING: some subset elements maybe empty (cell2mat will cause problem in case of being empty)
            end
            
            %% Calculation of I_l^k
            I_t2r=zeros(L_R,length(K_set)); % task k to RRH L_R indicator
            
            for l=1:L_R
                for j=K_set
                    if sum(j==cell2mat(subset(l))')
                        I_t2r(l,j)=1;
                    end
                end
            end
            %% calculation of h(p) g(p) nabla_g(p) nabla_h(p)
            nabla_h=zeros(length(K_set),length(K_set)); %h(p) in (26) %based on (22) for convex approximation
            %(26)
            for k=K_set
                for i=K_set
                    int=0;%accumulator for interference of k (in denominator of SINR formula (12))
                    numerator=0;%accumulator for numerator of  SINR Formula
                    for l=1:L_R
                        for j=cell2mat(subset(l))'
                            int=int+(p_var(j)*abs((H(:,l,k)'*H(:,l,j))^2))/(norm(H(:,l,j))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                        end
                        numerator=numerator+(I_t2r(l,i)*abs(H(:,l,k)'*H(:,l,i))^2)/norm(H(:,l,i))^2;%numerator of SINR formula  %dont know why it gets 0 most of the time
                    end
                    nabla_h(k,i)=numerator/(log(2)*(int+sigma^2));%(26) %based on (22) for convex approximation
                end
            end
            h=zeros(length(K_set),1);
            g=zeros(length(K_set),1);
            for k=K_set
                int_h=0;%accumulator for interference of k (in denominator of SINR formula (12))
                int_g=0;
                for l=1:L_R
                    for j=cell2mat(subset(l))'
                        int_h=int_h+(p_var(j)*abs(H(:,l,k)'*H(:,l,j))^2)/(norm(H(:,l,j))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                        if j~=k
                            int_g=int_g+(p_var(j)*abs(H(:,l,k)'*H(:,l,j))^2)/(norm(H(:,l,j))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                        end
                    end
                end
                h(k)=log2(int_h+sigma^2);%h(p) in (22)
                g(k)=log2(int_g+sigma^2);
            end
            
            %% g(p) in (26) based on (22) for convex approx.
            nabla_g=zeros(length(K_set),length(K_set));
            for k=K_set
                for i=K_set
                    int=0;%accumulator for interference of k (in denominator of SINR formula (12))
                    numerator=0;%accumulator for numerator of  SINR Formula
                    for l=1:L_R
                        for j=cell2mat(subset(l))'
                            if j~=k
                                int=int+(p_var(j)*abs((H(:,l,k)'*H(:,l,j))^2))/(norm(H(:,l,j))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                            end
                        end
                        numerator=numerator+(I_t2r(l,i)*abs(H(:,l,k)'*H(:,l,i))^2)/norm(H(:,l,i))^2;%numerator of SINR formula  %dont know why it gets 0 most of the time
                    end
                    if i==k
                        nabla_g(k,i)=0;
                    else
                        nabla_g(k,i)=numerator/(log(2)*(int+sigma^2));%(26) %based on (22) for convex approximation
                    end
                end
            end
            %% Calculation of h_tilde & Ie
            
            h_tilde=zeros(length(K_set),L_R,length(K_set));%Fraction of h and g in (22)
            
            for k=K_set
                for l=1:L_R
                    for j=K_set
                        h_tilde(k,l,j)=(I_t2r(l,j)*abs(H(:,l,k)'*H(:,l,j))^2)/abs(H(:,l,j)'*H(:,l,j));
                    end
                end
            end
            
            X=zeros(length(K_set),length(K_set));%matrices for h(p) - g(p) in C3 of (27)
            Y=zeros(length(K_set),1);
            
            for k=K_set
                for j=K_set
                    X(k,j)=sum(h_tilde(k,:,j));
                end
                Y(k,1)=sum(h_tilde(k,:,k));
            end
            Ie=zeros(N,N,length(K_set)); %parameter for inner summations of C3 in (27)
            for m=1:N
                for mm=1:N
                    for k=K_set
                        temp1=zeros(max(path),N);
                        temp2=zeros(N,max(path));
                        temp1(:,:)=I_l2p(m,mm,:,:);
                        temp2(:,:)=e(:,:,k);
                        Ie(m,mm,k)=sum(sum(temp1.*(temp2)'));
                    end
                end
            end
            nabla_Tau_tx=zeros(length(K_set),length(K_set));%
            for k=K_set
                nabla_Tau_tx(k,:)=(-D'.*(nabla_h(k,:)-nabla_g(k,:)))./((h(k)-g(k))^2);    %(23)
                %   nabla_Tau_tx(k,:)=(-D'.*(nabla_h(k,:)-nabla_g(k,:)));    %(23)
            end
            %% Calculation of p_var
            disp('POWER ALLOCATION PROBLEM');
            count_p=0;
            while true
                count_p=count_p+1;
                disp('count_p');
                disp(count_p);
                cvx_begin
                variable p(length(K_set))
                %  minimize sum(nabla_Tau_tx*(p-p_var))
                %maximize sum((log(X*p + ones(length(K_set),1)*(sigma^2) )/log(2))-g-nabla_g*(p - p_var))
                subject to
                
                for k=K_set
                    temp=zeros(L_R,length(K_set));
                    temp(:,:)=h_tilde(k,:,:);
                    log(sum(temp*p)+sigma^2)/log(2)-g(k)-nabla_g(k,:)*(p - p_var) >=  D(k)/(Tau_ran(k)+s(k)) %h(p) in (22)  %C1 in (27)
                    p(k)>=0 %c5 in (27)
                    p(k)<=0.5
                end
                
                
                for l=1:L_R
                    I_t2r(l,:)*(h+nabla_h*(p - p_var)-log(X*p - diag(Y)*p + ones(length(K_set),1)*(sigma^2) )/log(2)) <= C_front(l,1) %C4 in (27)
                end
                cvx_end
                disp(cvx_status)
                %             disp('cvx_optval')
                %             disp(cvx_optval)
                
                %   r_convex=h+nabla_h*(p-p_var)-g;
                r_concave=zeros(length(K_set),1);
                for k=K_set
                    temp=zeros(L_R,length(K_set));
                    temp(:,:)=h_tilde(k,:,:);
                    r_concave(k,1)=log(sum(temp*p)+sigma^2)/log(2)-g(k)-nabla_g(k,:)*(p - p_var);
                end
                Tau_tx=D./r_concave;
                disp('Tau_tx<=Tau_ran+s');
                disp((Tau_tx<=Tau_ran+s)');
                disp('Tau_tx-Tau_ran-s');
                disp((Tau_tx-Tau_ran-s)');
                %% Calculated concave approximation of r_k for obtained values of p (p-var)
                g(k)=log2(int_g+sigma^2);%g(p) in (22)
                for k=K_set
                    %  int_h=0;%accumulator for interference of k (in denominator of SINR formula (12))
                    int_g=0;
                    for l=1:L_R
                        for i=cell2mat(subset(l))'
                            %         int_h=int_h+(p(i)*abs((H(:,l,k)'*H(:,l,i))^2))/(norm(H(:,l,i))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                            if i~=k
                                int_g=int_g+(p(i)*abs((H(:,l,k)'*H(:,l,i))^2))/(norm(H(:,l,i))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                            end
                        end
                    end
                    g(k)=log2(int_g+sigma^2);%g(p) in (22)
                end
                
                r_convex=h+nabla_h*(p-p_var)-g;
                
                
                h=zeros(length(K_set),1);
                %g=zeros(length(K_set),1);
                r_original=zeros(length(K_set),1);
                
                parfor k=K_set
                    int_h=0;%accumulator for interference of k (in denominator of SINR formula (12))
                    %int_g=0;
                    for l=1:L_R
                        for i=cell2mat(subset(l))'
                            int_h=int_h+(p(i)*abs((H(:,l,k)'*H(:,l,i))^2))/(norm(H(:,l,i))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                            %       if i~=k
                            %          int_g=int_g+(p(i)*abs((H(:,l,k)'*H(:,l,i))^2))/(norm(H(:,l,i))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                            %     end
                        end
                    end
                    h(k)=log2(int_h+sigma^2);%h(p) in (22)
                    %g(k)=log2(int_g+sigma^2);%g(p) in (22)
                    r_original(k,1)=h(k)-g(k);
                    temp=zeros(L_R,length(K_set));
                    temp(:,:)=h_tilde(k,:,:);
                end
                
                nabla_h=zeros(length(K_set),length(K_set)); %h(p) in (26) %based on (22) for convex approximation
                %(26)
                for k=K_set
                    for i=K_set
                        int=0;%accumulator for interference of k (in denominator of SINR formula (12))
                        numerator=0;%accumulator for numerator of  SINR Formula
                        for l=1:L_R
                            for j=cell2mat(subset(l))'
                                int=int+(p(j)*abs((H(:,l,k)'*H(:,l,j))^2))/(norm(H(:,l,j))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                            end
                            numerator=numerator+(I_t2r(l,i)*abs(H(:,l,k)'*H(:,l,i))^2)/norm(H(:,l,i))^2;%numerator of SINR formula  %dont know why it gets 0 most of the time
                        end
                        nabla_h(k,i)=numerator/(log(2)*(int+sigma^2));%(26) %based on (22) for convex approximation
                    end
                end
                
                nabla_g=zeros(length(K_set),length(K_set));
                for k=K_set
                    for i=K_set
                        int=0;%accumulator for interference of k (in denominator of SINR formula (12))
                        numerator=0;%accumulator for numerator of  SINR Formula
                        for l=1:L_R
                            for j=cell2mat(subset(l))'
                                if j~=k
                                    int=int+(p(j)*abs((H(:,l,k)'*H(:,l,j))^2))/(norm(H(:,l,j))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                                end
                            end
                            numerator=numerator+(I_t2r(l,i)*abs(H(:,l,k)'*H(:,l,i))^2)/norm(H(:,l,i))^2;%numerator of SINR formula  %dont know why it gets 0 most of the time
                        end
                        if i==k
                            nabla_g(k,i)=0;
                        else
                            nabla_g(k,i)=numerator/(log(2)*(int+sigma^2));%(26) %based on (22) for convex approximation
                        end
                    end
                end
                %r_concave=h_p_var-g-nabla_g*(p-p_var);
                %   DC_objective(count_p,count)=sum(nabla_Tau_tx*(p-p_var));
                if norm(p_var-p)<=10e-3 || count_p==1
                    p_var=p;
                    break;
                end
                p_var=p;
                %p_var'
            end
            Tau_tx=D./r_concave;
            %% solving s(k)
            disp('ELASTIC VARIABLES PROBLEM');
            cvx_begin
            variable sss(length(K_set))
            minimize ( sum(sss) )
            subject to
            for k=K_set
                sss(k)>=Tau_tx(k)-Tau_ran(k)
                sss(k) >= 0
            end
            %         for k=K_reject
            %             sss(k)==0
            %         end
            cvx_end
            s=sss;
            disp(cvx_status)
            disp('elasticization variable ');
            disp(s);
            disp('Tau_tx<=Tau_ran+s');
            disp((Tau_tx<=Tau_ran+s)');
            disp('Tau_tx-Tau_ran-s');
            disp((Tau_tx-Tau_ran-s)');
            convergence(count_convergence)=sum(s);
            count_convergence=count_convergence+1;
            elastic_stack(:,count)=s;
            disp('The ratio of objective improvement')
            disp((convergence_check-sum(s))/convergence_check)
            if (convergence_check-sum(s))/convergence_check <1.0e-1|| count>=I_max || convergence_check-sum(s)<1.0e-2
                break;
            end
        end
        if max(s)>=.000001
            K_set=1:length(K_set)-1;
            [~,k_reject]=max(s);
            H(:,:,k_reject)=[];
            %c(k_reject)=[];
            p(k_reject)=[];
            p_var(k_reject)=[];
            e(:,:,k_reject)=[];
            v(:,k_reject)=[];
            s(k_reject)=[];
            L(k_reject)=[];
            D(k_reject)=[];
            Tau(k_reject)=[];
            Tau_prop(k_reject)=[];
            Tau_tx(k_reject)=[];
            Tau_ran(k_reject)=[];
            %Tau_exe(k_reject)=[];
        else
            break;
        end
    end
    Acceptance_ratio_power(count_tau_integer,1)=length(K_set)/K;
    %% POWER ALLOCATION ALGORITHM FOR DISJOINT METHOD - OPTIMIZATION
    disp('power problem ');
    count_p=0;
    while true
        count_p=count_p+1;
        disp('count_p');
        disp(count_p);
        cvx_begin
        variable p(length(K_set))
        %  minimize ( sum(nabla_Tau_tx*(p-p_var)) )
        %minimize ( -100*ones(1,length(K_set))*(p-p_var))
        minimize (sum(p))
        subject to
        
        for k=K_set
            temp=zeros(L_R,length(K_set));
            temp(:,:)=h_tilde(k,:,:);
            log(sum(temp*p)+sigma^2)/log(2)-g(k)-nabla_g(k,:)*(p - p_var) >=  D(k)/(Tau_ran(k)) %h(p) in (22)  %C1 in (27)
            p(k)>=0 %c5 in (27)
            p(k)<=0.5
        end
        
        for l=1:L_R
            I_t2r(l,:)*(h+nabla_h*(p - p_var)-log(X*p - diag(Y)*p + ones(length(K_set),1)*(sigma^2) )/log(2)) <= C_front(l,1) %C4 in (27)
        end
        cvx_end
        disp(cvx_status)
        disp('cvx_optval')
        disp(cvx_optval)
        cost_power=cvx_optval;
        
        %   r_convex=h+nabla_h*(p-p_var)-g;
        r_concave=zeros(length(K_set),1);
        parfor k=K_set
            temp=zeros(L_R,length(K_set));
            temp(:,:)=h_tilde(k,:,:);
            r_concave(k,1)=log(sum(temp*p)+sigma^2)/log(2)-g(k)-nabla_g(k,:)*(p - p_var);
        end
        Tau_tx=D./r_concave;
        disp('Tau_tx<=Tau_ran+s');
        disp((Tau_tx<=Tau_ran+s)');
        disp('Tau_tx-Tau_ran-s');
        disp((Tau_tx-Tau_ran-s)');
        %% Calculated concave approximation of r_k for obtained values of p (p-var)
        g(k)=log2(int_g+sigma^2);%g(p) in (22)
        for k=K_set
            %  int_h=0;%accumulator for interference of k (in denominator of SINR formula (12))
            int_g=0;
            for l=1:L_R
                for i=cell2mat(subset(l))'
                    %         int_h=int_h+(p(i)*abs((H(:,l,k)'*H(:,l,i))^2))/(norm(H(:,l,i))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                    if i~=k
                        int_g=int_g+(p(i)*abs((H(:,l,k)'*H(:,l,i))^2))/(norm(H(:,l,i))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                    end
                end
            end
            g(k)=log2(int_g+sigma^2);%g(p) in (22)
        end
        
        r_convex=h+nabla_h*(p-p_var)-g;
        
        
        h=zeros(length(K_set),1);
        %g=zeros(length(K_set),1);
        r_original=zeros(length(K_set),1);
        
        for k=K_set
            int_h=0;%accumulator for interference of k (in denominator of SINR formula (12))
            %int_g=0;
            for l=1:L_R
                for i=cell2mat(subset(l))'
                    int_h=int_h+(p(i)*abs((H(:,l,k)'*H(:,l,i))^2))/(norm(H(:,l,i))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                    %       if i~=k
                    %          int_g=int_g+(p(i)*abs((H(:,l,k)'*H(:,l,i))^2))/(norm(H(:,l,i))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                    %     end
                end
            end
            h(k)=log2(int_h+sigma^2);%h(p) in (22)
            %g(k)=log2(int_g+sigma^2);%g(p) in (22)
            r_original(k,1)=h(k)-g(k);
            temp=zeros(L_R,length(K_set));
            temp(:,:)=h_tilde(k,:,:);
        end
        
        nabla_h=zeros(length(K_set),length(K_set)); %h(p) in (26) %based on (22) for convex approximation
        %(26)
        for k=K_set
            for i=K_set
                int=0;%accumulator for interference of k (in denominator of SINR formula (12))
                numerator=0;%accumulator for numerator of  SINR Formula
                for l=1:L_R
                    for j=cell2mat(subset(l))'
                        int=int+(p(j)*abs((H(:,l,k)'*H(:,l,j))^2))/(norm(H(:,l,j))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                    end
                    numerator=numerator+(I_t2r(l,i)*abs(H(:,l,k)'*H(:,l,i))^2)/norm(H(:,l,i))^2;%numerator of SINR formula  %dont know why it gets 0 most of the time
                end
                nabla_h(k,i)=numerator/(log(2)*(int+sigma^2));%(26) %based on (22) for convex approximation
            end
        end
        
        nabla_g=zeros(length(K_set),length(K_set));
        for k=K_set
            for i=K_set
                int=0;%accumulator for interference of k (in denominator of SINR formula (12))
                numerator=0;%accumulator for numerator of  SINR Formula
                for l=1:L_R
                    for j=cell2mat(subset(l))'
                        if j~=k
                            int=int+(p(j)*abs((H(:,l,k)'*H(:,l,j))^2))/(norm(H(:,l,j))^2);%denominator of SINR formula    %%% .8 or .^ ?????? problem
                        end
                    end
                    numerator=numerator+(I_t2r(l,i)*abs(H(:,l,k)'*H(:,l,i))^2)/norm(H(:,l,i))^2;%numerator of SINR formula  %dont know why it gets 0 most of the time
                end
                if i==k
                    nabla_g(k,i)=0;
                else
                    nabla_g(k,i)=numerator/(log(2)*(int+sigma^2));%(26) %based on (22) for convex approximation
                end
            end
        end
        if norm(p_var-p)<=10e-3 || count_p==10
            p_var=p;
            break;
        end
        p_var=p;
        %p_var'
    end
    Tau_tx=D./r_concave;
    Tx_power(count_tau_integer,1)=cost_power;
    
    
    %% ADMISSION CONTROL AND TASK OFFLOADING ALGORITHM FOR DISJOINT METHID
    s=Tau_tx+Tau_prop+ones(length(K_set),1)*1.0e3*max(Tau);
    while true
        if isempty(s)==1
            
            break;
        end
        count=0;%counter of while loop for convergence of feasibility problem (16)
        elastic_stack=zeros(length(K_set),I_max);
        %rate_lower_bound=zeros(length(K_set),I_max);
        
        while true
            convergence_check=sum(s);
            count=count+1;
            disp('count');
            disp(count);
            
            
            %% cvx for (19)
            disp('computational capacity problem ');
            %constraint_bound=Tau+s-Tau_prop-Tau_tx %%%temporary
            cvx_begin
            variable c(length(K_set))
            %minimize ( sum(v*pow_p(c,3)) )
            minimize ( sum(L(K_set).*inv_pos(c)) )
            subject to
            for k=K_set
                c(k) >= (L(k))/(Tau(k)+s(k)-Tau_prop(k)-Tau_ran(k))
                c(k) >= 0
            end
            for n=1:N
                v(n,:)*c <= C_node(n)
            end
            cvx_end %end of (20)
            disp(cvx_status)
            Tau_exe=L(K_set)./c;% (Formula (7) for execution delay) (Load/capacity
            %     disp('Tau_tx+Tau_prop+Tau_exe<=Tau+s');
            %     disp(Tau_tx+Tau_prop+Tau_exe<=Tau+s);
            %     disp('Tau_exe ');
            %     disp(Tau_exe);
            %% MOSEK ... Optimization of binary variables
            Idelta=zeros(N,max(path));
            r_aug=zeros(N,max(path),length(K_set));
            for n=1:N
                for b=1:path(n)
                    Idelta(n,b)=sum(sum(I_l2p(:,:,b,n).*delta_link))/2;
                    r_aug(n,b,:)=r_convex;
                end
            end
            
            Idelta_aug=zeros(N,max(path),length(K_set));
            for k=K_set
                Idelta_aug(:,:,k)=Idelta;
            end
            
            I_l2p_aug=zeros(N,N,N,max(path),length(K_set));
            for m=1:N
                for mm=1:N
                    temp=zeros(max(path),N);
                    
                    temp(:,:)=I_l2p(m,mm,:,:);
                    for k=K_set
                        I_l2p_aug(m,mm,:,:,k)=temp';
                    end
                end
            end
            
            v_aug=zeros(N,max(path),length(K_set));
            for n=1:N
                for b=1:max(path)
                    for k=K_set
                        v_aug(n,b,k)=v(n,k);
                    end
                end
            end
            
            
            %% solving (30) for v,e,gamma
            disp('BINARY VARIABLES PROBLEM');
            cvx_begin
            cvx_solver MOSEK
            variable e_var(N,max(path),length(K_set)) binary
            variable gamma_var(N,max(path),length(K_set)) binary
            variable v_var(N,length(K_set)) binary
            
            minimize sum(sum(sum(e_var.*Idelta_aug)))
            subject to
            for k=K_set
                %%%%%%%%%%%%%%%%%%%%%%%%C1-a DCP error
                sum(sum(e_var(:,:,k).* Idelta)) <= Tau(k)+s(k)-Tau_ran(k)-Tau_exe(k) %C1-f
                sum(sum(e_var(:,:,k)))==1 %C8
                sum(v_var(:,k)) == 1 %C7
                sum(sum(gamma_var(:,:,k))) == 1 %C9
                
                for n=1:N
                    for b=1:path(n)
                        gamma_var(n,b,k) <= v_var(n,k)+1-e_var(n,b,k) %C15
                        v_var(n,k)       <= gamma_var(n,b,k)+1-e_var(n,b,k) %C16
                        gamma_var(n,b,k) <= e_var(n,b,k) %C17
                    end
                    if path(n)<max(path) %making unnecessary variables zero
                        for b=path(n)+1:max(path)
                            e_var(n,b,k)==0
                            gamma_var(n,b,k)==0
                        end
                    end
                end
            end
            for n=1:N
                v_var(n,:)*c <= C_node(n) %C2
            end
            
            for m=1:N
                for mm=1:N
                    temp=zeros(N,max(path),length(K_set));
                    temp(:,:,:)=I_l2p_aug(m,mm,:,:,:);
                    sum(sum(sum(r_aug.*temp.*e_var))) <= C_link(m,mm) %C3
                end
            end
            cvx_end
            disp(cvx_status)
            gamma_var=full(gamma_var); % MATLAB can not put sparse matrices in N>2 Dimensional arrays, which cause error
            e_var=full(e_var);
            v_var=full(v_var);
            %% Tau_prop for obtained v_var & e_var
            for k=K_set
                Tau_prop(k)=sum(sum(gamma_var(:,:,k).* Idelta));
            end
            v=v_var;
            e=e_var;
            
            for n=1:N
                for k=K_set
                    if abs(v(n,k)-1)<.5
                        v(n,k)=1;
                    else
                        v(n,k)=0;
                    end
                    for b=1:max(path)
                        if abs(e(n,b,k)-1)<.5
                            e(n,b,k)=1;
                        else
                            e(n,b,k)=0;
                        end
                    end
                end
            end
            
            
            %     disp('Tau_prop ');
            %     disp(Tau_prop);
            disp('Tau_tx+Tau_prop+Tau_exe<=Tau+s');
            disp((Tau_tx+Tau_prop+Tau_exe<=Tau+s)');
            disp('Tau_tx+Tau_prop+Tau_exe-Tau-s');
            disp((Tau_tx+Tau_prop+Tau_exe-Tau-s)');
            %% solving s(k)
            disp('elastic variables problem ');
            cvx_begin
            variable sss(length(K_set))
            minimize ( sum(sss) )
            subject to
            for k=K_set
                sss(k)>=Tau_prop(k)+Tau_exe(k)+Tau_ran(k)-Tau(k)
                sss(k) >= 0
            end
            for k=K_reject
                sss(k)==0
            end
            cvx_end
            s=sss;
            disp(cvx_status)
            %     for k=K_set
            %         s(k)=max(Tau_exe(k)+Tau_prop(k)+Tau_tx(k)-Tau(k),0);
            %     end
            %constraint_bound_updated=Tau+s-Tau_tx-Tau_prop
            disp('elasticization variable ');
            disp(s);
            %% ASM Modification Algorithm
            disp('ASM Modification Algorithm');
            s=s+(1.0e-9)*rand(length(s),1); % To differentiate the values of s(k) which may cause problems for finding sorted_indices entries
            sorted_tasks=sort(s);
            sorted_indices=zeros(length(K_set),1);
            for k=K_set
                sorted_indices(k,1)=find(sorted_tasks(k)==s);
            end
            C_tilde_node=zeros(N,1);
            for k=sorted_indices'
                parfor n=1:N
                    C_tilde_node(n,1)=C_node(n,1)-v(n,:)*c+c(k)*v(n,k)-.0001;
                end
                C_tilde_link=zeros(N,N);
                parfor m=1:N
                    for mm=1:N
                        temp=zeros(length(K_set),1);
                        temp1=zeros(max(path),N);
                        temp1(:,:)=I_l2p(m,mm,:,:);
                        for i=K_set
                            if i~=k
                                temp2=zeros(max(path),N);
                                temp2(:,:)=r_aug(:,:,i)';
                                temp(i,1)=sum(sum(temp1.*e(:,:,i)'.*temp2));
                            end
                        end
                        C_tilde_link(m,mm)=C_link(m,mm)-sum(temp);
                    end
                end
                Nodes_e=[];  % The set of nodes which have links with sufficient bandwidth
                feasiblepaths=zeros(N,max(path));
                for n=1:N
                    feasiblepaths_b=[];
                    for b=1:path(n)
                        flag=[];
                        for m=1:N
                            for mm=1:N
                                if I_l2p(m,mm,b,n)==1
                                    if r_convex(k,1)<=C_tilde_link(m,mm)
                                        flag=[flag;1];
                                    else
                                        flag=[flag;0];
                                    end
                                end
                            end
                        end
                        if prod(flag)==1
                            if ismember(n,Nodes_e)==0
                                Nodes_e=[Nodes_e,n] ;
                            end
                            feasiblepaths_b=[feasiblepaths_b,b];
                            %break;
                        end
                    end
                    feasiblepaths(n,1:length(feasiblepaths_b))=feasiblepaths_b;
                end
                Node_feasible=[];
                for n=Nodes_e
                    if C_tilde_node(n,1)>=c(k)-.001
                        Node_feasible=[Node_feasible,n];
                    end
                end
                Tau_exe_prop=nan(N,max(path));
                for n=Node_feasible
                    %   c(k)=C_tilde_node(n,1);
                    for b=setdiff(feasiblepaths(n,:),0)
                        Tau_exe_prop(n,b)=L(k)/C_tilde_node(n,1)+prop(n,b);
                    end
                end
                [n_star,b_star]=find(Tau_exe_prop==min(min(Tau_exe_prop)));
                Tau_prop(k,1)=prop(n_star,b_star);
                v(:,k)=zeros(N,1);
                v(n_star,k)=1;
                e(:,:,k)=zeros(N,max(path));
                e(n_star,b_star,k)=1;
                s_tilde=Tau_ran(k,1)+Tau_prop(k,1)+L(k)/C_tilde_node(n_star,1)-Tau(k,1);
                if s_tilde<0
                    c(k)=L(k)/(Tau(k,1)-Tau_ran(k,1)-Tau_prop(k,1));
                    s(k)=0;
                else
                    s(k)=s_tilde;
                    c(k)=C_tilde_node(n_star,1);
                end
                Tau_exe(k,1)=L(k)/c(k);
            end
            disp('Elastic Variables');
            disp(s);
            disp('Tau_tx+Tau_prop+Tau_exe<=Tau+s');
            disp((Tau_tx+Tau_prop+Tau_exe<=Tau+s)');
            disp('Tau_tx+Tau_prop+Tau_exe-Tau-s');
            disp((Tau_tx+Tau_prop+Tau_exe-Tau-s)');
            
            
            convergence(count_convergence)=sum(s);
            count_convergence=count_convergence+1;
            elastic_stack(:,count)=s;
            disp('The ratio of objective improvement')
            %  disp((obj_feas(count-1)-sum(s))/obj_feas(count-1))
            disp((convergence_check-sum(s))/convergence_check)
            %     if abs(obj_feas(count-1)-sum(s))<10e-1 || abs(obj_feas(count-1)-sum(s))/obj_feas(count-1)<10e-2 || count>=I_max
            if (convergence_check-sum(s))/convergence_check <10e-2|| count>=I_max || convergence_check-sum(s)<10e-1
                break;
            end
            
        end
        if max(s)>=.000001
            %K_set=setdiff(K_set,find(s==max(s)));
            %K_reject=[K_reject; find(s==max(s))];
            K_set=1:length(K_set)-1;
            [~,k_reject]=max(s);
            H(:,:,k_reject)=[];
            c(k_reject)=[];
            p(k_reject)=[];
            p_var(k_reject)=[];
            e(:,:,k_reject)=[];
            v(:,k_reject)=[];
            s(k_reject)=[];
            L(k_reject)=[];
            D(k_reject)=[];
            Tau(k_reject)=[];
            Tau_prop(k_reject)=[];
            Tau_tx(k_reject)=[];
            Tau_exe(k_reject)=[];
            Tau_ran(k_reject)=[];
            r_convex(k_reject)=[];
        else
            break;
        end
    end
    %% TASK OFFLOADING OPTIMIZATION
    cost_computation=sum(v*(c.^3));
    Igamma=zeros(N,max(path));
    for n=1:N
        for b=1:path(n)
            Igamma(n,b)=sum(sum(I_l2p(:,:,b,n).*cost_link))/2;
        end
    end
    Igamma_aug=zeros(N,max(path),length(K_set));
    for k=K_set
        Igamma_aug(:,:,k)=Igamma*r_convex(k,1);
    end
    cost_forwarding=sum(sum(sum(e.*Igamma_aug)));
    ovreall_cost=zeros(I_max,1);
    TO_cost=cost_computation+cost_forwarding;
    count=0;
    while true
        TO_cost=cost_computation+cost_forwarding;
        count=count+1;
        disp('count');
        disp(count);
        disp('computational capacity problem ');
        %% cvx for (19)
        %constraint_bound=Tau+s-Tau_prop-Tau_tx %%%temporary
        cvx_begin
        variable c(length(K_set))
        minimize ( sum(v*pow_p(c,3)) )
        % minimize ( sum(L(K_set).*inv_pos(c)) )
        subject to
        for k=K_set
            c(k) >= (L(k))/(Tau(k)-Tau_prop(k)-Tau_ran(k))
            c(k) >= 0
        end
        for n=1:N
            v(n,:)*c <= C_node(n)
        end
        cvx_end %end of (20)
        disp(cvx_status)
        cost_computation=cvx_optval;
        Tau_exe=L(K_set)./c;% (Formula (7) for execution delay) (Load/capacity
        %     disp('Tau_tx+Tau_prop+Tau_exe<=Tau+s');
        %     disp(Tau_tx+Tau_prop+Tau_exe<=Tau+s);
        %     disp('Tau_exe ');
        %     disp(Tau_exe);
        Exe_cost(count_tau_integer,1)=cost_computation;
        %% MOSEK ... Optimization of binary variables
        Idelta=zeros(N,max(path));
        Igamma=zeros(N,max(path));
        r_aug=zeros(N,max(path),length(K_set));
        for n=1:N
            for b=1:path(n)
                Idelta(n,b)=sum(sum(I_l2p(:,:,b,n).*delta_link))/2;
                Igamma(n,b)=sum(sum(I_l2p(:,:,b,n).*cost_link))/2;
                r_aug(n,b,:)=r_convex;
            end
        end
        
        Igamma_aug=zeros(N,max(path),length(K_set));
        for k=K_set
            Igamma_aug(:,:,k)=Igamma*r_convex(k,1);
        end
        
        I_l2p_aug=zeros(N,N,N,max(path),length(K_set));
        for m=1:N
            for mm=1:N
                temp=zeros(max(path),N);
                
                temp(:,:)=I_l2p(m,mm,:,:);
                for k=K_set
                    I_l2p_aug(m,mm,:,:,k)=temp';
                end
            end
        end
        
        v_aug=zeros(N,max(path),length(K_set));
        for n=1:N
            for b=1:max(path)
                for k=K_set
                    v_aug(n,b,k)=v(n,k);
                end
            end
        end
        %% solving (30) for v,e,gamma
        disp('binary variables problem ');
        cvx_begin
        cvx_solver MOSEK
        variable e_var(N,max(path),length(K_set)) binary
        variable gamma_var(N,max(path),length(K_set)) binary
        variable v_var(N,length(K_set)) binary
        
        minimize sum(sum(sum(e_var.*Igamma_aug)))
        subject to
        for k=K_set
            %%%%%%%%%%%%%%%%%%%%%%%%C1-a DCP error
            sum(sum(e_var(:,:,k).* Idelta)) <= Tau(k)+s(k)-Tau_ran(k)-Tau_exe(k) %C1-a
            sum(sum(e_var(:,:,k)))==1 %C8
            sum(v_var(:,k)) == 1 %C7
            sum(sum(gamma_var(:,:,k))) == 1 %C9
            
            for n=1:N
                for b=1:path(n)
                    gamma_var(n,b,k) <= v_var(n,k)+1-e_var(n,b,k) %C15
                    v_var(n,k)       <= gamma_var(n,b,k)+1-e_var(n,b,k) %C16
                    gamma_var(n,b,k) <= e_var(n,b,k) %C17
                end
                if path(n)<max(path) %making unnecessary variables zero
                    for b=path(n)+1:max(path)
                        e_var(n,b,k)==0
                        gamma_var(n,b,k)==0
                    end
                end
            end
        end
        for n=1:N
            v_var(n,:)*c <= C_node(n) %C2 in (30)
        end
        
        for m=1:N
            for mm=1:N
                temp=zeros(N,max(path),length(K_set));
                temp(:,:,:)=I_l2p_aug(m,mm,:,:,:);
                sum(sum(sum(r_aug.*temp.*e_var))) <= C_link(m,mm) %C3
            end
        end
        cvx_end
        disp(cvx_status)
        cost_forwarding=cvx_optval;
        gamma_var=full(gamma_var); % MATLAB can not put sparse matrices in N>2 Dimensional arrays, which cause error in obtaining Tau_prop(k) 
        e_var=full(e_var);
        v_var=full(v_var);
        
        %% Tau_prop for obtained v_var & e_var
        for k=K_set
            Tau_prop(k)=sum(sum(gamma_var(:,:,k).* Idelta));
        end
        v=v_var;
        e=e_var;
        
        for n=1:N
            for k=K_set
                if abs(v(n,k)-1)<.5
                    v(n,k)=1;
                else
                    v(n,k)=0;
                end
                for b=1:max(path)
                    if abs(e(n,b,k)-1)<.5
                        e(n,b,k)=1;
                    else
                        e(n,b,k)=0;
                    end
                end
            end
        end
        
        %     disp('Tau_prop ');
        %     disp(Tau_prop);
        disp('Tau_tx+Tau_prop+Tau_exe<=Tau+s');
        disp((Tau_tx+Tau_prop+Tau_exe<=Tau+s)');
        disp('Tau_tx+Tau_prop+Tau_exe-Tau-s');
        disp((Tau_tx+Tau_prop+Tau_exe-Tau-s)');
        %% Convergence criterion
        ovreall_cost(count,1)=cost_computation+cost_forwarding; %last iteration's sum(s) (feasibility problem)
        disp('ovreall_cost')
        disp(ovreall_cost(1:count))
            if abs(TO_cost-ovreall_cost(count))<10e-3 || count>=I_max
                break;
            end
    end
    Task_offloading_cost(count_tau_integer,1)=cost_computation+cost_forwarding;
    Acceptance_ratio(count_tau_integer,1)=length(K_set)/K;
end
% save('disjoint_vars.mat','Acceptance_ratio','Task_offloading_cost','count_tau','Tx_power','Acceptance_ratio_power','Exe_cost','H_initial');
%% Outputs
disp('Acceptance Ratio:')
disp(length(K_set)/K)
disp('Radio Transmission Latency:')
disp(Tau_tx)
disp('Propagation Latency:')
disp(Tau_prop)
disp('Execution Latency:')
disp(Tau_exe)