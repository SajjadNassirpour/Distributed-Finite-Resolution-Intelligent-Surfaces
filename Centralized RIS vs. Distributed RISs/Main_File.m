clc
clear
close all

K=4;
N=4;
M_budget=8;

iters=20;

Pt_dB=20;
sigma2=-80;

% pathloss exponents
alpha_d=3.5;
alpha_d1=2 ;
alpha_d2=2.1;

c0_dB=-30;

% Objective function of the optimization problem (centralized RIS)
% sum:      sum-rate
% max_min:  max-min rate
opt_objective='sum';


% Optimiztion parameters
ini_chance=1;
LN_check=10;
mu=10^(-1);
r=1e1;
divide_rate=10;

sum_rate_mat_perfect=zeros(1,iters);
sum_rate_mat_distrb_RIS=zeros(1,iters);
sum_rate_mat_center_RIS=zeros(iters,3);

for iter=1:iters
iter
% ---- Txi locations ---- %
Tx_x_pos_mat=[50  0    100   50];
Tx_y_pos_mat=[0   50   50   100];
Tx_loc=[Tx_x_pos_mat;Tx_y_pos_mat];

% ---- Rxi locations ---- %
a=0; 
b=100;
ini_Rx_x = (b-a).*rand(1,K) + a;
ini_Rx_y = (b-a).*rand(1,K) + a;

for rx_idx=1:K
    min_dist=inf;
    for ini_idx=1:length(ini_Rx_x)
        dist=sqrt((Tx_x_pos_mat(rx_idx)-ini_Rx_x(ini_idx))^2+(Tx_y_pos_mat(rx_idx)-ini_Rx_y(ini_idx))^2);
        if dist<min_dist
            min_dist=dist;
            ini_idx_opt=ini_idx;
            Rx_x_pos_mat(rx_idx) = ini_Rx_x(ini_idx);
            Rx_y_pos_mat(rx_idx) = ini_Rx_y(ini_idx);
        end
    end
    ini_Rx_x(ini_idx_opt)=[];
    ini_Rx_y(ini_idx_opt)=[];
end
Rx_loc=[Rx_x_pos_mat;Rx_y_pos_mat];


% Direct path channels
H_nonoise=Rayleigh_direct(K);


% sum rate No-RIS K parallel channels
sum_rate_mat_perfect(iter)=Perfect_alignment(K,Pt_dB,H_nonoise,sigma2,Tx_loc,Rx_loc,alpha_d,c0_dB);


%% Distributed RISs

% ---- RISi locations ---- %
RIS_x_pos_mat=[47  3    97   53];
RIS_y_pos_mat=[4   54   46   96];
RIS_loc=[RIS_x_pos_mat;RIS_y_pos_mat];

% ---- Reflected path channels ---- %
M=M_budget/K;
G1_nonoise=zeros(K,1,M);
G2_nonoise=zeros(K,1,M,K);
G1_far_nonoise=zeros(K,1,K,M);

for k=1:K
    d_Tx_RIS=sqrt((Tx_loc(1,k)-RIS_loc(1,k))^2+(Tx_loc(2,k)-RIS_loc(2,k))^2); 
    G1_nonoise(k,:,:)=LoS_h_multi_RIS(M,d_Tx_RIS);
    G2_nonoise(k,:,:,:)=Rayleigh_g_multi_RIS(K,M);
    G1_far_nonoise(k,:,:,:)=Rayleigh_h_far_multi_RIS(K,M);
end


% Our approach (distributed RIS)
[sum_rate_mat_distrb_RIS(iter), sum_rate_without_Tx, BF_our_ratio]=Our_approach_distrb_RISs(K,M,N,Pt_dB,H_nonoise,G1_nonoise,G1_far_nonoise,G2_nonoise,sigma2,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,mu,r,divide_rate,ini_chance,opt_objective,LN_check);



%% Centralized RIS
for loc_case=1:3

    % ---- RIS locations ---- %
    if loc_case==1
        RIS_x_pos_mat=[50  50   50  50];
        RIS_y_pos_mat=[50  50   50  50];
    elseif loc_case==2
        RIS_x_pos_mat=[47  47   47  47];
        RIS_y_pos_mat=[4   4    4   4];
    elseif loc_case==3
        a=0; 
        b=100;
        RIS_x_pos_mat = ((b-a).*rand(1,1) + a)*ones(1,K);
        RIS_y_pos_mat = ((b-a).*rand(1,1) + a)*ones(1,K);
    end
    RIS_loc=[RIS_x_pos_mat;RIS_y_pos_mat];

    % ---- Reflected path channels ---- %
    M=M_budget;
    G1_nonoise=zeros(1,K,M);
    for k=1:K
       d_Tx_RIS=sqrt((Tx_loc(1,k)-RIS_loc(1,k))^2+(Tx_loc(2,k)-RIS_loc(2,k))^2); 
       G1_nonoise(:,k,:)=LoS_h_asym(M,d_Tx_RIS);
    end
    G2_nonoise=Rayleigh_g_asym(K,M);

    % Our approach (centralized RIS)
    sum_rate_mat_center_RIS(iter,loc_case)=Our_approach_center_RIS(K,M,N,Pt_dB,H_nonoise,G1_nonoise,G2_nonoise,sigma2,mu,r,divide_rate,ini_chance,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective,LN_check);
end

end

% Average rates

Mean_centralized_RIS=mean(sum_rate_mat_center_RIS)
Mean_distributed_RISs=mean(sum_rate_mat_distrb_RIS)
Mean_parallel_channels=mean(sum_rate_mat_perfect)


% Outage capacity
Outage_capacity_distributed_RISs=outage(sum_rate_mat_distrb_RIS','dist')

Outage_capacity_centralized_RIS=outage(sum_rate_mat_center_RIS(:,2),'center')


