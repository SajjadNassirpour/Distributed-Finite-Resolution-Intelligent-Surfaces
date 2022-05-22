clc
clear
close all

K=4;
N=4;
M=8;

% A list of the optimiztion methods
% 1: Brute force
% 2: Our method
% 3: SR
% 4: M-SR
optimizer_list=[1];

% number of iterations
iters=3;


Pi_dB=20;
sigma2=-80;

% pathloss exponents
alpha_d=3.5;
alpha_d1=2 ;
alpha_d2=2.1;

c0_dB=-30; 
beta=2;     % Rician factor
est_SNR=20; % SNR of the noisy channels

% Noise status
% 1: Both noisy and noiseless channels
% 2: Noiseless channels
% 3: Noisy channels
noise_stat=1;


% Objective function of the optimization problem
% sum:      sum-rate
% max_min:  max-min rate
opt_objective='sum';

for iter=1:iters
    iter

% ---- Txi locations ---- %
Tx_x_pos_mat=[0  0   0   0];
Tx_y_pos_mat=[0  0   0   0];
Tx_loc=[Tx_x_pos_mat;Tx_y_pos_mat];

% ---- Rxi locations ---- %
Rx_x_pos_mat=[50 50 50 50];
Rx_y_pos_mat=[0  0  0  0];
Rx_loc=[Rx_x_pos_mat;Rx_y_pos_mat];

% ---- RISi locations ---- %
RIS_x_pos_mat=[3  3  3  3];
RIS_y_pos_mat=[4  4  4  4];
RIS_loc=[RIS_x_pos_mat;RIS_y_pos_mat];


% ---- Direct path channels ---- %
[hd_noisy,hd_nonoise]=Rayleigh_direct(K,est_SNR);


% No-RIS K parallel channels
rate_mat_perfect(iter,:)=Perfect_alignment(K,Pi_dB,hd_nonoise,sigma2,Tx_loc,Rx_loc,alpha_d,c0_dB,opt_objective);


% ---- h and g channels
    h_nonoise=zeros(1,K,M);
    h_noisy=zeros(1,K,M);
    for k=1:K
       d_Tx_RIS=sqrt((Tx_loc(1,k)-RIS_loc(1,k))^2+(Tx_loc(2,k)-RIS_loc(2,k))^2); 
       [h_noisy(:,k,:),h_nonoise(:,k,:)]=Rician_h_center_RIS(M,d_Tx_RIS,est_SNR,beta);
    end

    [g_noisy,g_nonoise]=Rayleigh_g_center_RIS(K,M,est_SNR);

    
    % Optimization parameters
    
    % SR method
    rand_ini_num=100;
    threshold=1e-4;
    
    % Our approach
    ini_chance=1;
    LN_check=10;
    mu=10^(-1);
    r=1e1;
    divide_rate=10;
 
    % Initial RIS configutaion
    theta_mat_Zhang=zeros(rand_ini_num,M);
    theta_mat_Zaid=zeros(rand_ini_num,M);
    theta_mat_GS=zeros(1,M);
    
%% Approaches

    % noisy channels
    if noise_stat==1 || noise_stat==3
        noisy_cond=false;

        if sum(optimizer_list==1)==1
            sum_rate_mat_BF_noisy(iter)=BF_search_center_RIS(K,M,N,Pi_dB,hd_nonoise,hd_noisy,h_nonoise,h_noisy,g_nonoise,g_noisy,sigma2,noisy_cond,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective);
        end

        if sum(optimizer_list==2)==1
            [sum_rate_mat_center_RIS_GS_noisy(iter), Eval_mat_BF_GS(iter,:), SINR_our,theta_mat_GS]=Our_approach_center_RIS(K,M,N,Pi_dB,hd_nonoise,hd_noisy,h_nonoise,h_noisy,g_nonoise,g_noisy,sigma2,mu,r,divide_rate,noisy_cond,ini_chance,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective,LN_check,theta_mat_GS);
        end

        if sum(optimizer_list==3)==1
            [sum_rate_mat_center_RIS_Zhang_noisy(iter), Eval_mat_Zhang(iter), theta_mat_Zhang]=SR_center_RIS(K,M,N,Pi_dB,hd_nonoise,hd_noisy,h_nonoise,h_noisy,g_nonoise,g_noisy,sigma2,rand_ini_num,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective,threshold,theta_mat_Zhang,noisy_cond);
        end

        if sum(optimizer_list==4)==1
            [sum_rate_mat_center_RIS_Zaid_noisy(iter), Eval_mat_Zaid(iter), theta_mat_Zaid]=MSR_center_RIS(K,M,N,Pi_dB,hd_nonoise,hd_noisy,h_nonoise,h_noisy,g_nonoise,g_noisy,sigma2,rand_ini_num,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective,theta_mat_Zaid,noisy_cond);
        end
    end
    
    
    % noiseless channels
    if noise_stat==1 || noise_stat==2
        noisy_cond=true;

        if sum(optimizer_list==1)==1
            sum_rate_mat_BF(iter)=BF_search_center_RIS(K,M,N,Pi_dB,hd_nonoise,hd_noisy,h_nonoise,h_noisy,g_nonoise,g_noisy,sigma2,noisy_cond,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective);
        end

        if sum(optimizer_list==2)==1
            [sum_rate_mat_center_RIS_GS(iter), BF_our_GS, SINR_our,theta_mat_GS]=Our_approach_center_RIS(K,M,N,Pi_dB,hd_nonoise,hd_noisy,h_nonoise,h_noisy,g_nonoise,g_noisy,sigma2,mu,r,divide_rate,noisy_cond,ini_chance,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective,LN_check,theta_mat_GS);
        end

        if sum(optimizer_list==3)==1
            [sum_rate_mat_center_RIS_Zhang(iter), eval_Zhang, theta_mat_Zhang]=SR_center_RIS(K,M,N,Pi_dB,hd_nonoise,hd_noisy,h_nonoise,h_noisy,g_nonoise,g_noisy,sigma2,rand_ini_num,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective,threshold,theta_mat_Zhang,noisy_cond);
        end

        if sum(optimizer_list==4)==1
            [sum_rate_mat_center_RIS_Zaid(iter), eval_Zaid, theta_mat_Zaid]=MSR_center_RIS(K,M,N,Pi_dB,hd_nonoise,hd_noisy,h_nonoise,h_noisy,g_nonoise,g_noisy,sigma2,rand_ini_num,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective,theta_mat_Zaid,noisy_cond);
        end
    end

end

% Average rate of No-RIS K parallel channel 
if strcmp(opt_objective,'sum')
    Mean_sum_rate_perfect=mean(rate_mat_perfect)
else
    Mean_min_rate_perfect=mean(rate_mat_perfect)
end


% ------- Average results of the noiseless channels ------- %
if noise_stat==1 || noise_stat==2

    if sum(optimizer_list==1)==1
        % Average rate
        if strcmp(opt_objective,'sum')
            Mean_sum_rate_BF_noiseless = mean(sum_rate_mat_BF)
        else
            Mean_min_rate_BF_noiseless = mean(sum_rate_mat_BF) 
        end
    end

    if sum(optimizer_list==2)==1
        % Average rate 
        if strcmp(opt_objective,'sum')
            Mean_sum_rate_our_noiseless = mean(sum_rate_mat_center_RIS_GS)
        else
            Mean_min_rate_our_noiseless = mean(sum_rate_mat_center_RIS_GS)
        end
    end

    if sum(optimizer_list==3)==1
        % Average rate 
        if strcmp(opt_objective,'sum')
            Mean_sum_rate_SR_noiseless=mean(sum_rate_mat_center_RIS_Zhang)
        else
            Mean_min_rate_SR_noiseless=mean(sum_rate_mat_center_RIS_Zhang)
        end
    end

    if sum(optimizer_list==4)==1
        % Average rate
        if strcmp(opt_objective,'sum')
            Mean_sum_rate_MSR_noiseless=mean(sum_rate_mat_center_RIS_Zaid)
        else
            Mean_min_rate_MSR_noiseless=mean(sum_rate_mat_center_RIS_Zaid)
        end
    end
end



% ------- Average results of the noisy channels ------- %
if noise_stat==1 || noise_stat==3

    if sum(optimizer_list==1)==1
        % Average rate
        if strcmp(opt_objective,'sum')
            Mean_sum_rate_BF_noisy = mean(sum_rate_mat_BF_noisy)
        else
            Mean_min_rate_BF_noisy = mean(sum_rate_mat_BF_noisy) 
        end
    end

    if sum(optimizer_list==2)==1
        % Average rate 
        if strcmp(opt_objective,'sum')
            Mean_sum_rate_our_noisy = mean(sum_rate_mat_center_RIS_GS_noisy)
        else
            Mean_min_rate_our_noisy = mean(sum_rate_mat_center_RIS_GS_noisy)
        end

        % Average complexity
        mean_complexity_BF = mean(Eval_mat_BF_GS(:,1))
        mean_complexity_our = mean(Eval_mat_BF_GS(:,2))
    end

    if sum(optimizer_list==3)==1
        % Average rate 
        if strcmp(opt_objective,'sum')
            Mean_sum_rate_SR_noisy=mean(sum_rate_mat_center_RIS_Zhang_noisy)
        else
            Mean_min_rate_SR_noisy=mean(sum_rate_mat_center_RIS_Zhang_noisy)
        end
        
        % Average complexity 
        mean_complexity_SR=mean(Eval_mat_Zhang)
    end

    if sum(optimizer_list==4)==1
        % Average rate
        if strcmp(opt_objective,'sum')
            Mean_sum_rate_MSR_noisy=mean(sum_rate_mat_center_RIS_Zaid_noisy)
        else
            Mean_min_rate_MSR_noisy=mean(sum_rate_mat_center_RIS_Zaid_noisy)
        end
        % Average complexity
        mean_complexity_MSR=mean(Eval_mat_Zaid)
    end
end







