% This function runs the M-SR optimization method

function [sum_rate, complexity,theta_mat] = MSR_center_RIS(K,M,N,Pt_dB,H_nonoise,H_noisy,G1_nonoise,G1_noisy,G2_nonoise,G2_noisy,sigma2,rand_ini_num,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective,theta_mat,noisy_cond)

% distantce matrices
dist_mat_Tx_Rx=zeros(K); 
for i=1:K
    for j=1:K
        dist_mat_Tx_Rx(i,j)=sqrt((Tx_loc(1,i)-Rx_loc(1,j))^2+(Tx_loc(2,i)-Rx_loc(2,j))^2);
    end
end

dist_mat_RIS_Rx=zeros(K);    
for i=1:K
    for j=1:K
        dist_mat_RIS_Rx(i,j)=sqrt((RIS_loc(1,i)-Rx_loc(1,j))^2+(RIS_loc(2,i)-Rx_loc(2,j))^2);
    end
end

dist_mat_Tx_RIS=zeros(K); 
for i=1:K
    for j=1:K
        dist_mat_Tx_RIS(i,j)=sqrt((Tx_loc(1,i)-RIS_loc(1,j))^2+(Tx_loc(2,i)-RIS_loc(2,j))^2);
    end
end

comp_val=0:(2*pi)/N:(2*pi-((2*pi)/N));
All_values=exp(1j*comp_val);

Pt=10^(Pt_dB/10);
Pt=Pt*10^(-3);

sigma2=10^(sigma2/10);
sigma2=sigma2*10^(-3);

c0=10^(c0_dB/10);

if noisy_cond()
    noisy_channel_H=sqrt(c0)*H_nonoise;

    noisy_channel_G1=sqrt(c0)*G1_nonoise;

    noisy_channel_G2=sqrt(c0)*G2_nonoise;
else
    noisy_channel_H=sqrt(c0)*H_noisy;

    noisy_channel_G1=sqrt(c0)*G1_noisy;

    noisy_channel_G2=sqrt(c0)*G2_noisy;
end

    eval_mat=0;
    sum_max_opt=0;
    
    for rand_idx=1:rand_ini_num

        if sum(sum(theta_mat))~=0
            theta_ini=theta_mat(rand_idx,:);
        else
            ar = randi([1,N],1,M);
            theta_ini = All_values(ar);
            theta_mat(rand_idx,:)=theta_ini; 
        end
        
        for i1=1:M
            sum_max=0;
            for item=1:N
                eval_mat=eval_mat+1;
                alpha_L=All_values(item);
                if sum(theta_ini==0)~=M
                    Theta_mat=diag(theta_ini);
                    SINR=zeros(1,K);
                    for User=1:K
                        H_noisy=transpose(squeeze(noisy_channel_H(1,:,User)));
                        H_noisy=H_noisy./sqrt(dist_mat_Tx_Rx(:,User).^alpha_d);

                        G1_noisy=squeeze(noisy_channel_G1(1,:,:));
                        for tx=1:K
                        G1_noisy(tx,:)=G1_noisy(tx,:)/sqrt(dist_mat_Tx_RIS(tx,tx)^alpha_d1);
                        end

                        G2_noisy=squeeze(noisy_channel_G2(1,:,User));
                        G2_noisy=G2_noisy/sqrt(dist_mat_RIS_Rx(User,User).^alpha_d2);
                        G2_noisy=transpose(G2_noisy);

                        Desired_signal=H_noisy(User)+G1_noisy(User,:)*Theta_mat*G2_noisy;
                        Desired_signal=Desired_signal+(alpha_L-theta_ini(1,i1))*G1_noisy(User,:)*G2_noisy;

                        H_intrf_noisy=H_noisy;
                        H_intrf_noisy(User)=[];

                        G1_intrf_noisy=G1_noisy;
                        G1_intrf_noisy(User,:)=[];

                        Intrf_signal=H_intrf_noisy+G1_intrf_noisy*Theta_mat*G2_noisy;
                        Intrf_signal=Intrf_signal+(alpha_L-theta_ini(1,i1))*G1_intrf_noisy*G2_noisy;

                        Desired_pow=Pt*(abs(Desired_signal)^2);
                        Intrf_pow=Pt*sum((abs(Intrf_signal).^2));
                        SINR(User)=Desired_pow/(sigma2+Intrf_pow);
                    end

                    if strcmp(opt_objective,'max_min')
                        if min(SINR)>sum_max
                            sum_max=min(SINR);
                            SINR_opt=SINR;
                            alpha_opt=alpha_L;
                        end
                    end

                    if strcmp(opt_objective,'sum')
                        if sum(log2(1+SINR))>sum_max
                            sum_max=sum(log2(1+SINR));
                            SINR_opt=log2(1+SINR);
                            alpha_opt=alpha_L;
                        end
                    end
                end
            end
            theta_ini(1,i1)=alpha_opt;
        end
        if sum_max>sum_max_opt
            sum_max_opt=sum_max;
            theta_opt=theta_ini;
        end
    end
    correct_channel_H=sqrt(c0)*H_nonoise;
    correct_channel_G1=sqrt(c0)*G1_nonoise;
    correct_channel_G2=sqrt(c0)*G2_nonoise;
            
    Theta_mat=diag(theta_opt);   
    for User=1:K                
        H_correct=transpose(squeeze(correct_channel_H(1,:,User)));
        H_correct=H_correct./sqrt(dist_mat_Tx_Rx(:,User).^alpha_d);

        G1_correct=squeeze(correct_channel_G1(1,:,:));
        for tx=1:K
        G1_correct(tx,:)=G1_correct(tx,:)/sqrt(dist_mat_Tx_RIS(tx,tx)^alpha_d1);
        end
        
        G2_correct=squeeze(correct_channel_G2(1,:,User));
        G2_correct=G2_correct/sqrt(dist_mat_RIS_Rx(User,User).^alpha_d2);
        G2_correct=transpose(G2_correct);

        Desired_signal=H_correct(User)+G1_correct(User,:)*Theta_mat*G2_correct;

        H_intrf_correct=H_correct;
        H_intrf_correct(User)=[];

        G1_intrf_correct=G1_correct;
        G1_intrf_correct(User,:)=[];

        Intrf_signal=H_intrf_correct+G1_intrf_correct*Theta_mat*G2_correct;

        Desired_pow=Pt*(abs(Desired_signal)^2);
        Intrf_pow=Pt*sum((abs(Intrf_signal).^2));
        
        SINR_opt(User)=Desired_pow/(sigma2+Intrf_pow);
    end

    
    if strcmp(opt_objective,'max_min')
        min_SINR=min(SINR_opt);
    end

    if strcmp(opt_objective,'max_min')
        sum_rate=mean(log2(1+min_SINR));
    else
        sum_rate=mean(sum(log2(1+SINR_opt),2));
    end
    complexity=mean(eval_mat);   
end












    