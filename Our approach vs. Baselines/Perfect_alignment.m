% This function obtains the results of No-RIS K parallel channels

function sum_rate_perfect=Perfect_alignment(K,Pt_dB,H_nonoise,sigma2,Tx_loc,Rx_loc,alpha,c0_dB,opt_objective)

    Pt=10^(Pt_dB/10);
    Pt=Pt*10^(-3);
    
    sigma2=10^(sigma2/10);
    sigma2=sigma2*10^(-3);
    
    c0=10^(c0_dB/10);
    
    dist_mat_Tx_Rx=zeros(K);
    for i=1:K
        for j=1:K
            dist_mat_Tx_Rx(i,j)=sqrt((Tx_loc(1,i)-Rx_loc(1,j))^2+(Tx_loc(2,i)-Rx_loc(2,j))^2);
        end
    end

    correct_channel=H_nonoise;
    SINR=zeros(1,K);
    for User=1:K
    H_correct=squeeze(correct_channel(1,User,:));

    SINR(User)=Pt*(abs(H_correct(User))^2)*c0/(dist_mat_Tx_Rx(User,User)^alpha*sigma2);
    end

    if strcmp(opt_objective,'max_min')
        min_SINR=min(SINR);
    end
    
    if strcmp(opt_objective,'max_min')
        sum_rate_perfect=mean(log2(1+min_SINR));
    else
        sum_rate_perfect=mean(sum(log2(1+SINR),2));
    end

end









    