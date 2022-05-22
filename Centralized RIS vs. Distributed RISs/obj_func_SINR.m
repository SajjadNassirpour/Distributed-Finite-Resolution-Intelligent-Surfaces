% -- Objective function of the optimzaion problem -- %

function snir_eval=obj_func_SINR(K,S,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective)

    M=size(S,2);
    
    SINR=zeros(1,K);
    zero_case=0;
    for User_idx=1:K
        if sum(S(User_idx,:)==0)==M
            zero_case=1;
        end
    end
    theta=S;
    if zero_case==0                    
        for User=1:K
            Theta_mat=diag(theta(User,:));
            H_correct=squeeze(H(User,:));
            H_correct=H_correct./sqrt(dist_mat_Tx_Rx(User,:).^alpha_d);

            G1_correct=squeeze(G1(User,:));
            G1_correct=G1_correct./sqrt(dist_mat_Tx_RIS(User,User)^alpha_d1);

            G2_correct=squeeze(G2(User,:,:));
            for k2=1:K
                G2_correct(:,k2)=G2_correct(:,k2)./sqrt(dist_mat_RIS_Rx(User,k2).^alpha_d2);
            end

            Desired_signal=H_correct(User)+G1_correct*Theta_mat*squeeze(G2_correct(:,User));


            H_intrf_correct=H_correct;
            H_intrf_correct(User)=[];

            G2_intrf_correct=G2_correct;
            G2_intrf_correct(:,User)=[];

            Intrf_signal=H_intrf_correct+G1_correct*Theta_mat*G2_intrf_correct;

            Desired_pow=Pt*(abs(Desired_signal)^2);


            Intrf_pow=Pt*sum((abs(Intrf_signal).^2));
            SINR(User)=-Desired_pow/(sigma2+Intrf_pow);
        end
        if strcmp(opt_objective,'max_min')
        snir_eval=max(SINR);
        end
        
        if strcmp(opt_objective,'sum')
        snir_eval=sum(SINR);
        end
            
    else
        snir_eval=0;
    end
end