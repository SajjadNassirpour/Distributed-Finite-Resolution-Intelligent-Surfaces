% -- Objective function of the optimzaion problem -- %

function snir_eval=obj_func_SINR_final_center_RIS(K,S,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0)
    M=length(S);
    Theta_mat=diag(S);
    SINR=zeros(1,K);
    if sum(S==0)~=M
        for User=1:K
            H_correct=squeeze(H(:,User));
            H_correct=H_correct./sqrt(dist_mat_Tx_Rx(:,User).^alpha_d);

            G1_correct=squeeze(G1(:,:));
            for tx=1:K
            G1_correct(tx,:)=G1_correct(tx,:)/sqrt(dist_mat_Tx_RIS(tx,tx)^alpha_d1);
            end

            G2_correct=squeeze(G2(:,User));
            G2_correct=G2_correct/sqrt(dist_mat_RIS_Rx(User,User).^alpha_d2);

            Desired_signal=H_correct(User)+G1_correct(User,:)*Theta_mat*G2_correct;

            H_intrf_correct=H_correct;
            H_intrf_correct(User)=[];

            G1_intrf_correct=G1_correct;
            G1_intrf_correct(User,:)=[];

            Intrf_signal=H_intrf_correct+G1_intrf_correct*Theta_mat*G2_correct;

            Desired_pow=Pt*(abs(Desired_signal)^2);
            Intrf_pow=Pt*sum((abs(Intrf_signal).^2));

            SINR(User)=Desired_pow/(sigma2+Intrf_pow);
        end
        snir_eval=SINR;
    else
        snir_eval=0;
    end
end