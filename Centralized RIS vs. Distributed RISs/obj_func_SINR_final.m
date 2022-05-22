% -- Objective function of the optimzaion problem -- %

function [SINR_with_Tx, SINR_without_Tx] = obj_func_SINR_final(K,S,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G1_far,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0)

    SINR_with_Tx=zeros(1,K);
    SINR_without_Tx=zeros(1,K);

    
    theta_opt=S;
                      
    for User=1:K
        H_correct=squeeze(H(:,User));
        H_correct=H_correct./sqrt(dist_mat_Tx_Rx(:,User).^alpha_d);

        G1_correct=squeeze(G1(:,:));
        
        for tx=1:K
            G1_correct(tx,:)=G1_correct(tx,:)/sqrt(dist_mat_Tx_RIS(tx,tx)^alpha_d1);
        end
        
        RIS_mat=1:K;
        RIS_mat(User)=[];
        
        G1_far_desired_correct=squeeze(G1_far(RIS_mat,User,:));
        
        for k2=1:K-1
            G1_far_desired_correct(k2,:)=G1_far_desired_correct(k2,:)./sqrt(dist_mat_Tx_RIS(User,RIS_mat(k2)).^alpha_d1);
        end

        G2_correct=squeeze(G2(:,:,User));
        for k2=1:K
            G2_correct(k2,:)=G2_correct(k2,:)./sqrt(dist_mat_RIS_Rx(k2,User).^alpha_d2);
        end

        
        theta_desired=theta_opt(User,:);

        Theta_mat=diag(theta_desired);
        Desired_signal=H_correct(User)+squeeze(G1_correct(User,:))*Theta_mat*transpose(squeeze(G2_correct(User,:)));
        
        Desired_signal_without_Tx=Desired_signal;

        for kk=1:K-1
            theta_desired1=theta_opt(RIS_mat(kk),:);
            Theta_mat1=diag(theta_desired1);
            Desired_signal=Desired_signal+squeeze(G1_far_desired_correct(kk,:))*Theta_mat1*transpose(squeeze(G2_correct(RIS_mat(kk),:)));
        end

        intf_idx=1:K;
        intf_idx(User)=[];

        Intrf_signal=zeros(K-1,1);

        for kk=1:K-1
            user_intrf_idx=intf_idx(kk);
            theta_intrf=theta_opt(user_intrf_idx,:);

            Theta_mat=diag(theta_intrf);
            Intrf_signal(kk)=H_correct(user_intrf_idx)+squeeze(G1_correct(user_intrf_idx,:))*Theta_mat*transpose(squeeze(G2_correct(user_intrf_idx,:)));
            
            Intrf_signal_without_Tx(kk)=Intrf_signal(kk);
            
            for RIS_idx=1:K-1 
                RIS_id=1:K;
                RIS_id(user_intrf_idx)=[];
                G1_far_intrf_correct=squeeze(G1_far(RIS_id,user_intrf_idx,:));
        
                for k2=1:K-1
                    G1_far_intrf_correct(k2,:)=G1_far_intrf_correct(k2,:)./sqrt(dist_mat_Tx_RIS(user_intrf_idx,RIS_id(k2)).^alpha_d1);
                    
                    theta_intrf1=theta_opt(RIS_id(RIS_idx),:);
                    Theta_mat2=diag(theta_intrf1);
                    Intrf_signal(kk)=Intrf_signal(kk)+squeeze(G1_far_intrf_correct(k2,:))*Theta_mat2*transpose(squeeze(G2_correct(RIS_id(RIS_idx),:)));
                end
                
            end
        
        end
        
        %--- with Tx
            
        Desired_pow=Pt*(abs(Desired_signal)^2);

        Intrf_pow=Pt*sum((abs(Intrf_signal).^2));

        SINR_with_Tx(User)=Desired_pow/(sigma2+Intrf_pow);
        
        
        %--- without Tx
        Desired_pow_without_Tx=Pt*(abs(Desired_signal_without_Tx)^2);

        Intrf_pow_without_Tx=Pt*sum((abs(Intrf_signal_without_Tx).^2));

        SINR_without_Tx(User)=Desired_pow_without_Tx/(sigma2+Intrf_pow_without_Tx);
        
    end
    
end