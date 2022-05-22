% -- Local search algorithm -- %

function [LM_string,eval]=LS_algo_center_RIS(K,N,ini_string,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,mode,r,opt_objective)
    
    comp_val=0:(2*pi)/N:(2*pi-((2*pi)/N));
    All_values=exp(1j*comp_val);

    eval=0;    
    M=length(ini_string);
    S_star=ini_string;
    
    % mode=1: Local search based on SINR
    % mode=2: Local search based on Filled function
    if mode==1
        score_max=obj_func_SINR_center_RIS(K,ini_string,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
    else
        score_max=Filled_func_center_RIS(K,S_star,S_star,r,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
    end
    
    % First round of local searching process
    theta_ini = ini_string;
    LM_string=ini_string;
    for i1=1:M
        for item=1:N
            if theta_ini(1,i1)~=All_values(item)
                theta=theta_ini;
                theta(1,i1)=All_values(item);
                if sum(theta==0)~=M
                    eval=eval+1;
                    if mode==1
                        score=obj_func_SINR_center_RIS(K,theta,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
                    else
                        score=Filled_func_center_RIS(K,S_star,theta,r,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
                    end
                    if score<=score_max 
                        score_max=score;
                        LM_string=theta;
                    end
                end
            end
        end
    end
    
    % Run the while loop until theta_ini will be the local optimum
    ok=0;
    while ok==0 
        theta_ini = LM_string;
        for i1=1:M
            for item=1:N
                if theta_ini(1,i1)~=All_values(item)
                    theta=theta_ini;
                    theta(1,i1)=All_values(item);
                    if sum(theta==0)~=M
                        eval=eval+1;
                        if mode==1
                            score=obj_func_SINR_center_RIS(K,theta,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
                        else
                            score=Filled_func_center_RIS(K,S_star,theta,r,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
                        end
                        if score<=score_max 
                            score_max=score;
                            LM_string=theta;
                        end
                    end
                end
            end
        end
        if sum(abs(LM_string-theta_ini)==0)==M  
            ok=1;
        end
    end

end




