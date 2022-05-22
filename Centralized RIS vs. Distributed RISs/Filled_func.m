% -- The sigmoid filled function -- %

function F_r=Filled_func(K,S_star,S,r,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective)

    f_S=obj_func_SINR(K,S,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
    f_S_star=obj_func_SINR(K,S_star,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
    t=f_S-f_S_star;
    
    if t>=0
        g=1;
    elseif t<=-r
        g=t+r;
    else
        g=1/(1+exp((-6/r)*(t+r/2)));
    end
    
    %----f_r(t)-----%
    if t>=0
        f=1;
    elseif t<=-r
        f=t+r;
    else
        f=1/(1+exp((-6/r)*(t+r/2)));
    end
    
    %----Filled function-----%
    
    if t<=-r
        cond=0;
    else
        cond=1;
    end
    
    F_r=(1/(cond*(norm(S-S_star)^2)+1))*g+f;
end