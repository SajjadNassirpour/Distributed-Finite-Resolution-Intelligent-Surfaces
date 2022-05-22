% -- Our optimization appraoch based on sigmoid filled function -- %

function [our_final_rate, BF_our_ratio, SINR_opt, theta_mat]=Our_approach_center_RIS(K,M,N,Pt_dB,H_nonoise,H_noisy,G1_nonoise,G1_noisy,G2_nonoise,G2_noisy,sigma2,mu,r,divide_rate,noisy_cond,ini_chance,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,opt_objective,LN_check,theta_mat)
    

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

    Pt=10^(Pt_dB/10);
    Pt=Pt*10^(-3);
    
    sigma2=10^(sigma2/10);
    sigma2=sigma2*10^(-3);
    
    c0=10^(c0_dB/10);
       
    % Possible values for each element at RIS
    comp_val=0:(2*pi)/N:(2*pi-((2*pi)/N));
    All_values=exp(1j*comp_val);
    
    
    % Optimization Parameters
    input_mu=mu;
    input_r=r;
    % \tau
    input_LN_check=LN_check;
    
    % Number of evaluations
    Eval=1;

    % The channels
    if (noisy_cond)
       H=sqrt(c0)*squeeze(H_nonoise(1,:,:));
       G1=sqrt(c0)*squeeze(G1_nonoise(1,:,:));
       G2=sqrt(c0)*squeeze(G2_nonoise(1,:,:));
    else
       H=sqrt(c0)*squeeze(H_noisy(1,:,:)); 
       G1=sqrt(c0)*squeeze(G1_noisy(1,:,:));
       G2=sqrt(c0)*squeeze(G2_noisy(1,:,:));
    end

    % Generate random initial switch configuration
    min_SINR=0;
    for ini_p=1:ini_chance

    if sum(sum(theta_mat))~=0
        theta_ini=theta_mat(1,:);
    else
        ar = randi([1,N],1,M);
        theta_ini = All_values(ar);
        theta_mat(1,:)=theta_ini;
    end


        %G1_ini = theta_ini;

        %===== Filled-function approach =====%
        mu=input_mu;
        r=input_r;

        r0=r;
        resume=1;

        step=0;
        ii=1;
        kk=0;
        S_star_mat=[];
        S_final=[];

        while resume==1
            step=step+1;
            jump=0;   
            if step==LN_check
                LN_check=LN_check+input_LN_check;
                mode=1;
                [S_star,eval]=LS_algo_center_RIS(K,N,theta_ini,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,mode,r,opt_objective);
            else
                eval=0;
                S_star=theta_ini;
            end
            Eval=Eval+eval;
            S_star_mat(kk+1,:)=S_star;

            if kk==0
                min_f=obj_func_SINR_center_RIS(K,S_star_mat(kk+1,:),H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
                S_final=S_star_mat(kk+1,:);
            end

            if kk>0
               if obj_func_SINR_center_RIS(K,S_star_mat(kk+1,:),H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective)>= min_f
                   jump=1;
               else
                   min_f=obj_func_SINR_center_RIS(K,S_star_mat(kk+1,:),H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
                   S_final=S_star_mat(kk+1,:);
                   r=r0;
               end
            end

            if jump==0
                ii=1;
                item=0;
                S_star=S_star_mat(kk+1,:);
                mode=2;
                [S_bar,eval]=LS_algo_center_RIS(K,N,S_star,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,mode,r,opt_objective);
                Eval=Eval+eval;
                theta_ini = S_bar;
                kk=kk+1;
            end

            if ii>=M+1
                if r>=mu
                    r=r/divide_rate;
                    ii=1;
                    item=0;
                    kk=kk-1;
                    S_star=S_star_mat(kk+1,:);
                    mode=2;
                    [S_bar,eval]=LS_algo_center_RIS(K,N,S_star,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,mode,r,opt_objective);
                    Eval=Eval+eval;
                    theta_ini = S_bar;
                    kk=kk+1;
                    jump=0;
                else
                    resume=0;
                end
            end    

            if ii<M+1 && jump==1
                kk=kk-1;
                S_star=S_star_mat(kk+1,:);
                kk=kk+1;
                extra=0;
                while extra==0 && ii<M+1
                    item=item+1;
                    if S_star(1,ii)~=All_values(item)
                        S_star_j0=S_star;
                        S_star_j0(1,ii)=All_values(item);
                        if sum(S_star_j0==0)~=M
                            extra=1;
                            mode=2;
                            [S_bar,eval]=LS_algo_center_RIS(K,N,S_star_j0,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,mode,r,opt_objective);
                            Eval=Eval+eval;
                            theta_ini = S_bar;
                        end
                    end
                    if item==N
                        ii=ii+1;
                        item=0;
                    end
                end
            end
        end
    end

    % Final evaluation based on optimal switch configuration
    H_correct=sqrt(c0)*squeeze(H_nonoise(1,:,:)); 
    G1_correct=sqrt(c0)*squeeze(G1_nonoise(1,:,:));
    G2_correct=sqrt(c0)*squeeze(G2_nonoise(1,:,:));
    SINR_opt=obj_func_SINR_final_center_RIS(K,S_final,H_correct,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1_correct,G2_correct,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0);


    % Number of evaluations
    if strcmp(opt_objective,'max_min')
        min_SINR(1,:)=min(SINR_opt);
    end
    
    BF_our_ratio=[N^M Eval];
    if strcmp(opt_objective,'max_min')
        our_final_rate=mean(log2(1+min_SINR));
    else
        our_final_rate=mean(sum(log2(1+SINR_opt),2));
    end
end


      
