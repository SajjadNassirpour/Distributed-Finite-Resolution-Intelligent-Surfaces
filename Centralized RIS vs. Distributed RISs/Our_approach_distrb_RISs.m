% -- Our optimization appraoch based on sigmoid filled function with K distributed RISs -- %

function [sum_rate_with_Tx, sum_rate_without_Tx, BF_our_ratio]=Our_approach_distrb_RISs(K,M,N,Pt_dB,H_nonoise,G1_nonoise,G1_far_nonoise,G2_nonoise,sigma2,Tx_loc,Rx_loc,RIS_loc,alpha_d,alpha_d1,alpha_d2,c0_dB,mu,r,divide_rate,ini_chance,opt_objective,LN_check)
      
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

    % Noise power
    sigma2=10^(sigma2/10);
    sigma2=sigma2*10^(-3);
    
    c0=10^(c0_dB/10);
  
       
    % Possible values for each element at RIS
    comp_val=0:(2*pi)/N:(2*pi-((2*pi)/N));
    All_values=exp(1j*comp_val);
    
    
    % Optimization Parameters
    input_mu=mu;
    input_r=r;
    input_LN_check=LN_check;
    
    % Number of evaluations
    Eval=1;
   
        
    % The channels
    H=sqrt(c0)*squeeze(H_nonoise(1,:,:));
    G1=sqrt(c0)*squeeze(G1_nonoise(:,1,:));
    G2=sqrt(c0)*squeeze(G2_nonoise(:,1,:,:));


    for ini_p=1:ini_chance
    % Generate random initial switch configuration
        for User=1:K
            ar = randi([1,N],1,M);
            theta_ini(User,:) = All_values(ar);
        end

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

                [S_star,eval]=LS_algo_diff_pow(K,N,theta_ini,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,mode,r,opt_objective);

            else
                eval=0;
                S_star=theta_ini;
            end

            Eval=Eval+eval;
            S_star_mat(kk+1,:,:)=S_star;

            if kk==0
                min_f=obj_func_SINR(K,squeeze(S_star_mat(kk+1,:,:)),H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
                S_final=squeeze(S_star_mat(kk+1,:,:));
            end


            if kk>0
               if obj_func_SINR(K,squeeze(S_star_mat(kk+1,:,:)),H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective)>= min_f
                   jump=1;

               else
                   min_f=obj_func_SINR(K,squeeze(S_star_mat(kk+1,:,:)),H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,opt_objective);
                   S_final=squeeze(S_star_mat(kk+1,:,:));
                   r=r0;
               end
            end

            if jump==0
                ii=1;
                item=1;
                theta_user_idx=0;
                S_star=squeeze(S_star_mat(kk+1,:,:));
                mode=2;
                [S_bar,eval]=LS_algo_diff_pow(K,N,S_star,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,mode,r,opt_objective);
                Eval=Eval+eval;
                theta_ini = S_bar;
                kk=kk+1;
            end

            if ii>=M+1
                if r>=mu
                    r=r/divide_rate;
                    ii=1;
                    item=1;
                    theta_user_idx=0;
                    kk=kk-1;
                    S_star=squeeze(S_star_mat(kk+1,:,:));
                    mode=2;
                    [S_bar,eval]=LS_algo_diff_pow(K,N,S_star,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,mode,r,opt_objective);
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
                S_star=squeeze(S_star_mat(kk+1,:,:));
                kk=kk+1;
                extra=0;
                while extra==0 && ii<M+1

                    theta_user_idx=theta_user_idx+1;

                    if S_star(theta_user_idx,ii)~=All_values(item)
                        S_star_j0=S_star;
                        S_star_j0(theta_user_idx,ii)=All_values(item);
                        if sum(S_star_j0(theta_user_idx,:)==0)~=M
                            extra=1;
                            mode=2;
                            [S_bar,eval]=LS_algo_diff_pow(K,N,S_star_j0,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0,mode,r,opt_objective);
                            Eval=Eval+eval;
                            theta_ini = S_bar;

                        end
                    end

                    if item==N
                        ii=ii+1;
                        item=1;
                    end

                    if theta_user_idx==K
                        item=item+1;
                        theta_user_idx=0;
                    end
                end
            end



        end
    end

    % Final evaluation based on optimal switch configuration
    H=sqrt(c0)*squeeze(H_nonoise(1,:,:)); 
    G1=sqrt(c0)*squeeze(G1_nonoise(:,1,:));
    G1_far=sqrt(c0)*squeeze(G1_far_nonoise(:,1,:,:));
    G2=sqrt(c0)*squeeze(G2_nonoise(:,1,:,:));
    [SINR_with_Tx, SINR_without_Tx] = obj_func_SINR_final(K,S_final,H,dist_mat_Tx_Rx,dist_mat_Tx_RIS,dist_mat_RIS_Rx,G1,G1_far,G2,Pt,sigma2,alpha_d,alpha_d1,alpha_d2,c0);

    
    % The tatio of (# of evaluations of Brute Force/ our approach)
    BF_our_ratio=(N^M)/mean(Eval);
    
    % --- Sum_rate with considering the channels between Txi and RISj --- %
    sum_rate_with_Tx=mean(sum(log2(1+SINR_with_Tx),2));


     % --- Sum_rate without considering the channels between Txi and RISj --- %
    sum_rate_without_Tx=mean(sum(log2(1+SINR_without_Tx),2));
end


      
