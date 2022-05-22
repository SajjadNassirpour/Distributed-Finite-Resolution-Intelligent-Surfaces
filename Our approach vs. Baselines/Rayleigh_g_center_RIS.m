% This function generates the channels between RIS and Rxi

function [H_noisy,H_nonoise]=Rayleigh_g_center_RIS(K,M,est_SNR)


    % H_nonoise = zeros(iterations,num_RIS_elements,num_Rx);
    H_nonoise = zeros(1,M,K);
    H_noisy= zeros(1,M,K);
    
    N_P=-est_SNR;
    N_P=10.^(N_P/10);
    

    for k1=1:K
        H_R=normrnd(0,1,[M,1])./sqrt(2);
        H_I=normrnd(0,1,[M,1])./sqrt(2);
        H=H_R+1i*H_I;

        noise_vec_R=normrnd(0,sqrt(N_P),[M,1])./sqrt(2);
        noise_vec_I=normrnd(0,sqrt(N_P),[M,1])./sqrt(2);
        noise_vec=noise_vec_R+1i*noise_vec_I;

        H_nonoise(1,:,k1) = H;
        H_noisy(1,:,k1)=H+noise_vec;
    end


end




