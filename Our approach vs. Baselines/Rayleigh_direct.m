% This function generates the direct path channels 
function [H_noisy,H_nonoise]=Rayleigh_direct(K,est_SNR)

    H_nonoise = zeros(1,K,K);
    H_noisy = zeros(1,K,K);
    
    N_P=-est_SNR;
    N_P=10.^(N_P/10);
    
    for k=1:K
        H_R=normrnd(0,1,[1,K])./sqrt(2);
        H_I=normrnd(0,1,[1,K])./sqrt(2);
        H=H_R+1i*H_I;

        noise_vec_R=normrnd(0,sqrt(N_P),[1,K])./sqrt(2);
        noise_vec_I=normrnd(0,sqrt(N_P),[1,K])./sqrt(2);
        noise_vec=noise_vec_R+1i*noise_vec_I;

        H_nonoise(1,k,:) = H;
        H_noisy(1,k,:)=H+noise_vec;
    end

end




