% This function generates the channels between Txi and RIS

function [H_noisy,H_nonoise]=Rician_h_center_RIS(M,d,est_SNR,beta)


    % Rician coefficients
    rician_var_coeff=sqrt(1/(1+beta));
    rician_mean_coeff=sqrt(beta*(1+beta));
    
    f=2.4*(10^9);
    c=3*(10^8);
    
    lambda=c/f;
        
    % H_nonoise = zeros(iterations,1,num_RIS_elements);
    H_nonoise = zeros(1,1,M);
    H_noisy = zeros(1,1,M);
    
    N_P=-est_SNR;
    N_P=10.^(N_P/10);
 
    H_R=normrnd(0,1,[1,M])./sqrt(2);
    H_I=normrnd(0,1,[1,M])./sqrt(2);
    H=rician_var_coeff*(H_R+1i*H_I)+(rician_mean_coeff*exp(-1i*2*pi*d/lambda)/sqrt(2));
    H_nonoise(1,1,:) = H;

    noise_vec_R=normrnd(0,sqrt(N_P),[1,M])./sqrt(2);
    noise_vec_I=normrnd(0,sqrt(N_P),[1,M])./sqrt(2);
    noise_vec=noise_vec_R+1i*noise_vec_I;
    noise_vec=rician_var_coeff*(noise_vec)+sqrt(N_P)*(rician_mean_coeff*exp(-1i*2*pi*d/lambda)/sqrt(2));
    H_noisy(1,1,:)=H+noise_vec;

end




