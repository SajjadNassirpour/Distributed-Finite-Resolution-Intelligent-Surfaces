function [H_noisy,H_nonoise]=Rayleigh_G1_multi_RIS(M,d,iters,est_SNR)




        
    
    H_mat = zeros(iters,1,M);
    for iter=1:iters
        
        H_R=normrnd(0,1,[1,M])./sqrt(2);
        H_I=normrnd(0,1,[1,M])./sqrt(2);
        H=H_R+1i*H_I;
        H_mat(iter,:) = H;
        
    end
    H_nonoise=H_mat;


    N_P=-est_SNR;
    N_P=10.^(N_P/10);
    

    H_noisy=[];
%     H_mat=[];
%     for iter=1:iters
% 
%         for kk=1:K
%         noise_vec1_R(kk,:)=normrnd(0,sqrt(N_P),[1,M])./sqrt(2);
%         noise_vec1_I(kk,:)=normrnd(0,sqrt(N_P),[1,M])./sqrt(2);
%         noise_vec1(kk,:)=noise_vec1_R(kk,:)+1i*noise_vec1_R(kk,:);
%         end
% 
% 
%         H_noisy(iter,:,:)=squeeze(H_nonoise(iter,:,:))+noise_vec1;
%     end
end




