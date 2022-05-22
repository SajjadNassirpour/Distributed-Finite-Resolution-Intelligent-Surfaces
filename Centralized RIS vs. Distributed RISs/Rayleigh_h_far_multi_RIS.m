function H_nonoise=Rayleigh_h_far_multi_RIS(K,M)


    H_nonoise = zeros(1,K,M);
    for k=1:K
        H_R=normrnd(0,1,[1,M])./sqrt(2);
        H_I=normrnd(0,1,[1,M])./sqrt(2);
        H=H_R+1i*H_I;
        H_nonoise(1,k,:) = H;
    end
end




