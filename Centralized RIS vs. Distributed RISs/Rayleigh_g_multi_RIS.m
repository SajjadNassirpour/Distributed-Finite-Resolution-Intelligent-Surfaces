function H_nonoise=Rayleigh_g_multi_RIS(K,M)

    H_nonoise = zeros(1,M,K);
    for k1=1:K
        H_R=normrnd(0,1,[M,1])./sqrt(2);
        H_I=normrnd(0,1,[M,1])./sqrt(2);
        H=H_R+1i*H_I;
        H_nonoise(1,:,k1) = H;
    end

end




