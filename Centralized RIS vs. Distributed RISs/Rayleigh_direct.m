% This function generated the direct path noiseless channels

function H_nonoise=Rayleigh_direct(K)
    H_nonoise = zeros(1,K,K);
    for k=1:K
        H_R=normrnd(0,1,[1,K])./sqrt(2);
        H_I=normrnd(0,1,[1,K])./sqrt(2);
        H=H_R+1i*H_I;
        H_nonoise(1,k,:) = H;
    end
end




