function H_nonoise=LoS_h_asym(M,d)

    f=2.4*(10^9);
    c=3*(10^8);
    
    lambda=c/f;
    
    H_nonoise=exp(-1i*2*pi*d/lambda)*ones(1,M);
end




