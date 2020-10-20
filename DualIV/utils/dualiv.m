function [ beta ] = dualiv(K,L,y,params)
%DUALIV estimates the coefficient vector beta

n = size(K,1); 
M = K*((make_psd(L) + n.*params(1)*eye(n))\L);
beta = (make_psd(M*K) + n.*params(2).*make_psd(K))\(M*y);

end

