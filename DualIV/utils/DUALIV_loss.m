function [ loss ] = DUALIV_loss(K2,L1,L2,y2,params)

n = size(K2,1);
beta = dualiv(K2,L2,y2,params);
alpha = (make_psd(L2) + 1e-8.*eye(n))\(K2*beta - y2);

loss = mean((L1*alpha).^2);

end