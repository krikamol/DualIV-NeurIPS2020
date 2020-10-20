function [ cvloss ] = loss(K1,K2,L1,L2,y1,y2,params)

m = size(K1,1);
beta = dualiv(K2,L2,y2,params);
%N = ((L1 + m.*params(1)*eye(m))\L1).';
%cvloss = (beta.'*K1.'*N*K1*beta)./(2*m) - (y1.'*N*K1*beta)./m + (y1.'*N*y1)./(2*m); % + (params(2)/2).*beta.'*K2*beta; 
cvloss = (beta.'*K1.'*K1*beta)./(2*m) - (y1.'*K1*beta)./m + (y1.'*y1)./(2*m); % + (params(2)/2).*beta.'*K2*beta; 

end