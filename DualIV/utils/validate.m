function [ optparams, optloss ] = validate(K,L,y,m,L1_params,L2_params)
%VALIDATE performs grid search to find the best parameter setting

% split the data into train and validation sets
N = size(K,1);
idx = randperm(N);
idxm = idx(1:m);
idxn = idx((m+1):N);
K1 = K(idxm,idxn);
K2 = K(idxn,idxn);
L1 = L(idxm,idxm);
L2 = L(idxn,idxn);
y1 = y(idxm);
y2 = y(idxn);

%%
optparams = [L1_params(1), L2_params(1)];
optloss = Inf;
for i=1:length(L1_params)
    for j=1:length(L2_params)
        curparams = [L1_params(i),L2_params(j)];
        curloss = loss(K1,K2,L1,L2,y1,y2,curparams);
        if curloss < optloss
            optloss = curloss;
            optparams = curparams;
        end
    end
end    
    
end