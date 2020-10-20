function [ optparams ] = validate_minimax(K,L,y,m,L1_params)
%VALIDATE performs grid search to find the best parameter setting

% split the data into train and validation sets
N = size(K,1);
n = N-m;
idx = randperm(N);
idxm = idx(1:m);
idxn = idx((m+1):N);
K1 = K(idxm,idxn);
K2 = K(idxn,idxn);
L1 = L(idxm,idxn);
L2 = L(idxn,idxn);
y1 = y(idxm);
y2 = y(idxn);

%%
x0 = log(0.05);
[optlambda2,fval] = fminimax(@loss_fun,x0);

[~,optidx1] = max(fval);
optparams = [L1_params(optidx1),exp(optlambda2)];

%% loss function
    function f = loss_fun(x)
        for k=1:length(L1_params)
            beta = dualiv(K2,L2,y2,[L1_params(k),exp(x)]);
            alpha = (make_psd(L2) + n.*L1_params(k).*eye(n))\(K2*beta - y2);
            f(k) = (beta.'*(K1.')*L1*alpha)./n - (y1.'*L1*alpha)./n - (alpha.'*(L1.')*L1*alpha)./(2*n);
        end
    end

end
