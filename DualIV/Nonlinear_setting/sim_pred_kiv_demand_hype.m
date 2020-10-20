function mse_val = sim_pred_kiv_demand_hype(design, N, rho)


[f,sim,x_vis,f_vis]=get_design(design);

% simulate data, split into stage 1 and 2 samples
[x,y,z]=sim(f,N, rho);


df=get_K_demand(x,y,z,x_vis);

% initialize hyperparameters for tuning
lambda_0=log(0.05); %hyp1=lambda. log(0.05) for NP
xi_0=log(0.05); %hyp2=xi. log(0.05) for NP
%  
% stage 1 tuning
KIV1_obj=@(lambda) KIV1_loss(df,exp(lambda)); %exp to ensure pos; fixed vx,vz
lambda_star=fminunc(KIV1_obj,lambda_0);
 
% stage 2 tuning
KIV2_obj=@(xi) KIV2_loss(df,[exp(lambda_star),exp(xi)]); %exp to ensure pos
xi_star=fminunc(KIV2_obj,xi_0);

% evaluate on full sample using tuned hyperparameters
y_vis=KIV_pred(df,[exp(lambda_star),exp(xi_star)],3); %exp to ensure pos
% y_vis=KIV_pred(df,[exp(lambda_0),exp(xi_0)],3); %exp to ensure pos

% mse
mse_val = mse(y_vis,f_vis);

end

