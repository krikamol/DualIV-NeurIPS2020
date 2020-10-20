function mse_value = sim_pred_dual_demand_hype(design, N, rho)

[f,sim,x_vis,f_vis]=get_design(design);

% simulate data, split into stage 1 and 2 samples
[x,y,z]=sim(f,N, rho);

p=x(:,1);
t=x(:,2);
s=x(:,3);

vp=median_inter(p);
vt=median_inter(t);
vs=median_inter(s);
vz=median_inter(z);
vy=median_inter(y);

ptest=x_vis(:,1);
ttest=x_vis(:,2);
stest=x_vis(:,3);

K_xx=get_K_matrix(p,p,vp).*get_K_matrix(t,t,vt).*get_K_matrix(s,s,vs);
K_Xtest=get_K_matrix(p,ptest,vp).*get_K_matrix(t,ttest,vt).*get_K_matrix(s,stest,vs);
K_zz = get_K_matrix(z, z, vz);
K_yy = get_K_matrix(y, y, vy);

K = K_xx;
%L = K_yy.*K_zz;

% alternative way of computing L
yz = [y, z];
vyz = 90000;
% vyz = 150000;

Vmat = [vy, vyz; vyz, vz];
L_yzyz = get_k_matrix_2d(yz, yz, Vmat);
L = L_yzyz;


% split the data into train and validation sets
N = size(K,1);
m = round(N.*0.5);
idx = randperm(N);
idxm = idx(1:m);
idxn = idx((m+1):N);
K1 = K(idxm,idxn);
K2 = K(idxn,idxn);
L1 = L(idxm,idxn);
L2 = L(idxn,idxn);
y1 = y(idxm);
y2 = y(idxn);

lambda_0=[log(0.001),log(0.001)];
% lambda_0=[log(0.05),log(0.05)];
DUALIV_obj=@(params) DUALIV_loss(K2,L1,L2,y2,exp(params));
[lambdas_star]=fminunc(DUALIV_obj,lambda_0);
beta = dualiv(K,L,y,exp(lambdas_star));

% ----------------

y_vis_dual = K_Xtest'*beta;

mse_value = mse(y_vis_dual, f_vis);


end