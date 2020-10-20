function mse_value = sim_pred_2sls_demand(design, N, rho, options)
[f,sim,x_vis,f_vis]=get_design(design);
% simulate data, split into stage 1 and 2 samples
[x,y,z]=sim(f,N, rho);

Pz = z.'*((z*z.')\z);
A = ((x*Pz*x.')\x)*Pz;
beta_2sls = A'*y;
y_vis_2sls = x_vis * beta_2sls;
mse_value = mse(y_vis_2sls,f_vis);
end