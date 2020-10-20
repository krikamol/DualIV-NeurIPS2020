clear;
addpath('../aux2'); %test;
addpath('../figures');
addpath('../results'); %test;

design='HLLT'; %NP, HLLT
[f,sim,x_vis,f_vis]=get_design(design);

% sample size
N=1000;

% simulate data, split into stage 1 and 2 samples
rng('default');
rho=0.5;
[x,y,z]=sim(f,N,rho);

% visualize design
t_vis=linspace(0,10,1000);
psi_vis=get_psi(t_vis);
plot(t_vis,psi_vis,'LineWidth',2,'Color',[0 0 0]);
hold on;
%scatter(x,y,36,[.5 .5 .5]);
hold off;
xlabel('t','FontSize',20)
ylabel('\psi(t)','FontSize',20)
hold off;
saveas(gcf,fullfile('../figures',strcat('design_',design)),'epsc');

%% KIV - IV, causal validation with lengths simply set

df=get_K(x,y,z,x_vis);

% initialize hyperparameters for tuning
lambda_0=log(0.05); %hyp1=lambda. log(0.05) for NP
xi_0=log(0.05); %hyp2=xi. log(0.05) for NP

% stage 1 tuning
KIV1_obj=@(lambda) KIV1_loss(df,exp(lambda)); %exp to ensure pos; fixed vx,vz
lambda_star=fminunc(KIV1_obj,lambda_0);

% stage 2 tuning
KIV2_obj=@(xi) KIV2_loss(df,[exp(lambda_star),exp(xi)]); %exp to ensure pos
xi_star=fminunc(KIV2_obj,xi_0);

% evaluate on full sample using tuned hyperparameters
y_vis=KIV_pred(df,[exp(lambda_star),exp(xi_star)],3); %exp to ensure pos

% mse
disp('mse:');
disp(mse(y_vis,f_vis));