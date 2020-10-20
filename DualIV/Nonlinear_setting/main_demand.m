clear all;
close all;
clc
%addpath('../aux'); %test;
%addpath('../figures');
%addpath('../results'); %test;
addpath('../utils');
addpath('../../aux2');
warning('off','all');
%% Generative process

design='HLLT'; %NP, SSG, HLLT
[f,sim,x_vis,f_vis]=get_design(design);

% sample size
N=50;
rho = 0.4;

% simulate data, split into stage 1 and 2 samples
% rng('default');
[x,y,z]=sim(f, N, rho);

% visualize design
t_vis=linspace(0,10,1000);
psi_vis=get_psi(t_vis);
plot(t_vis,psi_vis,'LineWidth',2,'Color',[0 0 0]);
hold on;
%scatter(x,y,36,[.5 .5 .5]);
hold off;
xlabel('t','FontSize',20)
ylabel('\psi(t)','FontSize',20)
hold off; grid on; axis tight;
legend({'Structural function','Data'},'Location','southeast','FontSize',20);
saveas(gcf,fullfile('../../results',strcat('/psi_function_',design)),'eps');


%% 2SLS
Pz = z.'*((z*z.')\z);
A = ((x*Pz*x.')\x)*Pz;
beta_2sls = A'*y;

y_vis_2sls = x_vis*beta_2sls;
disp('mse 2sls:');
disp(mse(y_vis_2sls,f_vis));

%% Sieve - IV
% df=get_K_demand(x,y,z,x_vis);
% lambda_0=log(0.05); %hyp1=lambda. log(0.05) for NP
% xi_0=log(0.05); %hyp2=xi. log(0.05) for NP
% % stage 1 tuning
% sieve1_obj = @(lambda) sieve1_loss(df,exp(lambda));
% lambda_star=fminunc(sieve1_obj,lambda_0);
% 
% sieve2_obj=@(xi) sieve2_loss(df,[exp(lambda_star),exp(xi)]); %exp to ensure pos
% xi_star=fminunc(sieve2_obj,xi_0);
% 
% y_vis=sieveIV_pred(df,[exp(lambda_star),exp(xi_star)],3); %exp to ensure pos
% % mse
% disp('mse kiv:');
% disp(mse(y_vis,f_vis));



%% KIV - IV, causal validation with lengths simply set
 
df=get_K_demand(x,y,z,x_vis);

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
disp('mse kiv:');
disp(mse(y_vis,f_vis));


%% Dual - IV
close all
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
Vmat = [vy, vyz; vyz, vz];
L_yzyz = get_k_matrix_2d(yz, yz, Vmat);
L = L_yzyz;
lambda1 = 0.001; %improves condition number of A for a better inv(A)
gamma = N * norm(L * L, 2) / norm(K * L, 2)^2;
A = L * L + 1 / N * gamma * L * (K * K) * L + lambda1 * eye(N);
Ainv = inv(A);
lambda2 = 0.001; % reqularized the solution to the following quadratic programming
Q = 2 * K' * L' * Ainv * L * K + lambda2 * eye(N);
R = - 2 * K' * L' * Ainv * L * y;
[beta, feval] = quadprog(Q, R);

y_vis_dual = K_Xtest'*beta;

% mse
disp('mse dual method 1:');
disp(mse(y_vis_dual,f_vis));
% 

%%
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

lambda_0=[log(0.05),log(0.05)];
DUALIV_obj=@(params) DUALIV_loss(K2,L1,L2,y2,exp(params));
[lambdas_star]=fminunc(DUALIV_obj,lambda_0);
beta = dualiv(K,L,y,exp(lambdas_star));

% ----------------

y_vis_dual = K_Xtest'*beta;


% mse
disp('mse dual method 2:');
disp(mse(y_vis_dual,f_vis));

% 

%% Multi trial simulations
% rng('default')
n_trials=20;
all_sample_sizes = [1000];
% all_sample_sizes = [1000, 5000];
rho = 0.9;
design='HLLT'; %NP, SSG, HLLT
for i=1:1:1
    fprintf('Sample size N=%d is done', N);
    N = all_sample_sizes(i);
    [f,sim,x_vis,f_vis]=get_design(design);
    % KIV - Multi trial simulations

    results_mse_kiv=zeros(n_trials,1);
    for i=1:n_trials
        fprintf('Trial %s for sample size %d\n', num2str(i), N);
        pause(0.05)
        results_mse_kiv(i)= sim_pred_kiv_demand_hype(design, N, rho);
    end

    % DualIV - Multi trial simulations

    % sample size
    results_mse_dualiv=zeros(n_trials,1);
    for i=1:n_trials
        fprintf('Trial %s for sample size %d\n', num2str(i), N);
        if mod(i,5)==0
            disp(num2str(i));
        end
%         results_mse_dualiv(i)= sim_pred_dual_demand(design, N, rho);
        results_mse_dualiv(i)= sim_pred_dual_demand(design, N, rho, true);
    end
    
    % 2SLS - Multi trial simulations
    results_mse_2sls=zeros(n_trials,1);
    for i=1:n_trials
        fprintf('Trial %s for sample size %d\n', num2str(i), N);
        if mod(i,5)==0
            disp(num2str(i));
        end
        results_mse_2sls(i)= sim_pred_2sls_demand(design, N, rho);
    end

    mse_methods = {log(results_mse_kiv)/log(10), log(results_mse_dualiv)/log(10), log(results_mse_2sls)/log(10)};
    sample_sizes = [N];
    
%     for i=1:1:length(mse_methods)
%         figure()
%         bp = boxplot(mse_methods{i}, sample_sizes,...
%             'BoxStyle', 'outline', 'Colors', colors{i},...
%             'Notch', 'off', 'MedianStyle', 'line',...
%             'Widths', 0.5, 'Symbol','o',...
%             'OutlierSize',5)
%         set(bp, 'LineWidth', 2)
%         hold on
%         errorbar(1, mean(mse_methods{i}), - std(mse_methods{i}) / sqrt(n_trials), std(mse_methods{i}) / sqrt(n_trials), 's')
%     end
    
%     for i=1:1:length(mse_methods)
% 
%     end
    
    
        

    mses_save_path = fullfile('../../results',sprintf('MSE_N=%d_trials=%d_rho=%0.2f.mat', N, n_trials,rho));
    save(mses_save_path, 'mse_methods', 'sample_sizes')

    %% plot barplots
    save_path = fullfile('../../results',sprintf('NoHype_MSE_Barplots_samplesize=%d_trials=%d_rho=%0.2f.pdf', N, n_trials,rho));
    colors = {[1 0 0], [0 0 1], [0 1 0]};
    make_box_plots(mse_methods(1:2), sample_sizes, rho, colors, save_path)
%     make_error_bars(mse_methods(1:2), sample_sizes, rho, colors, save_path)

end


fprintf("KIV Mean:%f\n", mean(mse_methods{1}))
fprintf("DualIV Mean:%f\n", mean(mse_methods{2}))
fprintf("2SLS Mean:%f\n", mean(mse_methods{3}))

fprintf("KIV STD:%f\n", std(mse_methods{1}))
fprintf("DualIV STD:%f\n", std(mse_methods{2}))
fprintf("2SLS STD:%f\n", std(mse_methods{3}))



