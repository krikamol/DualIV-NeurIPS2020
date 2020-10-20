clear all;
close all;
clc
addpath('../utils');
addpath('../../aux2');
warning('off','all');

%% Multi trial simulations
% rng('default')
n_trials=20;
all_sample_sizes = [20, 30, 50];
rho = 0.1;
design='HLLT'; % Only for demand experiment
[f,sim,x_vis,f_vis]=get_design(design);


for i=1:1:length(all_sample_sizes)
    N = all_sample_sizes(i);
 
    % KIV
    results_mse_kiv=zeros(n_trials,1);
    options = optimoptions('fminunc', 'Display', 'off');
    for i=1:n_trials
        fprintf('\n N=%s, KIV Trial: %% %s\n', num2str(N), num2str(i/n_trials*100));
        pause(0.05)
        results_mse_kiv(i)= sim_pred_kiv_demand(design, N, rho, options);
    end

    % DualIV
    results_mse_dualiv=zeros(n_trials,1);
    options = optimoptions('quadprog', 'Display', 'off');
    for i=1:n_trials
        fprintf('\n N=%s DualIV Trial: %% %s\n', num2str(N), num2str(i/n_trials*100));
        pause(0.05)
        results_mse_dualiv(i)= sim_pred_dual_demand(design, N, rho, options);
    end
    % 2SLS
    results_mse_2sls=zeros(n_trials,1);
    for i=1:n_trials
        fprintf('\n N=%s 2SLS Trial: %% %s\n', num2str(N), num2str(i/n_trials*100));
        pause(0.05)
        results_mse_2sls(i)= sim_pred_2sls_demand(design, N, rho, options);
    end
% 
    mse_methods = {log(results_mse_kiv)/log(10), log(results_mse_dualiv)/log(10), log(results_mse_2sls)/log(10)};
    sample_sizes = [N];
    mses_save_path = fullfile('MSE_results',sprintf('MSE_N=%d_trials=%d_rho=%0.2f.mat', N, n_trials,rho));
    save(mses_save_path, 'mse_methods', 'sample_sizes')
    fprintf('Sample size N=%d is done', N);
end

