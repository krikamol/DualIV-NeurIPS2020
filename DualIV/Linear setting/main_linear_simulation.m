clear all;
close all;
clc
% addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV');
% addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/aux'); %test;
% addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/figures');
% addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/results'); %test;
% addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/utils');

%% Simulation 1

n=200;              % sample size
n_eval=500;        % evaluation size
n_trials=100;      % number of trials
beta=0.8;          % structural parameter
%g=0.5:0.05:0.9;   % tradeoff parameter
g=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];  % tradeoff parameter
d=1;           % number of instruments

f_vis = @(x,b) b.*x;

n_tradeoffs = length(g);
n_instruments = length(d);

results_mse_2sls=zeros(n_trials,n_tradeoffs,n_instruments);
results_mse_dualiv=zeros(n_trials,n_tradeoffs,n_instruments);

rng('default') % For reproducibility
for i=1:n_tradeoffs
    gg = g(i);
    for j=1:n_instruments
        dd = d(j);
        for k=1:n_trials
            [x,y,z] = simulate(n, beta, gg, dd);
            beta_2sls = estimate_2sls(x,y,z);
            beta_dualiv = estimate_dualiv(x,y,z,1e-8,100);
            
            [xt,yt,zt] = simulate(n_eval, beta, gg, dd);
            results_mse_2sls(k,i,j)=mse(f_vis(xt,beta_2sls)',yt);
            results_mse_dualiv(k,i,j)=mse(f_vis(xt,beta_dualiv)',yt);
        end
    end
end

average_mse_2sls = squeeze(mean(results_mse_2sls,1));
average_mse_dualiv = squeeze(mean(results_mse_dualiv,1));

semilogy(g,average_mse_2sls,'-','LineWidth',3); hold on;
semilogy(g,average_mse_dualiv,'--','LineWidth',3);
xlabel('$\gamma$','FontSize',20,'Interpreter','latex');
ylabel('MSE','FontSize',20);
legend({'2SLS','DualIV'},'Location','northwest','FontSize',20);
axis square; axis tight;
grid on;
hold off;

%% Simulation 2 (The effect of weak instrument)

n=200;              % sample size
n_eval=500;         % evaluation size
n_trials=100;       % number of trials
beta=0.8;           % structural parameter
g=0.05:0.05:0.9;   % tradeoff parameter
d=2;                % number of instruments

f_vis = @(x,b) zscore(b.*x);

n_tradeoffs = length(g);

results_mse_2sls    = zeros(n_trials,n_tradeoffs);
results_mse_dualiv  = zeros(n_trials,n_tradeoffs);

rng('default') % For reproducibility
for i=1:n_tradeoffs
    gg = g(i);
    for k=1:n_trials
        [x,y,z] = simulate(n, beta, gg, d);
        beta_2sls = estimate_2sls(x,y,z);
        beta_dualiv = estimate_dualiv(x,y,z,1e-6,300);

        [xt,yt,zt] = simulate(n_eval, beta, gg, d);
        results_mse_2sls(k,i)=mse(f_vis(xt,beta_2sls)',yt);
        results_mse_dualiv(k,i)=mse(f_vis(xt,beta_dualiv)',yt);
    end
end

avg_mse_2sls = squeeze(mean(results_mse_2sls,1));
avg_mse_dualiv = squeeze(mean(results_mse_dualiv,1));
%results = results_mse_2sls - results_mse_dualiv;
%boxplot(results,'PlotStyle','compact');
plot(g,avg_mse_2sls); hold on;
plot(g,avg_mse_dualiv);

%% Simulation 3 (The effect of weak instrument and number of irrelevant instruments)

n=200;               % sample size
n_eval=500;         % evaluation size
n_trials=100;        % number of trials
beta=0.8;          % structural parameter
g=0.05:0.05:0.95;    % tradeoff parameter
d=1:1:5;           % number of instruments

f_vis = @(x,b) zscore(b.*x);

n_tradeoffs = length(g);
n_instruments = length(d);

results_mse_2sls=zeros(n_trials,n_tradeoffs,n_instruments);
results_mse_dualiv=zeros(n_trials,n_tradeoffs,n_instruments);

rng('default') % For reproducibility
for i=1:n_tradeoffs
    gg = g(i);
    for j=1:n_instruments
        dd = d(j);
        for k=1:n_trials
            [x,y,z] = simulate(n, beta, gg, dd);
            beta_2sls = estimate_2sls(x,y,z);
            beta_dualiv = estimate_dualiv(x,y,z,1e-5,100);
            
            [xt,yt,zt] = simulate(n_eval, beta, gg, dd);
            results_mse_2sls(k,i,j)=mse(f_vis(xt,beta_2sls)',yt);
            results_mse_dualiv(k,i,j)=mse(f_vis(xt,beta_dualiv)',yt);
        end
    end
end

results = squeeze(mean(results_mse_2sls - results_mse_dualiv,1));
plot(g,results(:,1)); hold on;
plot(g,results(:,2));
plot(g,results(:,3));
plot(g,results(:,4));
plot(g,results(:,5));