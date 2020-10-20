clear all;
close all;
clc
% addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV')
% addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/aux'); %test;
% addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/figures');
% addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/results'); %test;
% addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/utils')

%% Generative process

N=100;      % sample size
beta=0.8;   % structural parameter
g=0.8;      % tradeoff parameter
d=1;        % number of instruments

f_vis = @(x,b) b.*x;

rng('default') % For reproducibility
[x,y,z] = simulate(N, beta, g, d);
%% 

% visualize design
subplot(1,2,1);
plot(x,f_vis(x,beta),'LineWidth',5);
hold on;
scatter(x,y,36,[.5 .5 .5]);
xlabel('x','FontSize',20)
ylabel('y','FontSize',20)
legend({'Structural function','Data'},'Location','southeast','FontSize',20);
hold off;
%saveas(gcf,fullfile('./figures',strcat('/linear_model')),'png');

subplot(1,2,2);
scatter(z(1,:),x,36,[.5 .5 .5]);
xlabel('z','FontSize',20)
ylabel('x','FontSize',20)
legend({'Data'},'Location','southeast','FontSize',20);
hold off;

%% 2SLS

beta_2sls = estimate_2sls(x,y,z);

%% DualIV with generalized LSQR

beta_dualiv = estimate_dualiv(x,y,z,1e-6,100);

%% visualize estimator

figure()

plot(x, f_vis(x,beta),'LineWidth',5);
hold on;
scatter(x,y,36,[.5 .5 .5]);

plot(x, f_vis(x,beta_2sls),'--g','LineWidth',5);
plot(x, f_vis(x,beta_dualiv),'--r','LineWidth',5);

set(gca, 'FontSize', 20)
xlabel('x','FontSize',25)
ylabel('y','FontSize',25)
legend({'Structural function','Data', '2SLS', 'DualIV'},'Location','southeast','FontSize',20);
hold off;
% saveas(gcf,fullfile('figures',strcat('DualIV/estimated_design_',design)),'png');
% remove whitespace
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print(gcf, fullfile('../../figures',strcat('DualIV_','design','_',num2str(1))), '-dpdf')
% mse
disp('mse for DualIV:');
disp(mse(y_vis_dual, f_vis));
