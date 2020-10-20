clear all;
close all;
clc
%addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV')
%addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/aux2'); %test;
%addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/figures');
%addpath('/Volumes/krikamol/Works/Projects/DualIV/KernelIV_MATLAB/DualIV/results'); %test;
addpath('../utils');

%% Generative process

N=300;          % sample size
beta=0.7;       % structural parameter
est_beta=0.4;   % estimated parameter
g=0.2;          % tradeoff parameter
d=1;            % number of instruments

f_vis  = @(x,b) b.*x;
dual_f = @(cb,y,z) cb.*((1-g).*z) - y; 

rng('default') % For reproducibility
[x,y,z] = simulate(N, beta, g, d);

%% 
beta_ols = polyfit(x.',y.',1);

% visualize design
figure();
scatter(x,y,48,[.5 .5 .5],'filled');
hold on;
plot(x,f_vis(x,beta),'LineWidth',6);
hold on;
plot(x,f_vis(x,beta_ols(1))+beta_ols(2),'LineWidth',6,'LineStyle',':');
hold on;
plot(x,f_vis(x,est_beta),'--','LineWidth',6);
xlabel('$X$','FontSize',24,'Interpreter','latex')
ylabel('$Y$','FontSize',24,'Interpreter','latex')
legend({'Data','True','OLS','Estimated'},'Location','southeast','FontSize',24);
title('The second stage regression','FontSize',36,'Interpreter','latex');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',24);
axis square; 
axis tight;
grid on;
box on;
hold off;

yylim = ylim;

saveas(gcf, 'visualize_dual_1.eps', 'eps2c');

%subplot(1,3,2);
figure();
scatter(z(1,:),x,48,[.5 .5 .5],'filled');
hold on;
plot(z(1,:),f_vis(z,1-g),'LineWidth',6);
xlabel('$Z$','FontSize',24,'Interpreter','latex')
ylabel('$X$','FontSize',24,'Interpreter','latex')
legend({'Data','True'},'Location','southeast','FontSize',24);
title('The first stage regression','FontSize',36,'Interpreter','latex');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',24);
axis square; 
axis tight;
grid on;
box on;
hold off;

zzlim = xlim;

saveas(gcf, 'visualize_dual_2.eps', 'eps2c');

% heatmap of the dual function
[Z,Y] = meshgrid(zzlim(1):0.25:zzlim(2),yylim(1):0.25:yylim(2));
dual_evals = dual_f(est_beta,Y,Z);
%min_dual = min(dual_evals);
%max_dual = max(dual_evals);
%dual_evals = 2.*(dual_evals - min_dual)./(max_dual-min_dual) - 1;

%subplot(1,3,3);
figure();
contourf(Z,Y,dual_evals,20,'--','LineWidth',1);
colorbar;
shading interp;
xlabel('$Z$','FontSize',24,'Interpreter','latex');
ylabel('$Y$','FontSize',24,'Interpreter','latex');
title('The optimal dual function','FontSize',36,'Interpreter','latex')
xlim(zzlim); ylim(yylim);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',24);
axis square; 
axis tight;
grid on;
box on;
hold off;

saveas(gcf, 'visualize_dual_3.eps', 'eps2c');