clear all;
close all;
clc
%addpath('../aux'); %test;
%addpath('../figures');
%addpath('../results'); %test;
addpath('../utils');
addpath('../../aux2');

%% Generative process

design='NP'; %NP, SSG, HLLT
[f,sim,x_vis,f_vis]=get_design(design);

% sample size
N=500;

% simulate data, split into stage 1 and 2 samples
rng('default');
[x,y,z]=sim(f,N);


% visualize design
plot(x_vis,f_vis,'LineWidth',5);
hold on;
scatter(x,y,36,[.5 .5 .5]);
xlabel('x','FontSize',20)
ylabel('y','FontSize',20)
legend({'Structural function','Data'},'Location','southeast','FontSize',20);
hold off;
saveas(gcf,fullfile('../../figures',strcat('DualIV/generative_design_',design)),'png');




%% KIV
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
y_vis_kiv=KIV_pred(df,[exp(lambda_star),exp(xi_star)],3); %exp to ensure pos

% visualize estimator
plot(x_vis, f_vis,'LineWidth',5);
hold on;
scatter(x,y,36,[.5 .5 .5]);
plot(x_vis,y_vis_kiv,'--r','LineWidth',5);
xlabel('x','FontSize',20)
ylabel('y','FontSize',20)
legend({'Structural function','Data','KernelIV'},'Location','southeast','FontSize',20);
hold off;
saveas(gcf,fullfile('../figures',strcat('KIV_',design)),'epsc');

% mse
disp('mse for KIV:');
disp(mse(y_vis_kiv, f_vis));


%% Dual - IV
close all
vx = median_inter(x);
K_xx = get_K_matrix(x, x, vx);

vy = median_inter(y);
K_yy = get_K_matrix(y, y, vy);

vz = median_inter(z);
K_zz = get_K_matrix(z, z, vz);

K = K_xx;
L = K_yy.*K_zz;

%yz = [y, z];
%vyz = (median_inter_2d(yz))^2;
% L_yy = get_K_matrix(y, y, vy);
% L_zz = get_K_matrix(z, z, vz);
% L_yzyz = times(L_yy, L_zz);

%vyz = 2800000;
%Vmat = [vy, vyz; vyz, vz];
%L_yzyz = get_K_matrix_2d(yz, yz, Vmat);

%L = L_yzyz;
%lambda1 = 0.001; %improves condition number of A for a better inv(A)

%gamma = N * norm(L * L, 2) / norm(K * L, 2)^2;
%A = L * L + 1 / N * gamma * L * (K * K) * L + lambda1 * eye(N);
% A  = eye(N);

%Ainv = inv(A);
%lambda2 = 0.01; % reqularized the solution to the following quadratic programming
%Q = 2 * K' * L' * Ainv * L * K + lambda2 * eye(N);
%R = - 2 * K' * L' * Ainv * L * y;
%[beta, feval] = quadprog(Q, R);

% --- new version
L1_params = 10.^(-12:0.5:-1);
L2_params = 10.^(-5:1:5);

m = round(N.*0.3);
%[optparams,cvloss] = validate(K,L,y,m,L1_params,L2_params);
optparams = validate_minimax(K,L,y,m,L1_params);

%optparams = [1e-12,1e-1];
beta = dualiv(K,L,y,optparams);

% ----------------

vx_vis = median_inter(x_vis);
K_xx_vis = get_K_matrix(x, x_vis, vx_vis);
y_vis_dual = K_xx_vis' * beta;
% 
% [x_test, y_test, z_test]=sim(f,N);
% figure()
% plot(x_vis, y_pred,'LineWidth',5);
% hold on;
% scatter(x_test, y_test, 36, [.5 .5 .5]);


% visualize estimator
figure()

plot(x_vis, f_vis,'LineWidth',5);
hold on;
scatter(x,y,36,[.5 .5 .5]);
plot(x_vis, y_vis_kiv,'--g','LineWidth',5);
plot(x_vis, y_vis_dual,'--r','LineWidth',5);
set(gca, 'FontSize', 20)
xlabel('x','FontSize',25)
ylabel('y','FontSize',25)
legend({'Structural function','Data', 'KIV', 'DualIV'},'Location','southeast','FontSize',20);
hold off;
% saveas(gcf,fullfile('../figures',strcat('DualIV/estimated_design_',design)),'png');
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
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

%print(gcf, fullfile('../../figures',strcat('DualIV_',design,'_',num2str(1))), '-dpdf')
% mse
disp('mse for DualIV:');
disp(mse(y_vis_dual, f_vis));




%% KIV - 100 simulations

rng('default')
design='NP'; %NP, SSG, HLLT
[f,sim,x_vis,f_vis]=get_design(design);

% sample size
N=1000;
n_trials=10;
results=zeros(n_trials,N);
results_mse_kiv=zeros(n_trials,1);

figure()
plot(x_vis,f_vis,'LineWidth',2);
hold on;
for i=1:n_trials
    if mod(i,10)==0
        disp(num2str(i));
    end
    results(i,:)= sim_pred(design,N)';
    results_mse_kiv(i)=mse(results(i,:)',f_vis);
    plot(x_vis,results(i,:),'--r','LineWidth',2);
end

% visualize 100 samples
means_kiv=mean(results);
sorted=sort(results);
fifth_kiv=sorted(5,:);
nintetyfifth_kiv=sorted(9,:);

plot(x_vis, f_vis,'LineWidth',2);
hold on;
plot(x_vis,means_kiv,'--r','LineWidth',2);
plot(x_vis,fifth_kiv,':r','LineWidth',2);
plot(x_vis,nintetyfifth_kiv,':r','LineWidth',2);
xlabel('x');
ylabel('y');
hold off;
saveas(gcf,fullfile('../../figures',strcat('KIV_',design,'_',num2str(n_trials))),'epsc');
figure()
hist(results_mse_kiv, 100)
title('MSE for KIV')

%% DualIV - 100 simulations
rng('default')
design='NP'; %NP, SSG, HLLT
[f,sim,x_vis,f_vis]=get_design(design);

% sample size
N=1000;
n_trials=5;
results=zeros(n_trials,N);
results_mse_dualiv=zeros(n_trials,1);

figure()
plot(x_vis,f_vis,'LineWidth',2);
hold on;
for i=1:n_trials
    if mod(i,5)==0
        disp(num2str(i));
    end
    results(i,:)= sim_pred_dualiv(design,N)';
    results_mse_dualiv(i)=mse(results(i,:)',f_vis);
    plot(x_vis,results(i,:),'--r','LineWidth',2);
end

% visualize 100 samples
means_dualiv=mean(results);
sorted=sort(results);
fifth_dualiv=sorted(5,:);
nintetyfifth_dualiv=sorted(4,:);

figure()
plot(x_vis,f_vis,'LineWidth',2);
hold on;
plot(x_vis,means_dualiv,'--r','LineWidth',2);
plot(x_vis,fifth_dualiv,':r','LineWidth',2);
plot(x_vis,nintetyfifth_dualiv,':r','LineWidth',2);
xlabel('x');
ylabel('y');
hold off;
saveas(gcf,fullfile('../figures',strcat('KIV_',design,'_',num2str(n_trials))),'epsc');
figure()
hist(results_mse_dualiv, 100)
title('MSE for DualIV')


%% plot MSEs
mse_mean_dualiv = mean(results_mse_dualiv)
sorted=sort(results_mse_dualiv);
mse_fifth_dualiv=sorted(1);
mse_nintetyfifth_dualiv=sorted(3);



mse_mean_kiv = mean(results_mse_kiv)
sorted=sort(results_mse_kiv);
mse_fifth_kiv=sorted(1);
mse_nintetyfifth_kiv=sorted(3);


figure()
hE = errorbar([1, 2], [mse_mean_dualiv, mse_mean_kiv], [mse_fifth_dualiv, mse_nintetyfifth_dualiv], [mse_fifth_kiv, mse_nintetyfifth_kiv]);
set(hE, 'LineStyle', 'none', 'Marker', '.', 'Color', [.3 .3 .3])
set(hE, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor' , [.7 .7 .7])
xlim([0, 3])
xticks()

set(gca, 'FontSize', 20)

% Add labels
grid on
hTitle = title('Error of structural function estimation');
hXLabel = xlabel('');
hYLabel = ylabel('Mean Squared Error');

% Add text
hText = text(10, 800, ...
    sprintf(''));

% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hText], 'FontSize', 15)

set(hTitle, 'FontSize', 15)

% Adjust axes properties

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:0.03:0.1, ...
    'LineWidth', 1)

h=gca; 
h.XAxis.TickLength = [0 0];

xticks([0, 1, 2, 3])
xticklabels({'', 'DualIV', 'KIV', ''})
set([hYLabel], 'FontSize', 15)
% pbaspect([1 1 1])

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
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print(gcf, fullfile('../figures',strcat('DualIV_',design,'_',num2str(n_trials))), '-dpdf')


