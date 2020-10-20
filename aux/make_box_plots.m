% % 
% % %% Test
% clear all
% close all
% clc
% number_of_sample_sizes = 5;
% sample_sizes = [1e2, 2e2, 5e2, 1e3, 1e4];
% samples_per_setting = 100;
% rho = 0.2;
% design = 'demand';
% n_trials = 10;
% mse_method1 = normrnd(5,0.3,samples_per_setting,number_of_sample_sizes);
% mse_method2 = normrnd(7,0.6,samples_per_setting,number_of_sample_sizes);
% mse_method3 = normrnd(3.5,0.5,samples_per_setting,number_of_sample_sizes);
% colors = {[1 0 0], [0 0 1], [1 0 1]};
% mse_methods = {mse_method1, mse_method2, mse_method3}
% save_path = fullfile('../figures',strcat('DualIV_',design,'_',num2str(n_trials)))
% 
% make_box_plots1(mse_methods, sample_sizes, rho, colors, save_path)

%% Actual function
function a = make_box_plots(mse_methods, sample_sizes, rho, colors, save_path)
n_trials = length(mse_methods{1});
figure()
hold on
for i=1:1:length(mse_methods)
    bp = boxplot(mse_methods{i}, sample_sizes,...
        'BoxStyle', 'outline', 'Colors', colors{i},...
        'Notch', 'off', 'MedianStyle', 'line',...
        'Widths', 0.5, 'Symbol','o',...
        'OutlierSize',1)
    set(bp, 'LineWidth', 2)
%     eb = errorbar(1, mean(mse_methods{i}), - std(mse_methods{i}) / sqrt(n_trials), std(mse_methods{i}) / sqrt(n_trials), 's');
%     set(eb, 'LineWidth', 2, 'Color', colors{i})
end    

% xlim([0, 3])
ylim([2, 5]);
% xticks()

set(gca, 'FontSize', 20)

% Add labels
grid on
bpTitle = title(sprintf('$\\rho=%0.2f$', rho), 'Interpreter','latex');
bpXLabel = xlabel('sample size');
% bpYLabel = ylabel('Mean Squared Error');

% Adjust font
set(gca, 'FontName', 'Helvetica')
set([bpXLabel], 'FontSize', 20)

set(bpTitle, 'FontSize', 20)

% Adjust axes properties

set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gca,'box','off')
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

print(gcf, save_path, '-dpdf')

end





