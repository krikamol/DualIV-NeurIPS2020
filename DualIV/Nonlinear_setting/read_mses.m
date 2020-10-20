
n_trials=20;
N_vals=[50, 100, 1000];
rho=0.90;
generic_file_name = 'MSE_N=%d_trials=%d_rho=%s.mat';
all_kivs = cell(length(N_vals), 1);
all_dualivs = cell(length(N_vals), 1);
for i=1:length(N_vals)
    load_path = fullfile('MSE_results',sprintf(generic_file_name, N_vals(i), n_trials, num2str(rho, '%0.2f')))
    sample_size = load(load_path, 'sample_sizes');
    sample_size = sample_size.sample_sizes;
    temp_mses = load(load_path, 'mse_methods');
    temp_mses = temp_mses.mse_methods;
    all_kivs{i} = temp_mses{1};
    all_dualivs{i} = temp_mses{2};
    all_2sls{i} = temp_mses{3};
end
% 
kiv_mses = zeros(n_trials, length(N_vals));
dualiv_mses = zeros(n_trials, length(N_vals));
twosls_mses = zeros(n_trials, length(N_vals));

for i=1:length(N_vals)
    kiv_mses(:, i) = all_kivs{i};
    dualiv_mses(:, i) = all_dualivs{i};
    twosls_mses(:, i) = all_2sls{i};
end


%% DeepIV

n_trials_deepiv=10;

path_to_npy_files = '/Users/amehrjou/Projects/DualIV/Codes/DeepIV/ipynb/MSE_results';
generic_filename = 'MSE_N=%d_trials=%d_rho=%s.npy';
all_deepivs = zeros(n_trials, length(N_vals));
for i=1:length(N_vals)
    load_path = fullfile(path_to_npy_files, sprintf(generic_filename, N_vals(i), n_trials_deepiv, num2str(rho, '%0.2f')));
    temp_mses = readNPY(load_path);
    deepiv_mses(:, i) = log(temp_mses) / log(10);
end



%% Plot
% mse_methods = {kiv_mses, dualiv_mses};
% mse_methods = {kiv_mses, dualiv_mses, deepiv_mses};

% 
% colors = {[1 0 0], [0 0 1], [0 1 1]/2};
% save_path = fullfile('figures',strcat('boxplots_',num2str(rho, '%0.2f'), '.pdf'));
% make_box_plots(mse_methods, N_vals, rho, colors, save_path)

%% Merge with deepIV files
mse_methods = {deepiv_mses, deepiv_mses, deepiv_mses};
colors = {[0 1 1]/2, [0 1 1]/2, [0 1 1]/2};
save_path = fullfile('figures',strcat('boxplots_deepiv_',num2str(rho, '%0.2f'), '.pdf'));
make_box_plots(mse_methods, N_vals, rho, colors, save_path)

% 
% 



