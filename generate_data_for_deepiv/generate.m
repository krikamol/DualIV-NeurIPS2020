clear all;
close all;
clc
% %addpath('../aux'); %test;
% %addpath('../figures');
% %addpath('../results'); %test;
% addpath('../utils');
addpath('../aux2');
warning('off','all');
%% Generative process

design='HLLT'; %NP, SSG, HLLT
[f,sim,x_vis,f_vis]=get_design(design);

% sample size
ntrials = 20;
seeds = [1:20];
sample_sizes = [50, 100, 200, 500, 1000];
rho_values = [0.1, 0.25, 0.5, 0.75 0.9];
mkdir ./datasets
folder = "./datasets/";
for m=1:1:ntrials
    % simulate data, split into stage 1 and 2 samples
    rand('seed', m);
    for rho = rho_values
        for N = sample_sizes
            [x,y,g,z] = sim(f, N, rho);
            dataset_train = [x,y,g,z];
            filename = sprintf("seed_%d_trial_%d_samples_%d_rho_%0.2f.mat", m, m, N, rho);
            save(folder + filename, 'dataset_train')
        end
    end
end
dataset_test = [x_vis,f_vis];
filename = sprintf("test_dataset.mat");
save(folder + filename, 'dataset_test')
