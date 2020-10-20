function [price, y, z, x_latent] = simulate_demand(f, N, pnoise, psd, pmu, ysd, ymu, ynoise, ypcor)
% This function generate data from the demand experiment of DeepIV paper.
% Author: Arash Mehrjou
% Last update: 4 Oct 2019
% Returns [treatment, outcome, instrument, covariate]

% covariates: time and emotion
n_emotions = 7;
n = 100;
time = rand(n, 1) * 10;
emotion_id = randi([1, 7], n, 1);
one_hot_maker = @(i) full(ind2vec(i, n_emotions));
emotion_feature = one_hot_maker(emotion_id')';

% random instrument
z = randn(n, 1);

% z -> price
v = randn(n, 1) * pnoise;
price = sensf(time) .* (z + 3)  + 25.;
price = price + v;
price = (price - pmu) / psd;

% true observable demand function
x = cat(2, time, emotion_feature);
x_latent = cat(2, time, emotion_id);
% g = lambda x, z, p: storeg(x, p) % doesn't use z

% errors 
e = (ypcor * ynoise / pnoise) * v + randn(n, 1) * ynoise * sqrt(1 - ypcor ^ 2);

% response
y = f(x_latent, price, psd, pmu, ymu, ysd) + e;

