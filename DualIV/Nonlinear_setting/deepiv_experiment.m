clear all;
close all;
clc
% addpath('../aux'); %test;
% addpath('../figures');
% addpath('../results'); %test;

%%  Load data from the Generative process

clear all

N = 100;
seed = 6;
pnoise=1.;
psd = 3.7;
pmu = 17.779;
ysd = 158.;
ymu = -292.1;
ynoise=1.;
ypcor = 0.8;

[price, y, z, x] = demand(N, seed, pnoise, psd, pmu, ysd, ymu, ynoise, ypcor);
