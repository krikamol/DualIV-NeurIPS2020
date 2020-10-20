function [X,Y,Z] = simulate(N, beta, g, d)
%% [x,y,z] = simulate(N) generates N i.i.d samples from the linear model
%  Input:  N - the number of samples 
%          b - the true structural parameter
%          g - the gamma parameter controlling the strength of confounder 
%          d - the number of instruments
%
%  Output: [X,Y,Z] - matrices of treatments, outcomes, and instruments

e   = normrnd(0,2,[1,N]);
ep  = normrnd(0,0.15,[1,N]);
eta = normrnd(0,0.15,[1,N]);

Z   = normrnd(0,2,[d,N]);
X   = (1.0-g).*sin(Z(1,:)) + g.*e + eta;
Y   = beta.*X + e + ep;

end