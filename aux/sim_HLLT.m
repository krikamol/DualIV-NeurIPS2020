function [x,y,g,z] = sim_HLLT(f,N,rho)
%SIM_IV generates a data set of size N with according to Hartford, Lewis,
%Leyton-Brown, Taddy with structural function f and confoundedness rho

% simulate (z,t,s,v)

% z=normrnd(0,1,N,1);
% v=normrnd(0,1,N,1);
% s=unidrnd(7,N,1);
% t=rand(N,1).*10;
% e=normrnd(rho.*v,1-rho.^2,N,1);

z = randn(N, 1);
v = randn(N, 1);
s = randi([1, 7], N, 1);
t = rand(N,1) .* 10;

% simulate e from v
e = randn(N, 1) .* 1-rho.^2 + rho.*v;

% simulate p from z,t,v
p= 25 + (z+3) .* get_psi(t) + v;

% simulate y from f
g=f(p,t,s);
y=g+e;
x=[p t s];

end
