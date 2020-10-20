function y_vis = sim_pred_dualiv(design, N)
%simulate Newey Powell problem, return y_vis for KIV

% simulate data for this design
[f,sim,x_vis,~]=get_design(design);
[x,y,z]=sim(f,N);

vx = median_inter(x);
K_xx = get_K_matrix(x, x, vx);

vy = median_inter(y);
vz = median_inter(z);
yz = [y, z];
vyz = (median_inter_2d(yz))^2;
% L_yy = get_K_matrix(y, y, vy);
% L_zz = get_K_matrix(z, z, vz);
% L_yzyz = times(L_yy, L_zz);


vyz = 20000000;
Vmat = [vy, vyz; vyz, vz];
L_yzyz = get_K_matrix_2d(yz, yz, Vmat);




K = K_xx;
L = L_yzyz;
lambda1 = 0.01; %improves condition number of A for a better inv(A)

gamma = N * norm(L * L, 2) / norm(K * L, 2)^2;
A = L * L + 1 / N * gamma * L * (K * K) * L + lambda1 * eye(N);
% A  = eye(N);


Ainv = inv(A);
lambda2 = 0.01; % reqularized the solution to the following quadratic programming
Q = 2 * K' * L' * Ainv * L * K + lambda2 * eye(N);
R = - 2 * K' * L' * Ainv * L * y;
[beta, feval] = quadprog(Q, R);



vx_vis = median_inter(x_vis);
K_xx_vis = get_K_matrix(x, x_vis, vx_vis);
y_vis = K_xx_vis' * beta;

end

