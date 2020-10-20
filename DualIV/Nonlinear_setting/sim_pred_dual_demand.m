function mse_value = sim_pred_dual_demand(design, N, rho, grid_search)


[f,sim,x_vis,f_vis]=get_design(design);

% simulate data, split into stage 1 and 2 samples
[x,y,z]=sim(f,N, rho);

p=x(:,1);
t=x(:,2);
s=x(:,3);

vp=median_inter(p);
vt=median_inter(t);
vs=median_inter(s);
vz=median_inter(z);
vy=median_inter(y);

ptest=x_vis(:,1);
ttest=x_vis(:,2);
stest=x_vis(:,3);

K_xx=get_K_matrix(p,p,vp).*get_K_matrix(t,t,vt).*get_K_matrix(s,s,vs);
K_Xtest=get_K_matrix(p,ptest,vp).*get_K_matrix(t,ttest,vt).*get_K_matrix(s,stest,vs);
K_zz = get_K_matrix(z, z, vz);
K_yy = get_K_matrix(y, y, vy);

K = K_xx;
%L = K_yy.*K_zz;

% alternative way of computing L
yz = [y, z];
vyz = 90000;
Vmat = [vy, vyz; vyz, vz];
L_yzyz = get_k_matrix_2d(yz, yz, Vmat);
L = L_yzyz;


if grid_search
    format long
    lambda1_range = logspace(-10, -1, 10); %[1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
    lambda2_range = logspace(-10, -1, 10); %[1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
    best_mse = 1e10;
    for lambda1=lambda1_range
        for lambda2=lambda2_range
            disp(lambda1)
            gamma = N * norm(L * L, 2) / norm(K * L, 2)^2;
            A = L * L + 1 / N * gamma * L * (K * K) * L + lambda1 * eye(N);
            Ainv = inv(A);
            Q = 2 * K' * L' * Ainv * L * K + lambda2 * eye(N);
            R = - 2 * K' * L' * Ainv * L * y;
            [beta, feval] = quadprog(Q, R, [],[],[],[],[],[],[]);
            y_vis = K_Xtest'*beta;
            mse_value = mse(y_vis,f_vis);
            if mse_value < best_mse
                if mse_value < best_mse
                    best_mse = mse_value;
                end
            end
        end
    end
    mse_value = best_mse;
else
    % lambda1 = 0.0000001; %improves condition number of A for a better inv(A)
    gamma = N * norm(L * L, 2) / norm(K * L, 2)^2;
    A = L * L + 1 / N * gamma * L * (K * K) * L + lambda1 * eye(N);
    Ainv = inv(A);
    % lambda2 = 0.0001; % reqularized the solution to the following quadratic programming
    Q = 2 * K' * L' * Ainv * L * K + lambda2 * eye(N);
    R = - 2 * K' * L' * Ainv * L * y;
    [beta, feval] = quadprog(Q, R, [],[],[],[],[],[],[]);

    y_vis = K_Xtest'*beta;
    mse_value = mse(y_vis,f_vis);
end

% % mse
% disp('mse dual method 1:');
% disp(mse(y_vis_dual,f_vis));



