function [beta] = estimate_kernel_dualiv(x,y,z,kx,ky,kz,tol,iter)
    K   = kx(x,x);
    L   = ky(y,y).*kz(z,z);
    LL  = L*L;
    KL  = K*L;
    b   = L*y;
    M   = KL.'*KL;
    reg = norm(LL, 2)./norm(KL, 2)^2;
    B   = @(a) KL.'*a;
    C   = @(a) KL*a;
    iN  = @(a) (LL + reg.*M)\a;
    beta = glsqr(B, C, b, tol, iter, [], iN);
end