function [beta] = estimate_dualiv(x,y,z,tol,iter)
    N   = size(x,2);
    w   = [y; z];
    Cw  = (w*w.')./N;
    Cxw = (x*w.')./N;
    Cwx = Cxw.';
    b   = mean(bsxfun(@times,y,w),2);
    M   = Cwx*Cxw;
    reg = norm(Cw, 2) ./ norm(Cxw, 2)^2;
    B   = @(a) Cwx*a;
    C   = @(a) Cxw*a;
    iN  = @(a) (Cw + reg.*M)\a;
    beta = glsqr(B, C, b, tol, iter, [], iN);
end