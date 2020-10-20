function [beta] = estimate_2sls(x,y,z)
    Pz = z.'*((z*z.')\z);
    beta = ((x*Pz*x.')\x)*Pz*y.';
end

