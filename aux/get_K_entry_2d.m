function K_entry = get_K_entry_2d(x,z, Vmat)
%GET_K_ENTRY calculates an entry of the K matrix for inputs x and z. It
%   uses a radial basis function kernel with matrix bandwidth Vmat

K_entry=exp((x-z)' * Vmat * (x-z)/2);
    
end

