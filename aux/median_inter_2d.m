function vyz = median_inter_2d(yz)
%median interpoint distance
y = yz(:, 1);
z = yz(:, 2);
n=length(y);

A=repmat(y,1,n);
B=A';
dist_y=abs(A-B);
dist_y=reshape(dist_y,[],1);
% vy=median(dist_y);

A=repmat(z,1,n);
B=A';
dist_z=abs(A-B);
dist_z=reshape(dist_z,[],1);

vyz = median(dist_y .* dist_z);
end
