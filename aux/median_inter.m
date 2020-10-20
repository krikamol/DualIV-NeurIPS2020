function vx = median_inter(x)
%median interpoint distance

n=length(x);
A=repmat(x,1,n);
B=A';
dist=abs(A-B);
dist=reshape(dist,[],1);
vx=median(dist);

end