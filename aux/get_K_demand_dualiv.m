function df=get_K_demand_dualiv(x,y,z,x_vis)
% precalculate kernel matrices

p=x(:,1);
t=x(:,2);
s=x(:,3);

vp=median_inter(p);
vt=median_inter(t);
vs=median_inter(s);
vz=median_inter(z);


ptest=x_vis(:,1);
ttest=x_vis(:,2);
stest=x_vis(:,3);

df.y=y;

df.K_xx=get_K_matrix(p,p,vp).*get_K_matrix(t,t,vt).*get_K_matrix(s,s,vs);
df.K_Xtest=get_K_matrix(p,ptest,vp).*get_K_matrix(t,ttest,vt).*get_K_matrix(s,stest,vs);
