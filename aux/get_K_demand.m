function df=get_K_demand(x,y,z,x_vis)
% precalculate kernel matrices

p=x(:,1);
t=x(:,2);
s=x(:,3);

vp=median_inter(p);
vt=median_inter(t);
vs=median_inter(s);
vz=median_inter(z);

[x1, x2, y1, y2, z1, z2]=split(x,y,z,.5);

p1=x1(:,1);
t1=x1(:,2);
s1=x1(:,3);

p2=x2(:,1);
t2=x2(:,2);
s2=x2(:,3);

ptest=x_vis(:,1);
ttest=x_vis(:,2);
stest=x_vis(:,3);

df.y1=y1;
df.y2=y2;
df.y=y;

df.K_XX=get_K_matrix(p1,p1,vp).*get_K_matrix(t1,t1,vt).*get_K_matrix(s1,s1,vs);
df.K_xx=get_K_matrix(p2,p2,vp).*get_K_matrix(t2,t2,vt).*get_K_matrix(s2,s2,vs);
df.K_xX=get_K_matrix(p2,p1,vp).*get_K_matrix(t2,t1,vt).*get_K_matrix(s2,s1,vs);
df.K_Xtest=get_K_matrix(p1,ptest,vp).*get_K_matrix(t1,ttest,vt).*get_K_matrix(s1,stest,vs);

df.K_ZZ=get_K_matrix(z1,z1,vz);
df.K_Zz=get_K_matrix(z1,z2,vz);

end

