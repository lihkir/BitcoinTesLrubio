function [b, relres]=lsolve(a, u, b, h)
%int:n,m,gc,i,j,its,maxits
%double:a,h,nrm2_b,nrm2_r,relres
%vector:
%matrix:b,u,v,D,L,U,r

global gc;

m=size(u, 1);
n=size(u, 2)-2*gc;

v=zeros(m, n);

v=u(:, gc+1:n+gc);
D=zeros(m, m, n);
L=zeros(m, m, n);
U=zeros(m, m, n);

[D, L, U]=diffusion_matrix(v, h, a, D, L, U);
nrm2_b=0;
for j=1:n
  for i=1:m
    nrm2_b=nrm2_b+b(i, j)*b(i, j);
  end
end

r=zeros(m, n);
v=b;
b=trdsolve(D, L, U, b);
r=v;
r=residual(D, L, U, b, r);
nrm2_r=0;
for j=1:n
  for i=1:m
    nrm2_r=nrm2_r+r(i, j)*r(i, j);
  end
end

relres=sqrt(nrm2_r/nrm2_b);


