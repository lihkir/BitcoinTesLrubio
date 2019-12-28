function u=fpsolve(a, b, u, h, dt, maxits, tol)
%int:n,m,gc,i,j,its,maxits
%double:a,h,dt,tol,nrm2_b,nrm2_r,nrm2_r0
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
for its=1:maxits
  v=b;
  v=trdsolve(D, L, U, v);
  r=b;
  r=residual(D, L, U, v, r);
  nrm2_r0=0;
  for j=1:n
    for i=1:m
      nrm2_r0=nrm2_r0+r(i, j)*r(i, j);
    end
  end
  [D, L, U]=diffusion_matrix(v, h, a, D, L, U);
  r=b;
  r=residual(D, L, U, v, r);
  nrm2_r=0;
  for j=1:n
    for i=1:m
      nrm2_r=nrm2_r+r(i, j)*r(i, j);
    end
  end
  
  fprintf('\t%3d %e %e\n', its, sqrt(nrm2_r), sqrt(nrm2_r0));
  if ( nrm2_r < tol*tol*nrm2_b)
    fprintf('----------\n');
    break;
  end
end


u(:, gc+1:n+gc)=v;
u=bc(u);


