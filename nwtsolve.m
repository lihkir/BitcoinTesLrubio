function u=nwtsolve(a, b, u, h, dt, maxits, tol)

%int:n,m,gc,i,j,its,maxits
%double:a,h,dt,tol,nrm2_b,nrm2_r,nrm2_r0,beta,alpha,phi,phi0
%vector:
%matrix:b,u,v,D,L,U,r,r1,du

% volem trobar u tal que B(u)*u=b

global gc beta;

% beta s'introdueix per a evitar que els denominadors de les derivades
% de les pressions de capilaritats es facen cero massa aviat i N-R no convergisca

m=size(u, 1);
n=size(u, 2)-2*gc;

v=zeros(m, n);
v=u(:, gc+1:n+gc);

D=zeros(m, m, n);
L=zeros(m, m, n);
U=zeros(m, m, n);

r=zeros(m, n);
r1=zeros(m, n);
du=zeros(m, n);

[D, L, U]=diffusion_matrix(v, h, a, D, L, U);
r=b;
r=residual(D, L, U, v, r);

% r=b-(L|D|U)*v
nrm2_b=0;
for j=1:n
  for i=1:m
    nrm2_b=nrm2_b+b(i, j)*b(i, j);
  end
end

for its=1:maxits
  [D, L, U]=newton_matrix(v, h, a, D, L, U);
  du=r;
  du=trdsolve(D, L, U, du);
  r1=r;
  r1=residual(D, L, U, du, r1);
  nrm2_r0=0;
  for j=1:n
    for i=1:m
      nrm2_r0=nrm2_r0+r1(i, j)*r1(i, j);
    end
  end

  [D, L, U]=diffusion_matrix(v, h, a, D, L, U);
  r=b;
  r=residual(D, L, U, v, r);
  phi0=0;
  for j=1:n
    for i=1:m
      phi0=phi0+r(i, j)*r(i, j);
    end
  end
  alpha=1;
  phi=2*phi0+1;
  while ( alpha > 1e-4 && phi>=phi0)
    for j=1:n
      for i=1:m
	v(i, j)=v(i, j)+alpha*du(i, j);
      end
    end
    [D, L, U]=diffusion_matrix(v, h, a, D, L, U);
    r=b;
    r=residual(D, L, U, v, r);
    phi=0;
    for j=1:n
      for i=1:m
	phi=phi+r(i, j)*r(i, j);
      end
    end
    alpha=alpha/2;
  end
  
  % calculem r=b-B(v)*v
  [D, L, U]=diffusion_matrix(v, h, a, D, L, U);
  r=b;
  r=residual(D, L, U, v, r);
  nrm2_r=0;
  for j=1:n
    for i=1:m
      nrm2_r=nrm2_r+r(i, j)*r(i, j);
    end
  end
  
  % beta suposa una pertorbacio equivalent al quadrat de l'equacion a resoldre
  % que no impedeix la convergencia quadratica de NR
  beta=nrm2_r;
  if ( nrm2_r > 0 )
    fprintf('\t%3d %e %e\n', its, sqrt(nrm2_r), sqrt(nrm2_r0));
  end
  if ( nrm2_r < tol*tol*nrm2_b)
    if ( nrm2_r > 0 )
      fprintf('----------\n');
    end
    break;
  end
end

u(:, gc+1:n+gc)=v;
u=bc(u);