function u0=rkex(A, b, u0, Ta, cfl, L)

%int:i,n,gc,m,iter,wrt_freq,convec_type
%double:T,cfl,h,cs,dt,t,next_t,ds,cfl1,L,t0,t1
%vector:b,Ta
%matrix:u0,A,u,v

global gc wrt_freq convec_type;


m=size(u0, 1);
n=size(u0, 2);
u=zeros(m, n+2*gc);
u(:, gc+1:n+gc)=u0;

u=bc(u);

h=L/n;
cs=charspeed(u0);
ds=diffus_charspeed(u0);
dt=cfl/(cs/h+2*ds/h^2);
%dt=cfl*h;

t=0;
iter=0;
t0=cputime();
for i=1:length(Ta)
  T=Ta(i);
  while ( t <  T )
    next_t=t+dt;
    if ( next_t > T )
      dt=T-t;
      next_t=T;
    end

    u=do_rkex(A, b, u, h, dt);

    cs=charspeed(u);
    ds=diffus_charspeed(u);
    dt=cfl/(cs/h+2*ds/h^2);
    t=next_t;
    iter=iter+1;
    cfl1=dt*(cs/h+ds/h^2);
    fprintf('%4d % e % e %e %e\n', iter, t, cs, ds, cfl1);
    plot(linspace(0, 20, n), u(:, gc+1:n+gc), 'o-');
    drawnow
  end
  t1=cputime();
  fprintf('CPUTIME T=%.15e => %.15e\n', T, t1-t0);
  u0=u(:, gc+1:n+gc);

  save_data(u0, gc, T);
end

function u=do_rkex(A, b, u, h, dt)
%int:n,gc,i,m,s,j,p,q,l
%double:h,dt
%vector:b
%matrix:A,u,K0,K,ul

global gc;

s=size(A, 2);

n=size(u, 2)-2*gc;
m=size(u, 1);

ul=zeros(m, n+2*gc);
K0=zeros(m, n);

K=zeros(m, n, s);

for l=1:s
  for q=1:n
    for p=1:m
      ul(p, q+gc)=u(p, q+gc);
    end
  end
  
  for j=1:l-1
    for q=1:n
      for p=1:m
        ul(p, q+gc)=ul(p, q+gc)+dt*A(l, j)*K(p, q, j);
      end
    end
  end
  
  ul=bc(ul);
  K0=diffus(ul, h, K0);
  
  for j=1:n
    for i=1:m
      K(i, j, l)=K0(i, j);
    end
  end
  
  K0=convec(ul, h, K0);
  for j=1:n
    for i=1:m
      K(i, j, l)=K(i, j, l) + K0(i, j);
    end
  end
end


for j=1:s
  for q=1:n
    for p=1:m
      u(p, q+gc)=u(p, q+gc)+dt*b(j)*K(p, q, j);
    end
  end
end

u=bc(u);