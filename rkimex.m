function u0 = rkimex(A, b, Ah, bh, u0, Ta, cfl, L, convec_type0, imex_type)

%int:i,n,gc,m,iter,convec_type,convec_type0,wrt_freq,imex_type,test
%double:cfl,h,cs,dt,t,next_t,L,relres,t0,t1,T
%vector:b,bh,Ta 
%matrix:u0,A,Ah,u,v 

% convec_type=0 => PVM
% convec_type=1 => GLF

relres = 0; %% Avoid compiler warnings

global gc convec_type wrt_freq;

convec_type = convec_type0;

m = size(u0, 1);
n = size(u0, 2);

u = zeros(m, n+2*gc);
u(:, gc+1:n+gc) = u0;

u = bc(u);

h = L/n;
cs = charspeed(u);
dt = cfl*h/cs;

if (imex_type == 0)
  fprintf('%4s %-14s %-14s %-14s\n', 'iter', '   t', '   dt', '    cs');
else
  fprintf('%4s %-14s %-14s %-14s %-14s\n', 'iter', '   t', '   dt', '    cs', '||r||/||b|');
end

t = 0;
iter = 0;
t0 = cputime();

for i = 1:length(Ta)
    T = Ta(i);
    while (t <  T)
        next_t = t+dt;
        if (next_t > T)
            dt = T - t;
            next_t = T;
        end
        
        if (imex_type == 0)
            % NI: 100, 1e-10 son iteraciones maximas y tolerancia para resolvedor no lineal
            u = do_rkimex(A, b, Ah, bh, u, h, dt, 100, 1e-10);
        else
            % LI
            [u, relres] = do_lirkimex(A, b, Ah, u, h, dt);
        end
        
        t = next_t;
        iter = iter + 1;
        
        if (imex_type == 0)
            fprintf('iter = %4d | t = %.16f | dt = %.16f | cs = %.16f\n', iter, t, dt, cs);
        else
            fprintf('iter = %4d | t = %.16f | dt = %.16f | cs = %.16f | relres = %.16f\n', iter, t, dt, cs, relres);
        end
        
        cs = charspeed(u);
        dt = cfl*h/cs;
        
        if (wrt_freq > 0 && mod(iter, wrt_freq ) == 0)	
            dump(u);
        end
        
        plot(linspace(0, 1, n), u(:, gc+1:n+gc), 'o-');
        drawnow
        
    end
    t1 = cputime();
    fprintf('CPUTIME T=%.16f => %.16f\n', T, t1-t0);
    
    u0=u(:, gc+1:n+gc);
%   save_data(u0, gc, T);
end

function u = do_rkimex(A, b, Ah, bh, u, h, dt, maxits, tol)

%int:n,gc,i,m,s,j,p,q,l,maxits
%double:h,dt,tol
%vector:b,bh
%matrix:A,Ah,u,K0,K,Kh,B,ul

global gc;

s=size(A, 2);

n=size(u, 2)-2*gc;
m=size(u, 1);

ul=zeros(m, n+2*gc);
K0=zeros(m, n);

K=zeros(m, n, s);
Kh=zeros(m, n, s+1);

K0=convec(u, h, K0);

for j=1:n
  for i=1:m
    Kh(i, j, 1)=K0(i, j);
  end
end
%fprintf('convec 0\n');
%dump(K0, 0);

B=zeros(m, n);%per a allocation en C

for l=1:s
  for q=1:n
    for p=1:m
      B(p, q)=u(p, q+gc);
      for j=1:l-1
	B(p, q)=B(p, q)+dt*A(l, j)*K(p, q, j);
      end
    end
  end
  for j=1:l
    for q=1:n
      for p=1:m
	B(p, q)=B(p, q)+dt*Ah(l+1, j)*Kh(p, q, j);
      end
    end
  end
  ul=u;
  %fprintf('ul before\n');
  %dump(ul, 1);
  
  ul=nwtsolve(dt*A(l,l), B, ul, h, dt, maxits, tol); 
  
  %ul=_fpsolve(dt*A(l,l), B, ul, h, dt, maxits, tol);% afegeix condicions de frontera; bug matlab2cpp, cal llevar-li _ si es lleva %
  %fprintf('ul after\n');
  %dump(ul, 1);
  K0=diffus(ul, h, K0);
  %fprintf('diffus %d', l);
  %dump(K0, 0);
  
  for j=1:n
    for i=1:m
      K(i, j, l)=K0(i, j);
    end
  end
  
  if ( l < s | bh(s+1) ~= 0 )
    K0=convec(ul, h, K0);
    %fprintf('convec %d', l);
    %dump(K0, 0);
    for j=1:n
      for i=1:m
	Kh(i, j, l+1)=K0(i, j);
      end
    end
  end
end


for q=1:n
  for p=1:m
    for j=1:s
      u(p, q+gc)=u(p, q+gc)+dt*(b(j)*K(p, q, j)+bh(j)*Kh(p, q, j));
    end
  end
end


if ( bh(s+1) ~= 0 )
  for q=1:n
    for p=1:m
      u(p, q+gc)=u(p, q+gc)+dt*bh(s+1)*Kh(p, q, s+1);
    end
  end
end

%fprintf('after RKIMEX\n');
%dump(u, 1);

u=bc(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [u, relres]=do_lirkimex(A, b, At, u, h, dt)
%int:n,gc,i,m,s,j,p,q,l
%double:h,dt,relres
%vector:b
%matrix:A,At,u,K0,K,Kh,ut,uh,R

relres=0; %%% avoid compiler warnings
global gc;

s=size(A, 2);

n=size(u, 2)-2*gc;
m=size(u, 1);

K0=zeros(m, n);

K=zeros(m, n, s);


ut=zeros(m, n+2*gc);
uh=zeros(m, n+2*gc);
R=zeros(m, n);

for l=1:s
  for q=1:n
    for p=1:m
      uh(p, q+gc)=u(p, q+gc);
      ut(p, q+gc)=u(p, q+gc);
      for j=1:l-1
	uh(p, q+gc)=uh(p, q+gc)+dt*A (l, j)*K(p, q, j);
	ut(p, q+gc)=ut(p, q+gc)+dt*At(l, j)*K(p, q, j);
      end
    end
  end
  uh=bc(uh);
  ut=bc(ut);
  %%%%%%%%%%% K0=C(ut)
  K0=convec(ut, h, K0);
  %%%%%%%%%%% R=B(ut)*uh
  R=apply_diffus(ut, uh, h, R);
  for q=1:n
    for p=1:m
      K0(p, q)=K0(p, q)+R(p, q);
    end
  end
  
  %%%%%%%%%% K0=C(ut)+B(ut)*uh a la entrada
  [K0, relres]=lsolve(dt*A(l,l), ut, K0, h);
  %%%%%%%%%% K0=inv(I-dt*A(ll)*B(ut))*K0 a la salida
  
  for j=1:n
    for i=1:m
      K(i, j, l)=K0(i, j);
    end
  end
  
end


for q=1:n
  for p=1:m
    for j=1:s
      u(p, q+gc)=u(p, q+gc)+dt*b(j)*K(p, q, j);
    end
  end
end

%fprintf('after RKIMEX\n');
%dump(u, 1);

u=bc(u);