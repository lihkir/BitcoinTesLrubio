function u0=pvm_muscl(u0, T, cfl, L, idx_q)

%int: i, j, gc, M, N, iter, idx_q, n_coef, wrt_freq
%double: h, t, next_t, T, cfl,lambda, cs, alpha, dt, L
%matrix: u0, u, v, un, fh, S, ul, ur

global gc n_coef wrt_freq;

if (idx_q == 2)
  n_coef=1;
elseif ( idx_q == 3 )
  n_coef=2;  
elseif ( idx_q == 4 )
  n_coef=3;  
elseif ( idx_q == 5 )
  n_coef=3;  
elseif ( idx_q == 6 )
  n_coef=5;  
elseif ( idx_q == 7 )
  n_coef=5;
end
M=size(u0, 2);
N=size(u0, 1);

gc=2;

S=zeros(2, M+2*gc);
ul=zeros(N, M+2*gc);
ur=zeros(N, M+2*gc);
u=zeros(N, M+2*gc);
un=zeros(N, M+2*gc);
v=zeros(N, M+2*gc);
fh=zeros(N, M+2*gc);

u(:, gc+1:M+gc)=u0;

h=L/M;

t=0;
iter=1;
u=bc(u);
cs=soundspeed(u);
dump(u, 0);


while ( t < T )
  alpha=cs;
  dt=h*cfl/cs;
  next_t=t+dt;
  if ( next_t >= T )
    next_t=T;
    dt=T-t;
  end

  lambda=dt/h;
  un=u;
  v=u;
  cs0=cs;
  [v, cs, fh]=euler(v, alpha, lambda, cs, fh, ul, ur, S, idx_q);
  [v, cs, fh]=euler(v, alpha, lambda, cs, fh, ul, ur, S, idx_q);
  for i=gc+1:M+gc
    for j=1:N 
      u(j, i)=0.5*(un(j, i)+v(j, i));
    end
  end
  cs=charspeed(u);
  fprintf('%4d %e %e %e\n', iter, next_t, dt, cs0);
  
  if ( wrt_freq > 0 && mod(iter, wrt_freq ) == 0 )	
    dump(u, iter);
  end
  t=next_t;
  iter=iter+1;
  %cla
  %plot(u(gc+1:M+gc, :), '+-');
  %hold on
  %plot(sum(u(gc+1:M+gc, :)'), 'k+-');
  %drawnow
end

u0=u(:, gc+1:M+gc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v, cs, fh]=euler(v, alpha, lambda, cs, fh, ul, ur, S, idx_q)

%int: gc, M, N, i, j, idx_q, N5, n_coef
%double: alpha, lambda, cs, csi
%matrix: v, fh, S, ul, ur, A, Al
%vector: fh1, uli, uri, Sl, Sr, S, Qc, Qurmul, fl, fr, ua, ul, ur, dif, Si

global gc n_coef;
M=size(v, 2)-2*gc;
N=size(v, 1);


fh1=zeros(N, 1);
uli=zeros(1,N);
uri=zeros(1,N);
Sl=zeros(1,N);
Sr=zeros(1,N);

N5=max(N, 5);
A=zeros(N,N5); % Se toma un minimo de 5 columnas por si acaso se pasa decons.
Al=zeros(N,N); 

fl=zeros(1,N);
fr=zeros(1,N);
ua=zeros(1,N);
dif=zeros(N,1);
Qurmul=zeros(N,1);
Qc=zeros(1, n_coef);
Si=zeros(1,2);

v=bc(v);

S=minmax_charspeed(v, S);

cs=0;
for i=gc+1:M+gc
  csi=max(abs(S(1, i)), abs(S(2, i)));
  cs=max(cs, csi);
end

[ul, ur]=muscl(v, ul ,ur);

for i=gc+1:M+gc-1
  for j = 1:N
    uli(j) = ul(j, i);
    uri(j) = ur(j, i);
  end
  for j = 1:2
    Sl(j) = S(j, i);
    Sr(j) = S(j, i+1);
  end

  if ( idx_q > 0 )
    fh1=numflux_pvm(uli, uri, Sl, Sr, idx_q, fh1, A, Al, fl, fr, ua, dif, Qurmul, Qc, Si);
  else
    fh1=numflux_hll(uli, uri, Sl, Sr, fh1, fl, fr, Si);
  end  
  for j=1:N
    fh(j, i)=fh1(j);
  end    
end
for j=1:N
  fh(j, gc)=0;
  fh(j, M+gc)=0;
end    

for i=gc+1:M+gc
  for j=1:N
    v(j, i)=v(j, i)-lambda*(fh(j, i)-fh(j, i-1));
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cs=soundspeed(v)

%int: i, j, gc, M, N
%double: cs, csi
%vector: vi, vip1, S
%matrix: v

global gc;
M=size(v, 2)-2*gc;
N=size(v, 1);


cs=0;
vi=zeros(N, 1);
vip1=zeros(N, 1);
S=zeros(2, 1);

for i=1:M+2*gc-1
  for j=1:N
    vi(j)=v(j, i);
  end
  for j=1:N
    vip1(j)=v(j, i+1);
  end
  S=lminmax_charspeed(vi, vip1, S);
  csi=max(abs(S(1)), abs(S(2)));
  cs=max(cs, csi);
end

