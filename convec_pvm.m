function K=convec_pvm(v, h, K)

%int: gc, M, N, i, j, idx_q, N5, n_coef
%double:  cs, csi, h
%matrix: v, K, fh, S, ul, ur, A, Al
%vector: fh1, uli, uri, Sl, Sr, S, Qc, Qurmul, fl, fr, ua, ul, ur, dif, Si

global gc idx_q;

M=size(v, 2)-2*gc;
N=size(v, 1);

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
else
  n_coef=5;
end

S=zeros(2, M+2*gc);
ul=zeros(N, M+2*gc);
ur=zeros(N, M+2*gc);
fh=zeros(N, M+2*gc);

fh1=zeros(N, 1);
uli=zeros(1,N);
uri=zeros(1,N);
Sl=zeros(1,N);
Sr=zeros(1,N);

if ( N== 1 )
    A=zeros(N, N);
else
    N5=max(N, 5);
    A=zeros(N,N5); % Se toma un minimo de 5 columnas por si acaso se pasa decons.
end

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
  elseif ( idx_q == -2)
      fh1= fhroe(uli, uri, fl, fr);
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
    K(j, i-gc)=-(fh(j, i)-fh(j, i-1))/h;
  end
end