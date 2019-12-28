function K=convec_glf(v, h, K)
%int:n,m,gc,i,j
%double:h
%matrix:v,K,fh

global gc;

n=size(v, 2)-2*gc;
m=size(v, 1);

fh=zeros(m, n+1);
fh=numflux_glf(v, fh);
for i=1:m
  for j=1:n
    % fh(*, j) => fh en j+2.5 => K(*, j)= Dfh en j+3
    K(i, j)=-(fh(i, j+1)-fh(i, j))/h;
  end
end


function fh=numflux_glf(v, fh)
%int:n,m,gc,i,j,k
%double:cs,fp0,fm0
%vector:f1,fi,vi
%matrix:v,fh,f,fp,fm

global gc;

cs=charspeed(v);

n=size(v, 2)-2*gc;
m=size(v, 1);

fi=zeros(m, 1);
vi=zeros(m, 1);
f=zeros(m, n+2*gc);
for i=1:n+2*gc
  vi=v(:, i);
  fi=fn_flux(vi, fi);
  for j=1:m
    f(j, i)=fi(j);
  end
end


fp=zeros(m, n+2*gc);
fm=zeros(m, n+2*gc);

for i=1:m
  for j=1:n+2*gc
    fp(i, j)=0.5*(f(i, j)+cs*v(i, j));
    fm(i, j)=0.5*(f(i, j)-cs*v(i, j));
  end
end

f1=zeros(5,1);
for k=1:m
  for i=2:n
    for j=1:5
      % stencil: i...i+4 flujo+ en i+2.5
      f1(j)=fp(k, i+j-1); 
    end
    fp0=weno5(f1, 1e-100);

    for j=1:5
      f1(j)=fm(k, i+6-j); % flujo- en i+2.5
    end
    fm0=weno5(f1, 1e-100);
    % fh(*, i) => flujo en i+2.5
    fh(k, i)=fp0+fm0;
  end
end

fh(:, n+1)=0;
fh(:, 1)=0;




