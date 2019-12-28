function S=lminmax_charspeed(phil, phir, S)

%int:  j, N, i_buf
%double: phi0, phi1, phi2, phi, eta1, eta2, eta, dphi, deta, w1, w2, w, A, B, C, D, alpha,min_0, max_0, min_n, max_n, nexp, phimax
%vector: delta, S, dphi, deta, phil, phir, buf_opt

global delta nexp phimax;

N=length(phil);

phi0=((nexp-2)*phimax-1)/(nexp-3);

phi1=0;
eta1=0;

for j=1:N
  phi1=phi1+phil(j);
  eta1=eta1+delta(j)*phil(j);
end

phi2=0;
eta2=0;
for j=1:N
  phi2=phi2+phir(j);
  eta2=eta2+delta(j)*phir(j);
end

dphi=phi1-phi2;
deta=eta1-eta2;
w1=sed_hsf(phi1);
w2=sed_hsf(phi2);

%fprintf('vinf=%e phimax=%e phi0=%e nexp=%e phi1=%e w1=%e phi2=%e w2=%e\n', vinf, phimax, phi0, nexp, phi1, w1, phi2, w2);

i_buf=1;

buf_opt=zeros(10,1);
buf_opt(i_buf)=w1*(delta(1)-eta1); i_buf=i_buf+1;
buf_opt(i_buf)=w2*(delta(1)-eta2); i_buf=i_buf+1;

if ( dphi ~= 0 && deta ~= 0 ) 
  % Hay minimo local en segmento con phi < phi0 ?
  alpha=((nexp-1)*(delta(1)-eta2)*dphi+(1-phi2)*deta)/(nexp*dphi*deta);
  if ( alpha > 0 && alpha < 1 )
    phi=phi2+alpha*dphi;
    if ( phi<phi0 )
      eta=eta2+alpha*deta;
      buf_opt(i_buf)=sed_hsf(phi)*(delta(1)-eta);  i_buf=i_buf+1;
    end
  end
  
  A=3*dphi*dphi*deta;
  B=-(2*dphi*deta*(1+phimax-2*phi2)+2*dphi*dphi*(delta(1)-eta2));
  C=dphi*(1+phimax-2*phi2)*(delta(1)-eta2)-(phi2-phimax)*(1-phi2)*deta;
  % Hay minimo local en segmento con phi >= phi0 ?
  D=B*B-4*A*C;

  if ( D >= 0 )
    alpha=(-B+sqrt(D))/(2*A);
    
    if ( alpha > 0 && alpha < 1 )
      phi=phi2+alpha*dphi;
      if ( phi >= phi0 )
	eta=eta2+alpha*deta;
	buf_opt(i_buf)=sed_hsf(phi)*(delta(1)-eta);  i_buf=i_buf+1;
      end
    end
    
    alpha=(-B-sqrt(D))/(2*A);
    
    if ( alpha > 0 && alpha < 1 )
      phi=phi2+alpha*dphi;
      if ( phi >= phi0 )
	eta=eta2+alpha*deta;
	buf_opt(i_buf)=sed_hsf(phi)*(delta(1)-eta);  i_buf=i_buf+1;
      end
    end
  end
end
%  singularidad en segmento ? 
alpha=-(phi2-phi0)/dphi;
if ( alpha > 0 && alpha < 1 ) 
  phi=phi2+alpha*dphi;
  eta=eta2+alpha*deta;
  buf_opt(i_buf)=sed_hsf(phi)*(delta(1)-eta);  i_buf=i_buf+1;
end

[min_0, max_0]=find_extrema(buf_opt, i_buf-1);

i_buf=1;
% Extremos segmento  
buf_opt(i_buf)=w1*delta(N)+eta1*(sed_hsf_der(phi1)*(1-phi1)-2*w1);
i_buf=i_buf+1;
buf_opt(i_buf)=w2*delta(N)+eta2*(sed_hsf_der(phi2)*(1-phi2)-2*w2);
i_buf=i_buf+1;
if ( dphi ~= 0 && deta ~= 0 ) 
  % Hay minimo local en segmento con phi < phi0 ?
  alpha=((nexp-1)*(delta(N)-(nexp+1)*eta2)*dphi+(1-phi2)*(nexp+1)*deta)/(nexp*(nexp+1)*dphi*deta);
  if ( alpha > 0 && alpha < 1  ) 
    phi=phi2+alpha*dphi;
    if  ( phi<phi0 )
      eta=eta2+alpha*deta;
      w=sed_hsf(phi);
      buf_opt(i_buf)=w*delta(N)+eta*(sed_hsf_der(phi)*(1-phi)-2*w);
      i_buf=i_buf+1;
    end
  end
  
  % Hay minimo local en segmento con phi >= phi0 ?
  C=deta *(1- 5*phi2+ 3*phimax+ 4*phi2*phi2- 3*phi2*phimax)+  dphi*(delta(N)- 5*eta2- 2*delta(N)*phi2+ delta(N)*phimax+ 8*eta2*phi2- 3*eta2*phimax); 
  B=2*dphi*(dphi*(4*eta2 - delta(N)) +deta*(- 5 + 8*phi2 - 3*phimax));
  A=12*deta*dphi*dphi;
  D=B*B-4*A*C;
  if ( D >= 0 )
    alpha=(-B+sqrt(D))/(2*A);
    
    if ( alpha > 0 && alpha < 1 )
      phi=phi2+alpha*dphi;
      if ( phi >= phi0)
	eta=eta2+alpha*deta;
	w=sed_hsf(phi);
	buf_opt(i_buf)=w*delta(N)+eta*(sed_hsf_der(phi)*(1-phi)-2*w);
	i_buf=i_buf+1;
      end
    end
    alpha=(-B-sqrt(D))/(2*A);
    if ( alpha > 0 && alpha < 1 )
      phi=phi2+alpha*dphi;
      if ( phi >= phi0)
	eta=eta2+alpha*deta;
	w=sed_hsf(phi);
	buf_opt(i_buf)=w*delta(N)+eta*(sed_hsf_der(phi)*(1-phi)-2*w);
	i_buf=i_buf+1;
      end
    end
  end
end
%  singularidad en segmento ? 
alpha=-(phi2-phi0)/dphi;
if ( alpha > 0 && alpha < 1 ) 
  phi=phi2+alpha*dphi;
  eta=eta2+alpha*deta;
  w=sed_hsf(phi);
  buf_opt(i_buf)=w*delta(N)+eta*(sed_hsf_der(phi)*(1-phi)-2*w);
  i_buf=i_buf+1;
end
[min_n, max_n]=find_extrema(buf_opt, i_buf-1);

S(1)=min_n;
S(2)=max_0;
  
function [vmin, vmax]=find_extrema(v, l)
%vector: v
%int: i, l
%double: vmin, vmax

vmin=v(1);
vmax=v(1);

for i=2:l
  vmin=min(vmin, v(i));
  vmax=max(vmax, v(i));
end