function fh =numflux_pvm(ul, ur, Sl, Sr, idx_q, fh, A, Al, fl, fr, ua, dif, Qurmul,Qc, S)

%int: i, j, idx_q, int_form, N
%double: ds, s
%vector: S, Qc, Qurmul, fl, fr, fh, Sl, Sr, ua, ul, ur, dif
%matrix: A, Al


global int_form;

N=length(ul);
  

for i = 1:N
    dif(i)=ur(i)-ul(i);
end

if (int_form == 1)
  % Midpoint Rule
  ua=average(ul, ur, 0.5, ua);
  A=jacobiana(ua, A);
elseif (int_form == 5)
  % Midpoint Rule deconstructed Jacobian
  ua=average(ul, ur, 0.5, ua);
  A=jacobiana_dec(ua, A);
elseif (int_form == 2)
  % Gaussian Rule 2 nodes
  
  %0.211324865405187=0.5*(1-1/sqrt(3))
  ua=average(ul, ur, 0.211324865405187, ua);
  Al=jacobiana(ua, Al);
  
  for i = 1:N
    for j = 1:N
      A(i,j)=0.5*Al(i,j);
    end
  end
  %0.788675134594813=0.5*(1+1/sqrt(3))
  ua=average(ul, ur, 0.788675134594813, ua);
  Al=jacobiana(ua, Al);
  
  for i = 1:N
    for j = 1:N
      A(i,j)=A(i, j)+0.5*Al(i,j);
    end
  end
  
elseif (int_form == 3)
  % Gaussian Rule 3 nodes

  ua=average(ul, ur, 0.112701665379258, ua);
  Al=jacobiana(ua, Al);
  
  for i = 1:N
    for j = 1:N
      A(i,j)=5/18*Al(i,j);
    end
  end
  ua=average(ul, ur,  0.887298334620742, ua);
  Al=jacobiana(ua, Al);
  
  for i = 1:N
    for j = 1:N
      A(i,j)=A(i, j)+5/18*Al(i,j);
    end
  end
  ua=average(ul, ur, 0.5, ua);
  Al=jacobiana(ua, Al);
  
  for i = 1:N
    for j = 1:N
      A(i,j)=A(i, j)+8/18*Al(i,j);
    end
  end

elseif (int_form == 4)

  for i = 1:N
    for j = 1:N
      A(i,j)=0;
    end
  end
  ds=1e-4;
  for s=ds/2:ds:1
    ua=average(ul, ur, s, ua);
    Al=jacobiana(ua, Al);
    
    for i = 1:N
      for j = 1:N
	A(i,j)=A(i, j)+Al(i,j);
      end
    end
  end
  for i = 1:N
    for j = 1:N
      A(i,j)=ds*A(i, j);
    end
  end
end

S(1)=min(Sl(1), Sr(1));
S(2)=max(Sl(2), Sr(2));

fl=fn_flux(ul,fl);
fr=fn_flux(ur,fr);

if ( S(1) >= 0)
  for i=1:N
    fh(i)=fl(i);
  end
elseif ( S(2) <= 0)
  for i=1:N
    fh(i)=fr(i);
  end
else
  Qc=Qcoeff(S, idx_q, Qc);
  %pasamos fh para espacio en memoria
  Qurmul=hornerm(A, dif, Qc, fh, Qurmul);        
  for i=1:N
    fh(i)=0.5*(fl(i)+fr(i)-Qurmul(i));
  end
end  