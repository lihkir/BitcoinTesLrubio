function fh =numflux_hll(ul, ur, Sl, Sr, fh, fl, fr, S)

%int: i, j, N
%double:  gp, gm, ids
%vector: fl, fr, fh, Sl, Sr, ul, ur, S


N=length(ul);
  


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
  ids=1/(S(2)-S(1));
  gp=S(2)*ids;
  gm=-S(1)*ids;
  for i=1:N
    fh(i)=gp*(fl(i)-S(1)*ul(i))+gm*(fr(i)-S(2)*ur(i));
  end
end
  