function cs=charspeed(u)

%int:i,j,M,gc
%double:cs,csi
%matrix:u,S

global gc;

M=size(u, 2)-2*gc;

S=zeros(2, M+2*gc);
S=minmax_charspeed(u, S);

cs=0;
for i=gc+1:M+gc
  csi=max(abs(S(1, i)), abs(S(2, i)));
  cs=max(cs, csi);
end


