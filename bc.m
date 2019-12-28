function u=bc(u)

%int:n,m,i,j,gc
%matrix:u

global gc;

m=size(u, 1);
n=size(u, 2)-2*gc;

for j=1:gc
  for i=1:m
%     u(i, j)=1e10;
%     u(i, gc+n+j)=1e10;
    u(i, j)=u(i, 2*gc+1-j);
    u(i, gc+n+j)=u(i, gc+1+n-j);
  end
end