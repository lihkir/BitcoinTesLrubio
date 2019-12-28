function [ul, ur]=muscl(u, ul ,ur)

%int: gc, M, N, i, j
%double: du
%matrix: ul, ur, u

global gc;
M=size(u, 2)-2*gc;
N=size(u, 1);


for i=gc:M+gc+1
  for j=1:N
    du=minmod(u(j, i)-u(j, i-1), u(j, i+1)-u(j, i));
    %ur(i-1)=u(j, i);
    %ul(i)=u(j, i);
    ur(j, i-1)=u(j, i)-0.5*du;
    ul(j, i)=u(j, i)+0.5*du;
  end
end