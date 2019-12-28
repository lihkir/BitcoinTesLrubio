function u=prints(u, pad)
%int: i, j, m, n,gc,pad,gc0
%matrix: u

global gc;
m=size(u, 1);
if ( pad == 1 )
  n=size(u, 2)-2*gc;
  gc0=gc;
else
  n=size(u, 2);
  gc0=0;
end


for i=1:m
  for j=1:n
    fprintf('% .15e ', u(i, j+gc0));
  end
  fprintf('\n');
end