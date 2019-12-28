function p=matmult(A, x, p)
% calcula p=p+A*x

%int: i, j, k, N, K
%matrix: A, x, p

N=size(x, 1);
K=size(x, 2);

for k=1:K
  for i=1:N
    for j=1:N
      p(i, k) = p(i, k)+A(i, j)*x(j, k);
    end
  end
end