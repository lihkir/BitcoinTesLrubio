function p=matmult_dec(Adec, x, p)
% calcula p=p+A*x, A=diag(Adec(:,1))+Adec(:, 2)*Adec(:, 3)'

%int: i, j, k, N, K
%matrix: Adec, x, p
%double: s

N=size(x, 1);
K=size(x, 2);

if ( N== 1 )
    p=matmult(Adec, x, p);
else
for k=1:K
  % p=p+diag(Adec(:,1))*x
  for i=1:N
    p(i, k) = p(i, k)+Adec(i,1)*x(i, k);
  end
  
  % s=sum(x); p=p+Adec(:,2)*s
  s=0;
  for i=1:N
    s=s+x(i, k);
  end
  for i=1:N
    p(i, k)=p(i, k)+Adec(i, 2)*s;
  end
end
end