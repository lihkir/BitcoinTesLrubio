function  r=residual(D, L, U, v, r)
%int:n,m,i,j,k
%matrix:v,D,L,U,r


m=size(r, 1);
n=size(D, 3);

for i=1:m
  for k=1:m
    r(i, 1)=r(i, 1)-D(i, k, 1)*v(k, 1) - U(i, k, 1)*v(k, 2);
  end
  for j=2:n-1
    for k=1:m
      r(i,j)=r(i,j)-L(i,k,j-1)*v(k,j-1)-D(i,k,j)*v(k,j)-U(i,k,j)*v(k,j+1);
    end
  end
  for k=1:m
    r(i, n)=r(i, n)-L(i, k, n-1)*v(k, n-1)-D(i, k, n)*v(k, n);
  end
end

