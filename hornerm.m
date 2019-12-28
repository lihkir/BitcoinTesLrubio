function pAx=hornerm(A, x, a, p, pAx)

%int: k, n, N, i, j, int_form, K
%vector: a
%matrix: A, pAx, x, p 

global int_form;

N=length(x);

K=size(x, 2);
n=length(a); 

for k=1:K
  for i=1:N
    pAx(i, k)=a(n)*x(i, k);
  end
end

for k=1:K
  for j=n-1:-1:1
    % p=a*x
    if ( a(j) ~= 0 )
      for i = 1:N
	p(i, k)=a(j)*x(i, k);
      end
    else
      for i = 1:N
	p(i, k)=0;
      end
    end
    
    %p=p+A*pAx = a*x+A*pAx
    if ( int_form ~= 5 )
      p=matmult(A, pAx, p);
    else
      p=matmult_dec(A, pAx, p);
    end
    %pAx= p = a*x+A*pAx
    for i = 1:N
      pAx(i, k)=p(i, k);
    end
  end
end