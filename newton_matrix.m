function [D, L, U]=newton_matrix(u, h, a, D, L, U)

%int:n,m,i,j,p,q,r
%double:a,h,f,z
%vector:
%matrix:B,Bp,u,D,L,U

n=size(u, 2);
m=size(u, 1);

f=-a*0.5/h^2;
B=zeros(m, m, n);
Bp=zeros(m, m, n);
B=diffusion_tensor(u, B);

for j=1:n
  for i=1:m
    for p=1:m
      D(p, i, j)=0;
      L(p, i, j)=0;
      U(p, i, j)=0;
    end
  end
end

for r=1:m
  Bp=der_diffusion_tensor(u, Bp, r);
  for p=1:m
    for i=2:n
      z=B(p, r, i)+B(p, r, i-1);
      D(p, r, i)=D(p, r, i)-f*z;
      L(p, r, i-1)=L(p, r, i-1)+f*z;
      for q=1:m
	D(p, r, i)=D(p, r, i)-f*Bp(p, q, i)*(u(q, i)-u(q, i-1));
      end
      for q=1:m
	L(p, r, i-1)=L(p, r, i-1)-f*Bp(p, q, i-1)*(u(q, i)-u(q, i-1));
      end
    end
    for i=1:n-1
      z=B(p, r, i+1)+B(p, r, i);
      U(p, r, i)=U(p, r, i)+f*z;
      D(p, r, i)=D(p, r, i)-f*z;
      for q=1:m
	D(p, r, i)=D(p, r, i)+f*Bp(p, q, i)*(u(q, i+1)-u(q, i));
      end
      for q=1:m
	U(p, r, i)=U(p, r, i)+f*Bp(p, q, i+1)*(u(q, i+1)-u(q, i));
      end

    end
  end
end
for i=1:n
  for r=1:m
    D(r, r, i)=D(r, r, i)+1;
  end
end