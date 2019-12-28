function K=apply_diffus(ut, uh, h, K)

%int:n,m,gc,i,j,p
%double:h
%matrix:u,ut,uh,K,B,ut0

global gc;

n=size(uh, 2)-2*gc;
m=size(uh, 1);

u=zeros(m, n);
u=uh(:, gc+1:n+gc);
ut0=zeros(m, n);
ut0=ut(:, gc+1:n+gc);
B=zeros(m, m, n);
B=diffusion_tensor(ut0, B);


for p=1:m
  for j=1:n
    K(p, j)=0;
    if ( j > 1 )
      for i=1:m
	K(p, j)=K(p, j)-(B(p, i, j)+B(p, i, j-1))*(u(i, j)-u(i, j-1));
      end
    end
    if ( j <  n )
      for i=1:m
	K(p, j)=K(p, j)+(B(p, i, j+1)+B(p, i, j))*(u(i, j+1)-u(i, j));
      end
    end
    K(p, j)=0.5/h^2*K(p, j);
  end
end