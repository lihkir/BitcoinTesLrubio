function K=diffus(u1, h, K)
%int:n,m,gc,i,j,p
%double:h
%matrix:u,u1,K,B


global gc;

n=size(u1, 2)-2*gc;
m=size(u1, 1);

u=zeros(m, n);
u=u1(:, gc+1:n+gc);
B=zeros(m, m, n);
B=diffusion_tensor(u, B);


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



