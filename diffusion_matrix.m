function [D, L, U]=diffusion_matrix(u, h, a, D, L, U)

%int:n,m,i,j,p
%double:a,h,f,z
%vector:
%matrix:B,u,D,L,U

n = size(u, 2);
m = size(u, 1);

f = -a*0.5/h^2;
B = zeros(m, m, n);
B = diffusion_tensor(u, B);

for j = 1:n
    for i = 1:m
        for p = 1:m
            D(p, i, j) = 0;
            L(p, i, j) = 0;
            U(p, i, j) = 0;
        end
    end
end

for p = 1:m
    for j = 1:n
        if (j > 1)
            for i = 1:m
                %K(p, j) = K(p, j)-(B(p, i, j)+B(p, i, j-1))*(u(i, j)-u(i, j-1));
                z = B(p, i, j)+B(p, i, j - 1);
                D(p, i, j) = D(p, i, j) - f*z;
                L(p, i, j-1) = L(p, i, j-1) + f*z;
            end
        end
        if (j <  n)
            for i = 1:m
                %K(p, j)=K(p, j)+(B(p, i, j+1)+B(p, i, j))*(u(i, j+1)-u(i, j));
                z = B(p, i, j + 1)+B(p, i, j);
                U(p, i, j) = U(p, i, j) + f*z;
                D(p, i, j) = D(p, i, j) - f*z;
            end
        end
    end
end

for j = 1:n
    for i = 1:m
        D(i, i, j) = D(i, i, j) + 1;
    end
end