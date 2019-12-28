function B = diffusion_tensor(u, B)

%int:n,m,i,j,k
%double:phit,phii,D0,nexp
%vector:u
%matrix:u,B

global D0 nexp diff_idx delta mu_g;

m = size(u, 1);
n = size(u, 2);

phik = zeros(m, 1);
for k = 1:n
    
    phit = 0;
    p2 = 0;
    for i = 1:m
        phii = u(i, k);
        phit = phit + phii;
        p2 = p2 + phii*delta(i);
        phik(i) = phii;
    end
    
    if (phit > 1)
        error('phit > 1');
    end
    
    if (diff_idx ~= 0)
        quot  = 1/phit*(1 - phit);
        siged = solid_stress_der(phit);
        sige  = solid_stress(phit);
        wphit = sed_hsf(phit);
        for i = 1:m
            for j = 1:m
                B(i, j, k) = mu_g*wphit*(phik(i)*(delta(i) - p2)*siged - (delta(i)*kronecker(i, j) - delta(j)*phik(i) - phik(i)*quot*(delta(i) - p2))*sige)/phit;
            end
        end
    else
        for i = 1:m
            B(i, i, k) = D0*(1-phit)^nexp;
        end
    end
end