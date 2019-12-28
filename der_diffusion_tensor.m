function B=der_diffusion_tensor(u, B, r)

%int:n,m,i,j,k,r
%double:phit,phii,D0,nexp
%vector:u
%matrix:u,B

global D0 nexp diff_idx delta mu_g;

n = size(u, 2);
m = size(u, 1);

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
        quot1  = 1/phit;
        quot2  = quot1^2;
        quot3  = 1/phit*(1 - phit);
        quot4  = quot3^2;
        sigedd = solid_stress_der_der(phit);
        siged  = solid_stress_der(phit);
        sige   = solid_stress(phit);
        wphit  = sed_hsf(phit);
        wphitd = sed_hsf_der(phit);
        for i = 1:m
            for j = 1:m
                B(i, j, k) = mu_g*((phit*wphitd - wphit)*quot2*(phik(i)*(delta(i) - p2)*siged - (delta(i)*kronecker(i, j) - delta(j)*phik(i) - phik(i)*quot3*(delta(i) - p2))*sige) + wphit*quot1*( phik(i)*(delta(i) - p2)*sigedd - ((1 - 2*phit)*quot4*phik(i)*(delta(i) - p2)*sige + siged*(delta(i)*kronecker(i, j) - delta(j)*phik(i) - phik(i)*(delta(i) - p2)*quot3))));
            end
        end
    else
        for i=1:m
            B(i, i, k) = -D0*nexp*(1-phit)^(nexp-1);
        end
    end
end