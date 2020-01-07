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
        fprintf("\nphit > 1 in the column %d!!\n\n", k);
        error('phit > 1!!');
    end
    
    if (diff_idx ~= 0)
        quot   = 1/phit*(1 - phit);
        sigedd = solid_stress_der_der(phit);
        siged  = solid_stress_der(phit);
        sige   = solid_stress(phit);
        wphit  = sed_hsf(phit);
        wphitd = sed_hsf_der(phit);
        for i = 1:m
            theta = phik(i)*(delta(i) - p2);
            for j = 1:m
                psih = delta(i)*kronecker(i, j) - delta(j)*phik(i) - theta/phit;
                B(i, j, k) = mu_g*((phit*wphitd - wphit)*(theta*siged - psih*sige*quot)/phit^2 + ... 
                             wphit*(theta*sigedd - (theta*sige/phit^2 + psih*((1 - phit)*siged + sige)/(1 - phit))/(1 - phit))/phit);
            end
        end
    else
        for i=1:m
            B(i, i, k) = -D0*nexp*(1-phit)^(nexp-1);
        end
    end
end