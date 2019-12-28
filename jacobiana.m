function jac = jacobiana(phi, jac)

%int: N, j, i
%double: rho, p2, u, v
%vector: phi, delta
%matrix: jac 

global  delta;

N = length(phi);

rho = 0;
p2 = 0;
for j = 1:N
    rho = rho + phi(j);
    p2 = p2 + phi(j)*delta(j);
end

u = sed_hsf(rho);
v = sed_hsf_der(rho); 

for i = 1:N
    for j = 1:N
        jac(i, j) = phi(i)*(v*(delta(i) - p2) - u*delta(j));
    end
end
    
for i = 1:N
    jac(i, i) = jac(i, i) + u*(delta(i) - p2);
end