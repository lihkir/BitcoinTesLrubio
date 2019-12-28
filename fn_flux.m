function f = fn_flux(phi, f)

%int: j, N
%double: rho, vrho, p2
%vector: phi, f, delta

global delta;

N = length(phi);

rho = 0;
p2 = 0;
for j = 1:N
    rho = rho + phi(j);
    p2 = p2 + phi(j)*delta(j);
end

vrho = sed_hsf(rho);

for j = 1:N
    f(j) = phi(j)*vrho*(delta(j) - p2);
end