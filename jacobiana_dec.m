function Jdec = jacobiana_dec(phi, Jdec)

%int: N, j, i, jac
%double: rho, p2, wrho, w1rho
%vector: phi, delta
%matrix: Jdec 

global  delta;
N=length(phi);

rho = 0;
p2 = 0;
for j = 1:N
  rho = rho+phi(j);
  p2 = p2+phi(j)*delta(j);
end

wrho = sed_hsf(rho);
w1rho = sed_hsf_der(rho);

% J=diag(Jdec(:,1))+Jdec(:, 2:3)*Jdec(:, 4:5)'
for i = 1:N
  Jdec(i, 1) = wrho*(delta(i)-p2);
  Jdec(i, 2) = phi(i)*w1rho*(delta(i)-p2);
  Jdec(i, 3) = phi(i)*wrho;
  Jdec(i, 4) = 1;
  Jdec(i, 5) = -delta(i);
end