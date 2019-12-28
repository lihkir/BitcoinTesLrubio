function u=rk122imex(n, T, cfl, L, convec_type, imex_type, test)

%int:n,gc,i
%double:cfl,D0,T,beta,nexp
%vector:b,bh,delta
%matrix:u0,A,Ah,u

global delta gc;

if (convec_type == 0)
  %MUSCL - PVM
  gc = 2;
else
  % WENO5 - COMPGLF
  gc = 3;
end

if (test == 1)
    m = 3;
else
    quit
end

u0 = zeros(m, n);
        
if (test == 1)
    delta(1) = 1;
    delta(2) = 0.5;
    delta(3) = 0.25;
    for j = 1:n
        u0(:, j) = 0.04;
    end
else
    error("Undefined test");
end

A = zeros(1,1);
b = zeros(1,1);
A(1,1) = 0.5;
b(1) = 1;

Ah = zeros(2,2);
bh = zeros(1,2);
Ah(2,1) = 1;
bh(1) = 0.5;
bh(2) = 0.5;

u = rkimex(A, b, Ah, bh, u0, T, cfl, L, convec_type, imex_type);      