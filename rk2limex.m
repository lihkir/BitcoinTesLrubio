function u=rk2limex(u0, T, cfl, L, convec_type)
%int:n,gc,i,imex_type
%double:cfl,D0,T,beta,nexp
%vector:b,bh,delta
%matrix:u0,A,Ah,u

global gc beta delta;

beta=1e-15;

if ( convec_type == 0 )
  %MUSCL - PVM
  gc=2;
else
  % WENO5 - COMPGLF
  gc=3;
end


A=zeros(2,2);
A(2,1)=0.5;
A(2,2)=0.5;

Ah=zeros(2,2);
Ah(2,1)=1;

b=zeros(2,1);
b(1)=0.5;
b(2)=0.5;

imex_type=1;

u=rkimex(A, b, Ah, b, u0, T, cfl, L, convec_type, imex_type);      
      


