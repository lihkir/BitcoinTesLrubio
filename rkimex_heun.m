function u=rkimex_heun(u0, T, cfl, L, convec_type)
%int:n,gc,i
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

  

A=zeros(1,1);
b=zeros(1,1);
A(1,1)=0.5;
b(1)=1;

Ah=zeros(2,2);
bh=zeros(1,2);
Ah(2,1)=1;
bh(1)=0.5;
bh(2)=0.5;

%for i=1:length(delta)
%  fprintf('delta(%d)=%.15e\n', i, delta(i));
%end

u=rkimex(A, b, Ah, bh, u0, T, cfl, L, convec_type);      
      


