function u=rk2ex(n, T, cfl, L, convec_type, test)

%int:n,gc,i,convec_type0,convec_type
%double:cfl,D0,T,beta,nexp
%vector:b,bh,delta
%matrix:u0,A,Ah,u

global vinf nexp delta wrt_freq gc beta int_form D0 idx_q;

beta=1e-15;
wrt_freq=1;
D0=0;
idx_q=9;

if ( convec_type == 0 )
  %MUSCL - PVM
  gc=2;
else
  % WENO5 - COMPGLF
  gc=3;
end

g = 9.81;
mu_c = 6.5e-3;
rho_d = 1090;
rho_c = 880;
vinf=g*(rho_d-rho_c)/(18*mu_c);
nexp=4.65;
int_form=5;

if ( test == 3 )
   delta(1)=1.745725584039578e-07;
   delta(2)=8.502483987726233e-08;
   delta(3)=4.114961917286227e-08;
   delta(4)=2.055894601643213e-08;
   delta(5)=1.002369639342369e-08;
   delta(6)=4.709986484700493e-09;
   delta(7)=2.328921953264102e-09;
   delta(8)=1.148277502313191e-09;
   u0=zeros(8, n);
   for j=1:n
       u0(1, j)=3.288742053268562e-03;
       u0(2, j)=1.138015958334105e-01;
       u0(3, j)=2.500988514694888e-01;
       u0(4, j)=9.920555854169366e-02;
       u0(5, j)=2.305478842696419e-02;
       u0(6, j)=8.206216521787478e-03;
       u0(7, j)=5.020573462087361e-03;
       u0(8, j)=1.834930402387243e-03;
   end
elseif ( test == 4 )
    delta=zeros(1,1);
    delta(1)=5e-5;
    u0=zeros(1, n);
    for j=1:n
        u0(1, j)=0.2;
    end
else
    quit
end

A=zeros(2,2);
b=zeros(1,2);
A(2,1)=0.5;
b(1)=0;
b(2)=1;

u=rkex(A, b, u0, T, cfl, L);      