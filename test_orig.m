format short e

global D0 vinf nexp delta int_form wrt_freq idx_q imex_type;

wrt_freq=1;
M=50;
T=10;
cfl=0.6;


g = 9.81;
mu_c = 6.5e-3;
rho_d = 1090;
rho_c = 880;
vinf=g*(rho_d-rho_c)/(18*mu_c);
L=0.02;
nexp=4.65;


delta =flipud([ 6.1011e-06
   2.3810e-05
   3.4185e-05
   4.8391e-05
   6.8986e-05
   9.9751e-05
   1.4020e-04
   2.0143e-04
   ]).^2;
u00=flipud([   1.7583e-01
   5.7205e-01
   1.5712e+00
   4.7065e+00
   7.9280e+00
   4.4309e+00
   6.4096e-01
   8.5905e-02
   ])/100;

u0_sed=u00*ones(1,M);


int_form=5;
idx_q=9;


D0=0;
convec_type=0; %PVM
imex_type=0;
disp('================ D0=0 ============================');

um=rkimex_heun(u0_sed, T, cfl, L, convec_type);
disp('==================================================');
u=pvm_muscl(u0_sed, T, cfl, L, idx_q);
disp(norm(um-u, 'fro')/M);
disp('==================================================');


