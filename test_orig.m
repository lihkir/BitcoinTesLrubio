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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Tests Droplets: Lihki Rubio
% 
% format long e
% 
% global rho_s rho_f g mu_f vinf nexp phimax phi0 wrt_freq beta int_form D0 idx_q phic sigma0 kexp sigedc sigeddc mu_g diff_idx epsilon;
% 
% n = 100;
% T = 50;
% cfl = 0.6;
% L = 0.3;
% convec_type = 0;
% imex_type = 0;
% d1 = 4.96e-4;
% test = 1;
% 
% beta = 1e-15;
% wrt_freq = 1;
% D0 = 1e-7;
% idx_q = 9;
% int_form = 1;
% g = 9.81; 
% mu_f = 0.02416; 
% mu_g = mu_f/g;
% rho_s = 2790; 
% rho_f = 1208; 
% vinf = g*d1*d1*(rho_s-rho_f)/(18*mu_f);
% phi0 = ((nexp-2)*phimax-1)/(nexp-3);
% phimax = 0.66;
% nexp = 4.7;
% phic = 0.2;
% sigma0 = 180; 
% kexp = 2;
% sigedc = sigma0*kexp/phic;
% sigeddc = sigma0*kexp*(kexp - 1)/phic^2;
% epsilon = 1e-5;
% diff_idx = 0;
% 
% u = rk122imex(n, T, cfl, L, convec_type, imex_type, test);
% 