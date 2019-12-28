%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests Sedimentation: Lihki Rubio

format long e

global rho_s g mu_f vinf nexp phimax phi0 wrt_freq beta int_form D0 idx_q phic sigma0 kexp sigedc sigeddc diff_idx epsilon mu_g;

n = 100;
T = 0.07;
cfl = 0.6;
L = 1;
convec_type = 0;
imex_type = 0;
d1 = 4.96e-4;
test = 1;

beta = 1e-15;
wrt_freq = 1;
D0 = 1e-7;
idx_q = 9;
int_form = 1;
g = 9.81; 

mu_f = 10e-3; 
rho_s = 1800; 
vinf = mu_f*rho_s;
mu_g = 1/g*rho_s;

phi0 = ((nexp-2)*phimax-1)/(nexp-3);
phimax = 0.66;
nexp = 4.7;
phic = 0.2;
sigma0 = 180; 
kexp = 2;
sigedc = sigma0*kexp/phic;
sigeddc = sigma0*kexp*(kexp - 1)/phic^2;
epsilon = 1e-5;
diff_idx = 0;

u = rk122imex(n, T, cfl, L, convec_type, imex_type, test);

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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%