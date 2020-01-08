%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tests Sedimentation: Lihki Rubio

clear
format short

%%
global rho_s g mu_f vinf nexp phimax phi0 wrt_freq beta int_form D0 idx_q phic sigma0 kexp sigedc sigeddc diff_idx epsilon mu_g;

n = 200;
T = 10000;
cfl = 0.6;
L = 1;
convec_type = 0;
imex_type = 1;
d1 = 1.19e-5;
test = 1;
%%

beta = 1e-15;
wrt_freq = 1;
D0 = 1e-7;
idx_q = 9;
int_form = 1;
g = 9.81; 
mu_f = 1e-3; 
rho_s = 1800; 
vinf = g*d1*d1*rho_s/(18*mu_f);
mu_g = vinf/rho_s;

phi0 = ((nexp - 2)*phimax-1)/(nexp - 3);
phimax = 0.66;
nexp = 4.7;
phic = 0.2;
sigma0 = 180; 
kexp = 2;
sigedc = sigma0*kexp/phic;
sigeddc = sigma0*kexp*(kexp - 1)/phic^2;
epsilon = 1e-5;
diff_idx = 2;

%%
u = rk122imex(n, T, cfl, L, convec_type, imex_type, test);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%