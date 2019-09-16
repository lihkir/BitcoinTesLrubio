#include <stdio.h>
#include <math.h>
#include "sed_hsf.h"
#include "sed_hsf_der.h"

void lminmax_charspeed(std::vector<double> &phil, std::vector<double> &phir, std::vector<double> &S)
{
    N=length(delta);

s_phil=0;                                                                                                                                   
etal=0;
for (j=1; j<=N; j++ ) {
s_phil=s_phil+phil(j);
etal=etal+delta(j)*phil(j);
}

s_phir=0;
etar=0;
for (j=1; j<=N; j++ ) {
s_phir=s_phir+phir(j);
etar=etar+delta(j)*phir(j);
}

s_delta=s_phil-s_phir;
eta_delta=etal-etar;

w1=sed_hsf(s_phil);
w2=sed_hsf(s_phir);

maxv1=delta(1)*w1;
maxv2=delta(1)*w2;
max_0=MAX(maxv1, maxv2);

// if  ( idx_q == 6 && N > 1 ) {
//     smaxv1=delta(2)*w1;
//     smaxv2=delta(2)*w2;
//     S(3)=MAX(smaxv1, smaxv2); //Segundo Maxima Cota de Autovalores
// }

i_buf=1;

grid * pbuf_opt=zeros(3,1); grid & buf_opt=*pbuf_opt;
buf_opt(i_buf)=w1*delta(N)+etal*sed_hsf_der(s_phil); i_buf=i_buf+1;
buf_opt(i_buf)=w2*delta(N)+etar*sed_hsf_der(s_phir); i_buf=i_buf+1;

if  ( eta_delta !=0 && s_delta != 0 ) {
Q=nexp*eta_delta*s_delta+((s_delta)*(s_delta))*delta(N);
alpha=((1-s_phir)*(delta(N)*s_delta+eta_delta)-(nexp-1)*etar*s_delta)/Q;
    if  ( alpha > 0 && alpha < 1) {
vphi=s_phir+alpha*s_delta;
etav=etar+alpha*eta_delta;
w=sed_hsf(vphi);
buf_opt(i_buf)=w*delta(N)+etav*sed_hsf_der(vphi); i_buf=i_buf+1;
    }
}

min_n=find_extrema(buf_opt, i_buf-1);

S(1)=min_n;
S(2)=max_0;
//printf("lminmax_charspeed: %.12e %.12e %.12e %.12e || % .12e % .12e\n", s_phil, s_phir, etal, etar, S(1), S(2));
 fflush(stdout);
 delete(pbuf_opt);

}

double find_extrema(grid& v, int l )
{
  double vmin;
  int i;
  
  
  //vector: v
//int: i, l
//double: vmin

vmin=v(1);

for (i=2; i<=l; i++ ) {
vmin=MIN(vmin, v(i));
}
 return vmin;
}
