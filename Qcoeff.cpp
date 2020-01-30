#include <math.h>
#include "Qcoeff.h"

void Qcoeff(Matrix<double> &S, int idx, Matrix<double> &Qc)
{
  double Sm, SM, S_rus, iden, iS0, SI, Denom1, Denom2, S0;
  
  if (idx == 2)
    {/** Q = p(A), p(x)=S_rus **/
      S_rus = max(abs(S(1)), abs(S(2)));
      Qc(1) = S_rus;
    }
  else if (idx == 3)
    {/** Q = p(A), p(x)=a x + b, a=(|v2|-|v1|)/(v2-v1), b=(v2*|v1|-v1*|v2|)/(v2-v1) **/
      if (abs(S(2) - S(1)) < 1e-8)
        printf("Problemas en HLL %e %e\n", S(1), S(2));
      
      iden = 1/(S(2)-S(1));
      Qc(2) = (abs(S(2))-abs(S(1)))*iden;
      Qc(1) = (S(2)*abs(S(1))-S(1)*abs(S(2)))*iden; 
    }
  else if (idx == 4)
    {/** Q = p(A), p(x)=a x^2 + b, a = 1/2S_0, b = S_0/2 **/
      S_rus = max(abs(S(1)), abs(S(2)));
      Qc(3) = 1/(2*S_rus);
      Qc(2) = 0;
      Qc(1) = S_rus/2;
    }
  else if (idx == 5)
    {/** Q = p(A), p(x) = a + bx + cx^2 **/
      if (abs(S(1)) < abs(S(2)))
	{
	  Sm=S(1);
	  SM=S(2);
	}
      else
	{
	  Sm=S(2);
	  SM=S(1);
	}
      iS0 = 1/pow(Sm - SM, 2);
      Qc(1)=(pow(SM, 2)*Sm*(sgn(Sm)-sgn(SM)))*iS0;
      Qc(2)=(SM*(abs(SM)-abs(Sm))+Sm*(sgn(SM)*Sm-sgn(Sm)*SM))*iS0;
      Qc(3)=(Sm*(sgn(Sm)-sgn(SM)))*iS0;
    }
  else if (idx == 6)
    {/** Q = p(A), p(x) = a + cx^2 + ex^4 **/
      if (abs(S(1)) < abs(S(2)))
	{
	  SI = max(abs(S(1)), abs(S(2)));
	  SM = max(S(1),S(2));
	}   
      else
	{
	  SI = max(abs(S(1)),abs(S(2))); 
	  SM = max(S(1), S(2));  
	}
      Denom1 = 1/pow(abs(SI)+abs(SM), 2);
      Denom2 = 1/abs(SM);
      Qc(1) = 0.5*abs(SM)*abs(SI)*(2*abs(SM)+abs(SI))*Denom1;
      Qc(2) = 0;
      Qc(3) = 0.5*Denom2+abs(SM)*Denom1;
      Qc(4) = 0;
      Qc(5) = -0.5*Denom1*Denom2;
    }
  else if (idx == 7)
    {/** Q = p(A), p(x) = a + cx^2 + ex^4 **/
      S_rus = max(abs(S(1)), abs(S(2)));
      Qc(1) = (3*S_rus)/8;
      Qc(2) = 0;
      Qc(3) = 3/(4*S_rus);
      Qc(4) = 0;
      Qc(5) = -1/(8*(pow(S_rus,3)));
    }
  else if (idx == 8)
    {/** Q = p(A), p(x) = a + cx^2 + ex^4 **/
      S0 = max(abs(S(1)), abs(S(2)));
      iS0=1/S0;
      Qc(1) = 5/32*S0;
      Qc(2) = 0;
      Qc(3) = 54/32*iS0;
      Qc(4) = 0;
      Qc(5) = -27/32*pow(iS0, 3);
    }   
  else if (idx == 9)
    {
      S0 = max(abs(S(1)), abs(S(2)));
      iS0=1/S0;
      Qc(1) = 0.178550162085811*S0;
      Qc(2) = 0;
      Qc(3) = 1.492278834485132*iS0;
      Qc(4) = 0;
      Qc(5) =-0.670828996570943*pow(iS0, 3);
    }
}
