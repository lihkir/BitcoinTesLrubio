#include <math.h>
#include "Qcoeff.h"
#include "utilities.h"

std::vector<double> Qcoeff(std::vector<double>& S, int idx)
{
	std::vector<double> Qc;

	if (idx == 2) 
	{
		//Q=p(A), p(x)=S_rus
		double S_rus = fabs(S[1]);
		Qc.push_back(S_rus);
	}
	else if (idx == 3) 
	{
		//Q=p(A), p(x)=a x + b, a=(||v2||-||v1||)/(v2-v1), b=(v2*||v1||-v1*||v2||)/(v2-v1)
		if (fabs(S[1] - S[0]) < 1e-8)
			printf("Problemas en HLL %e %e\n", S[0], S[1]);

		double iden = 1.0 / (S[1] - S[0]);

		Qc.push_back((S[1] * fabs(S[0]) - S[0] * fabs(S[1])) * iden);
		Qc.push_back((fabs(S[1]) - fabs(S[0])) * iden);
	}
	else if (idx == 4) 
	{	
		//Q=p(A), p(x)=a ((x)*(x)) + b, a=1.0/2S_0, b= S_0/2
		double S_rus = fabs(S[1]);

		Qc.push_back(S_rus / 2);
		Qc.push_back(0);
		Qc.push_back(1.0 / 2 * S_rus);
	}
	else if (idx == 5) 
	{
		// Q=p(A), p(x)= a + bx + ((cx)*(cx))
		double Sm, SM;
		if (fabs(S[0]) < fabs(S[1])) 
		{
			Sm = S[0];
			SM = S[1];
		}
		else 
		{
			Sm = S[1];
			SM = S[0];
		}
		double iS0 = 1.0 / ((((Sm - SM)) * ((Sm - SM))));
		Qc.push_back(((((SM)) * ((SM))) * Sm * (sgn(Sm) - sgn(SM))) * iS0);
		Qc.push_back((SM * (fabs(SM) - fabs(Sm)) + Sm * (sgn(SM) * Sm - sgn(Sm) * SM)) * iS0);
		Qc.push_back((Sm * (sgn(Sm) - sgn(SM))) * iS0);
	}
	else if (idx == 6) 
	{
		//Q=p(A), p(x)= a + ((cx)*(cx)) + ((ex)*(ex)*(ex)*(ex))
		double SI, SM;
		if (fabs(S[0]) < fabs(S[1])) 
		{
			SI = fabs(S[0]);
			SM = fabs(S[1]);
		}
		else 
		{
			SI = fabs(S[1]);
			SM = fabs(S[0]);
		}

		double Denom = (((SI + SM)) * ((SI + SM)));
		Qc.push_back((SM * SI * (2 * SI + SM)) / (2 * Denom));
		Qc.push_back(0);
		Qc.push_back(1.0 / (2 * SI) + SI / Denom);
		Qc.push_back(0);
		Qc.push_back(-1.0 / (2 * SI * Denom));
	}
	else if (idx == 7) 
	{
		//Q=p(A), p(x)= a + ((cx)*(cx)) + ((ex)*(ex)*(ex)*(ex))
		double S_rus = fabs(S[1]);
		Qc.push_back((3 * S_rus) / 8);
		Qc.push_back(0);
		Qc.push_back(3.0 / (4 * S_rus));
		Qc.push_back(0);
		Qc.push_back(-1.0 / (8 * (((S_rus) * (S_rus) * (S_rus)))));
	}
	else if (idx == 8) 
	{
		//Q=p(A), p(x)= a + ((cx)*(cx)) + ((ex)*(ex)*(ex)*(ex))
		double S0 = max(fabs(S[0]), fabs(S[1]));
		double iS0 = 1.0 / S0;
		Qc.push_back(5.0 / 32 * S0);
		Qc.push_back(0);
		Qc.push_back(54.0 / 32 * iS0);
		Qc.push_back(0);
		Qc.push_back(-27.0 / 32 * ((iS0) * (iS0) * (iS0)));
	}
	else if (idx == 9) 
	{
		double S0 = max(fabs(S[0]), fabs(S[1]));
		double iS0 = 1.0 / S0;
		Qc.push_back(0.178550162085811 * S0);
		Qc.push_back(0);
		Qc.push_back(1.492278834485132 * iS0);
		Qc.push_back(0);
		Qc.push_back(-0.670828996570943 * ((iS0) * (iS0) * (iS0)));
	}
	return Qc;
}