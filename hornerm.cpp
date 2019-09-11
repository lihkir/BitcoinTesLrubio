#include "hornerm.h"
#include "globals.h"
#include "matmult.h"
#include "stdio.h"

void hornerm(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &x, std::vector<double> &a, std::vector<std::vector<double>> &p, std::vector<std::vector<double>> &pAx)
{
	global_data* pt_data;

	for (unsigned int i = 0; i < pAx.size(); i++)
		for (unsigned int j = 0; j < pAx[i].size(); j++)
			pAx[i][j] = a[a.size()-1]*x[i][j];

	for (unsigned int k = 0; k < x[0].size(); k++)
    {
		for (unsigned int j = a.size()-1; j-- > 0 ;)
		{
			if (a[j] != 0)
			{
				for (unsigned int i = 0; i < x.size(); i++) 
					p[i][k] = a[j] * x[i][k];				
			}
			else
			{
				for (unsigned int i = 0; i < x.size(); i++ ) 
					p[i][k] = 0;
			}
			matmult(A, pAx, p);
			for (unsigned int i = 0; i < x.size(); i++ ) 
				pAx[i][k] = p[i][k];
		}
	}
}