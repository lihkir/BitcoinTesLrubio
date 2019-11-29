#include <stdio.h>
#include "hornerm.h"
#include "matmult.h"
#include "utilities.h"

void hornerm(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& dif, std::vector<double>& Qc, std::vector<std::vector<double>>& fh, std::vector<std::vector<double>>& Qurmul)
{
	for (unsigned int j = 0; j < dif[0].size(); j++)
		for (unsigned int i = 0; i < dif.size(); i++)
			Qurmul[i][j] = Qc[Qc.size() - 1] * dif[i][j];
		
	for (unsigned int k = 0; k < dif[0].size(); k++)
	{
		for (int j = Qc.size() - 2; j >= (0); j-=1)
		{
			if (Qc[j] != 0)
			{
				for (unsigned int i = 0; i < dif.size(); i++)
					fh[i][k] = Qc[j] * dif[i][k];
			}
			else
			{
				for (unsigned int i = 0; i < dif.size(); i++)
					fh[i][k] = 0;
			}
			matmult(A, Qurmul, fh);
			for (unsigned int i = 0; i < dif.size(); i++)
				Qurmul[i][k] = fh[i][k];
		}
	}
}