#include "bc.h"
#include "globals.h"

void bc(std::vector<std::vector<double>> &u)
{
	global_data* pt_data;
	int N = u[0].size() - 2 * pt_data->gc;
	
	for (unsigned int i = 0; i < u.size(); i++)
	{
		for (auto j = 0; j != pt_data->gc; j++)
		{
			if (pt_data->gc_id == 1)
			{
				u[i][j] = u[i][N + j];
				u[i][pt_data->gc + N + j] = u[i][pt_data->gc + j];
			}
			else 
			{
				u[i][j] = 1e10;
				u[i][pt_data->gc + N + j] = 1e10;
			}
		}
	}
}