#include "save_solution.h"
#include <fstream>

void save_solution(Matrix<double> &u0, int N_rows, int N_cols, int imex_type, int convec_type, double T)
{
    std::string init_name = "out";
    std::string strN_rows = std::to_string(N_rows);
    std::string strN_cols = std::to_string(N_cols);
    std::string strImex_type = std::to_string(imex_type);
	std::string strConvec_type = std::to_string(convec_type);
    std::string strTime = std::to_string((int)T);

    std::string nameh = init_name + "_m" + strN_rows + "_n" +  strN_cols + "_imex" +  strImex_type + "_convec" +  strConvec_type + "_time" +  strTime; 
    std::ofstream fileout(nameh);
  	for (int i = 1; i <= N_rows; i++)
    	for (int j = 1; j <= N_cols; j++) 
      		fileout << std::fixed << std::setprecision(16) << std::scientific << u0(i, j) << std::endl;
  	fileout.close();
}