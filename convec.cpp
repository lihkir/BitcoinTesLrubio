#include "convec.h"
#include "convec_glf.h"
#include "convec_pvm.h"
#include "test_cases.h"

void convec(Matrix<double>& u, double h, Matrix<double>& K)
{
	struct test_cases* pt_test = get_tests();
	
	if (pt_test->convec_type == 0) {
		convec_pvm(u, h, K);
	} else {
		convec_glf(u, h, K);
	} 
}