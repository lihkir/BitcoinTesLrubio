#include <math.h>
#include "convec_glf.h"
#include "convec_pvm.h"
#include "convec.h"
#include "test_cases.h"
#include "convec_glf.h"

void convec(std::vector<std::vector<double>>& u, double h, std::vector<std::vector<double>>& K)
{
	struct test_cases* pt_test = get_tests();

	if (pt_test->convec_type == 0)
		convec_pvm(u, h, K);
	else 
		convec_glf(u, h, K);
}