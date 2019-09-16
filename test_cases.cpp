#include "test_cases.h"
#include <pthread.h>
#include <cstdio>

struct test_cases *get_tests()
{
    struct test_cases *pt_test = new test_cases;    

	extern pthread_mutex_t lock;
	extern int idx_test;

	pthread_mutex_lock(&lock);
	if (idx_test == 1)
	{
        pt_test->beta.push_back(60.0);
		pt_test->beta.push_back(67.5);
		pt_test->beta.push_back(75.0);
		pt_test->beta.push_back(82.5);
		pt_test->beta.push_back(90.0);
		pt_test->beta.push_back(97.5);
		pt_test->beta.push_back(105.0);
		pt_test->beta.push_back(112.5);
		pt_test->beta.push_back(120.0);
		pt_test->delta.push_back(4.057404490000000e-08);
  		pt_test->delta.push_back(1.965604000000000e-08);
  		pt_test->delta.push_back(9.950262001000001e-09);
  		pt_test->delta.push_back(4.759068196000001e-09);
  		pt_test->delta.push_back(2.341688881000000e-09);
  		pt_test->delta.push_back(1.168614225000000e-09);
  		pt_test->delta.push_back(5.669161000000000e-10);
  		pt_test->delta.push_back(3.722342121000000e-11);
		pt_test->M_rows = pt_test->beta.size();
		pt_test->phimax = 50;
		pt_test->int_form = 5;
		pt_test->gc =2;
		pt_test->gc_id = 1;
		pt_test->g = 9.81;
		pt_test->mu_c = 6.5e-3;
		pt_test->rho_d = 1090;
		pt_test->rho_c = 880;
		pt_test->vinf = pt_test->g*(pt_test->rho_d-pt_test->rho_c)/(18*pt_test->mu_c);
		pt_test->L = 0.02;
		pt_test->nexp = 4.65;
	}
	pthread_mutex_lock(&lock);
	return pt_test;
}