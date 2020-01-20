CXXFLAGS=-DNDEBUG -O3 -Wall -Wno-unused-variable -lm -std=c++11
#CXXFLAGS= -g -Wall -Wno-unused-variable
LDFLAGS=
LIBS= -L/usr/local/atlas/lib -llapack -lf77blas -latlas -lgfortran -lm

OBJS_CODE=global.o test_cases.o sed_hsf.o sed_hsf_der.o average.o bc.o fn_flux.o matmult.o matmult_dec.o hornerm.o lminmax_charspeed.o minmax_charspeed.o der_diffusion_tensor.o solid_stress.o solid_stress_der.o solid_stress_der_der.o charspeed.o apply_diffus.o diffusion_tensor.o apply_diffus.o diffus.o diffus_charspeed.o diffusion_matrix.o jacobiana.o jacobiana_dec.o weno5.o convec.o convec_glf.o convec_pvm.o numflux_hll.o muscl.o fhroe.o minmod.o numflux_pvm.o Qcoeff.o residual.o trdsolve.o lsolve.o nwtsolve.o newton_matrix.o rkimex.o update_utilities.o save_solution.o

pvm:  main_imex.o $(OBJS_CODE)
	$(CXX) $(LDFLAGS) -o $@ $^

.PHONY: clean

clean:
		rm   *.o pvm *.h~ *.cc~ Makefile~ *.gch
