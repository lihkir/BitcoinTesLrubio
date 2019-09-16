CXXFLAGS=-DNDEBUG -O3 -Wall -Wno-unused-variable -lm -std=c++11
#CXXFLAGS= -g -Wall -Wno-unused-variable
LDFLAGS=
LIBS= -L/usr/local/atlas/lib -llapack -lf77blas -latlas -lgfortran -lm

OBJS_CODE=test_cases.o average.o bc.o fn_flux.o matmult.o hornerm.o extern.o lminmax_charspeed.o 

pvm:  main_pvm_sed.o $(OBJS_CODE)
	$(CXX) $(LDFLAGS) -o $@ $^

.PHONY: clean

clean:
		rm   *.o pvm *.h~ *.cc~ Makefile~ *.gch
