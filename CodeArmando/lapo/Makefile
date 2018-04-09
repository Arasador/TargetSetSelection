CC = gcc
CPP = g++
CARGS = -O3 -fPIC -fexceptions -DNDEBUG -DIL_STD  -lm
CPLXINC = -I/home/lapo/cplex/cplex125/cplex/include 
CNCRTINC = -I/home/lapo/cplex/cplex125/concert/include 
CPLXLIB = -L/home/lapo/cplex/cplex125/cplex/lib/x86-64_sles10_4.1/static_pic 
CNCRTLIB = -L/home/lapo/cplex/cplex125/concert/lib/x86-64_sles10_4.1/static_pic 
CPLXARGS = -lilocplex -lcplex
CNRTARGS = -lconcert -lpthread

pdtspms:
	$(CPP) $(CARGS) $(CPLXINC) $(CNCRTINC) ../main.cpp ../TSPPDMS.cpp ../Dinic.cpp ../TSPPDMS.h ../MaxFlow.h ../Dinic.h -o pdtspms $(CNCRTLIB) $(CNRTARGS) $(CPLXLIB) $(CPLXARGS)

clean:
	rm -f *.o *~ pdtspms
