CC = gcc
CPP = g++
CARGS = -O3 -fPIC -fexceptions -DIL_STD -DNDEBUG -lm #-g
CPLXINC = -I/opt/ibm/ILOG/CPLEX_Studio1271/cplex/include
CNCRTINC = -I/opt/ibm/ILOG/CPLEX_Studio1271/concert/include
CPLXLIB = -L/opt/ibm/ILOG/CPLEX_Studio1271/cplex/lib/x86-64_linux/static_pic
CNCRTLIB = -L/opt/ibm/ILOG/CPLEX_Studio1271/concert/lib/x86-64_linux/static_pic
CPLXARGS = -lilocplex -lcplex
CNRTARGS = -lconcert -lpthread

pci_test:
	$(CPP) $(CARGS) $(CPLXINC) $(CNCRTINC) -std=c++11 main.cpp PCI_solver.h -o pci_test $(CNCRTLIB) $(CNRTARGS) $(CPLXLIB) $(CPLXARGS)

clean:
	rm -f *.o *~ pci_test
