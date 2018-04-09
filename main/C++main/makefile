CC = gcc
CPP = g++
CARGS = -g -fPIC -fexceptions -DIL_STD -lm #DNDEBUG
CPLXINC = -I/home/lapo/alanapanizzi/Documents/cplex/cplex/include
CNCRTINC = -I/home/lapo/alanapanizzi/Documents/cplex/concert/include
CPLXLIB = -L/home/lapo/alanapanizzi/Documents/cplex/cplex/lib/x86-64_linux/static_pic
CNCRTLIB = -L/home/lapo/alanapanizzi/Documents/cplex/concert/lib/x86-64_linux/static_pic
CPLXARGS = -lilocplex -lcplex
CNRTARGS = -lconcert -lpthread

pci_test:
	$(CPP) $(CARGS) $(CPLXINC) $(CNCRTINC) -std=c++11 main.cpp PCI_solver.h -o pci_test $(CNCRTLIB) $(CNRTARGS) $(CPLXLIB) $(CPLXARGS)

clean:
	rm -f *.o *~ pci_test
