CC = gcc
CPP = g++
CARGS = -O3 -fPIC -fexceptions -DIL_STD -DNDEBUG -lm
CPLXINC = -I/opt/ibm/ILOG/CPLEX_Studio1271/cplex/include
CNCRTINC = -I/opt/ibm/ILOG/CPLEX_Studio1271/concert/include
CPLXLIB = -L/opt/ibm/ILOG/CPLEX_Studio1271/cplex/lib/x86-64_linux/static_pic
CNCRTLIB = -L/opt/ibm/ILOG/CPLEX_Studio1271/concert/lib/x86-64_linux/static_pic
CPLXARGS = -lilocplex -lcplex
CNRTARGS = -lconcert -lpthread

#objects = preprocessing.o s_cutter.o separation.o pci_solver.o

pci_test: preprocessing.o pci_solver.o
	$(CPP) $(CARGS) $(CPLXINC) $(CNCRTINC) \
	-std=c++11 -o pci_test preprocessing.o pci_solver.o  s_cutter.o separation.o main.cpp \
	$(CNCRTLIB) $(CNRTARGS) $(CPLXLIB) $(CPLXARGS)

pci_solver.o: s_cutter.o separation.o includes.h
	$(CPP) $(CARGS) $(CPLXINC) $(CNCRTINC) \
	-c -o pci_solver.o PCI_solver.cpp \
	$(CNCRTLIB) $(CNRTARGS) $(CPLXLIB) $(CPLXARGS)

preprocessing.o: preprocessing.cpp preprocessing.h
	$(CPP) -c -o preprocessing.o preprocessing.cpp

s_cutter.o: S_cutter.cpp S_cutter.h
	$(CPP) -c -o s_cutter.o S_cutter.cpp

separation.o: Separation.cpp Separation.h ClusterHCS.cpp
	$(CPP) -c -o separation.o Separation.cpp

clean:
	rm -f *.o *~ pci_test
