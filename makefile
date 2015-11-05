IDIR=-I./required_procedures/
CC=mpicxx
FORT=gfortran
LIBS=$(shell pkg-config --cflags --libs gsl)
ifndef LIBS
	$(error The Gnu scientific library needs to be installed to run this code. Please install GSL.)
endif

SRCFNSLOCFLG=-D 'SRCFNSLOC="./source_functions/"'
TABDATALOCFLG=-D 'SRCFNSLOC="./source_functions/"'

simulate_coronal_scan:
	$(FORT) -O3 -c ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -c simulate_coronal_scan.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) $(TABDATALOCFLAG) -O3 -o simulate_coronal_scan.o
	$(CC) simulate_coronal_scan.o ./required_procedures/ipbackgroundCFR_fun.o -lgfortran $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -O3 -o simulate_coronal_scan.x

simulate_coronal_scan_debug:
	$(FORT) -g -c -O0 ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -g -c simulate_coronal_scan.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) $(TABDATALOCFLAG) -O0 -o simulate_coronal_scan.o
	$(CC) -g simulate_coronal_scan.o ./required_procedures/ipbackgroundCFR_fun.o -lgfortran $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -O0 -o simulate_coronal_scan.x

master_fit:
	$(FORT) -O3 -c ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -c master_fit.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) $(TABDATALOCFLAG) -O3 -o master_fit.o
	$(CC) master_fit.o ./required_procedures/ipbackgroundCFR_fun.o -lgfortran $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -O3 -o master_fit.x

master_fit_debug:
	$(FORT) -g -c -O0 ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -g -c master_fit.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) $(TABDATALOCFLAG) -O0 -o master_fit.o
	$(CC) -g master_fit.o ./required_procedures/ipbackgroundCFR_fun.o -lgfortran $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -O0 -o master_fit.x

