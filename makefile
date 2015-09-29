IDIR=-I./required_procedures/
CC=mpicxx
LIBS=$(shell pkg-config --cflags --libs gsl)
ifndef LIBS
	$(error The Gnu scientific library needs to be installed to run this code. Please install GSL)
endif

SRCFNSLOCFLG=-D 'SRCFNSLOC="./source_functions/"'
TABDATALOCFLG=-D 'SRCFNSLOC="./source_functions/"'

simulate_coronal_scan:
	$(CC) simulate_coronal_scan.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) $(TABDATALOCFLAG)-O3 -o simulate_coronal_scan.x

generate_S:
	$(CC) generate_S.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) $(TABDATALOCFLAG) -O3 -o generate_S.x
