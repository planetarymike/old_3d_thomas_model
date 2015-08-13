 HOSTNAME=$(shell hostname)
 ifndef HOSTNAME
    $(error no hostname found)
 endif
ifeq ($(HOSTNAME),mike-MacBookAir)
  HOSTNAMERECOGNIZED=true
endif
ifeq ($(HOSTNAME),maven-iuvs-itf.lasp.colorado.edu)
  HOSTNAMERECOGNIZED=true
endif
ifndef HOSTNAMERECOGNIZED
   $(error The machine you are running on, $(HOSTNAME), is not recognized. Please alter the master_fit makefile so that the H coronal fitting code is pointed to the correct model source function location, or disable H coronal fitting)
endif

IDIR=-I./required_procedures/
CC=mpicxx
LIBS=$(shell pkg-config --cflags --libs gsl)
SRCFNSLOCFLG=-D 'SRCFNSLOC="./source_functions/"'

master_fit:
	$(CC) master_fit.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -O3 -o master_fit.x

simulate_coronal_scan:
	$(CC) simulate_coronal_scan.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -O3 -o simulate_coronal_scan.x

generate_S:
	$(CC) generate_S.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -O3 -o generate_S.x
