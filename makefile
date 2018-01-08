IDIR=-I./required_procedures/ -I$(CURC_GSL_ROOT)/include/ -L$(CURC_GSL_ROOT)/lib/ -L$(CURC_INTEL_ROOT)/intel64/ -I$(CURC_INTEL_ROOT)/include/
CC=icpc
FORT=ifort
CFORTLINK=ifort
CLIBS= -lgsl -lgslcblas -lm -lfor -fPIC
LLIBS=-lstdc++ -fPIC -nofor-main -lgsl -lgslcblas -lm
#LIBS=-lgsl -lgslcblas -lm -lgfortran -fPIC
#The Gnu scientific library needs to be installed to run this code. Please install GSL.

python_corona_multi_sim:
	CC=icpc python setup_multi_summit.py build_ext --inplace

python_corona_sim:
	$(FORT) -fPIC -g -c -O3 ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	python setup.py build_ext --inplace

simulate_coronal_scan:
	$(CC) -c simulate_coronal_scan.cpp $(IDIR) $(CLIBS) -O3 -o simulate_coronal_scan.o
	$(FORT) ./required_procedures/ipbackgroundCFR_fun.f simulate_coronal_scan.o $(IDIR) $(LLIBS) -O3 -o simulate_coronal_scan.x

simulate_coronal_scan_debug:
	$(CC) -pg -c simulate_coronal_scan.cpp $(IDIR) $(CLIBS) -O0 -o simulate_coronal_scan.o
	$(FORT) -pg ./required_procedures/ipbackgroundCFR_fun.f simulate_coronal_scan.o $(IDIR) $(LLIBS) -O0 -o simulate_coronal_scan.x

master_fit_all_cpp:
	$(CC) master_fit.cpp -lf2c -u MAIN__ $(IDIR) $(CLIBS) -O3 -o master_fit.x

master_fit:
	$(CC) -c master_fit.cpp $(IDIR) $(CLIBS) -O3 -o master_fit.o
	$(FORT) ./required_procedures/ipbackgroundCFR_fun.f master_fit.o $(LLIBS) $(IDIR) -O3 -o master_fit.x

master_fit_debug:
	$(CC) -pg -c master_fit.cpp $(IDIR) $(CLIBS) -O0 -o master_fit.o
	$(FORT) -pg -c ./required_procedures/ipbackgroundCFR_fun.f master_fit.o $(LLIBS) $(IDIR) -O0 -o master_fit.x

multi_corona_test:
	$(FORT) -g -O0 -c ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -g -c test_multi_corona.cpp $(IDIR) $(CLIBS) $(MPIFLAGS) -O0 -o test_multi_corona.o
	$(CFORTLINK) -g test_multi_corona.o ./required_procedures/ipbackgroundCFR_fun.o $(IDIR) $(FLIBS) $(MPIFLAGS) -O0 -o test_multi_corona.x

source_function_init:
	$(CC) -c source_function_init.cpp $(IDIR) $(CLIBS) -O3 -o source_function_init.o
	$(FORT) ./required_procedures/ipbackgroundCFR_fun.f source_function_init.o $(LLIBS) $(IDIR) -O3 -o source_function_init.x
