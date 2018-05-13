IDIR=-I./required_procedures/ -I/usr/include/eigen3/
CC=clang++
FORT=gfortran
LIBS=-lgsl -lgslcblas -lm -lgfortran -fPIC
MPFLAGS=-fopenmp
#The Gnu scientific library needs to be installed to run this code. Please install GSL.

SRCFNSLOCFLG=-D 'SRCFNSLOC="./source_functions/"'

python_corona_multi_sim:
	$(FORT) -fPIC -g -c -O3 ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	python setup_multi_janus.py build_ext --inplace

python_corona_sim:
	$(FORT) -fPIC -g -c -O3 ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	python setup.py build_ext --inplace

simulate_coronal_scan_pgi:
#	pgfortran -p -g -O0 -fPIC -c ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	pgc++ -acc -Minfo -p -g -std=c++14 -fPIC $(IDIR) -lgsl -lgslcblas -fast -DEIGEN_NO_CPUID simulate_coronal_scan.cpp -o simulate_coronal_scan.x
#	pgc++ -p -g -shared -O0 simulate_coronal_scan.o ./required_procedures/ipbackgroundCFR_fun.o -o simulate_coronal_scan.x

simulate_coronal_scan:
	$(FORT) -Ofast -c ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -std=c++14 -O3 -c simulate_coronal_scan.cpp $(IDIR) $(LIBS) $(MPFLAGS) $(SRCFNSLOCFLG) -o simulate_coronal_scan.o
	$(CC) -std=c++14 -Ofast simulate_coronal_scan.o ./required_procedures/ipbackgroundCFR_fun.o -lgfortran $(IDIR) $(LIBS) $(MPFLAGS) $(SRCFNSLOCFLG) -o simulate_coronal_scan.x

simulate_coronal_scan_debug:
	$(FORT) -fPIC -p -g -c -O0 ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -p -g -c -std=c++14 simulate_coronal_scan.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -O0 -o simulate_coronal_scan.o
	$(CC) -p -g -std=c++14 simulate_coronal_scan.o ./required_procedures/ipbackgroundCFR_fun.o -lgfortran $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -O0 -o simulate_coronal_scan.x

master_fit_all_cpp:
	$(CC) master_fit.cpp -lf2c -u MAIN__ $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -O3 -o master_fit.x

master_fit:
	$(FORT) -O3 -c ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -c master_fit.cpp $(IDIR) $(LIBS) $(MPIFLAGS) $(SRCFNSLOCFLG) -O3 -o master_fit.o
	$(CC) master_fit.o ./required_procedures/ipbackgroundCFR_fun.o -lgfortran $(IDIR) $(LIBS) $(MPIFLAGS) $(SRCFNSLOCFLG) -O3 -o master_fit.x

master_fit_debug:
	$(FORT) -p -g -c -O0 ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -p -g -c master_fit.cpp $(IDIR) $(LIBS) $(MPIFLAGS) $(SRCFNSLOCFLG) -O0 -o master_fit.o
	$(CC) -p -g master_fit.o ./required_procedures/ipbackgroundCFR_fun.o -lgfortran $(IDIR) $(LIBS) $(MPIFLAGS) $(SRCFNSLOCFLG) -O0 -o master_fit.x

multi_corona_test:
	$(FORT) -g -O0 -c ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -g -c test_multi_corona.cpp $(IDIR) $(LIBS) $(MPIFLAGS) $(SRCFNSLOCFLG) -O0 -o test_multi_corona.o
	$(CC) -g test_multi_corona.o ./required_procedures/ipbackgroundCFR_fun.o -lgfortran $(IDIR) $(LIBS) $(MPIFLAGS) $(SRCFNSLOCFLG) -O0 -o test_multi_corona.x

source_function_init:
	$(FORT) -O3 -c ./required_procedures/ipbackgroundCFR_fun.f -o ./required_procedures/ipbackgroundCFR_fun.o
	$(CC) -c source_function_init.cpp $(IDIR) $(LIBS) $(MPIFLAGS) $(SRCFNSLOCFLG) -O3 -o source_function_init.o
	$(CC) source_function_init.o ./required_procedures/ipbackgroundCFR_fun.o -lgfortran $(IDIR) $(LIBS) $(MPIFLAGS) $(SRCFNSLOCFLG) -O3 -o source_function_init.x
