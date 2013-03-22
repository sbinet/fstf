CXXFLAGS=-O2 -fopenmp -g -mavx -march=native -Iinclude/eigen3 -I/opt/intel/composer_xe_2013.1.117/mkl/include/ -Iinclude/root
LDFLAGS = -L/opt/intel/composer_xe_2013.1.117/compiler/lib/intel64/ -lsvml -lgomp -L/opt/intel/composerxe-2013/mkl/lib/intel64/ -lmkl_core -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_avx

all: main.o testAVX.o testAVX.h
	g++ $(LDFLAGS) testAVX.o main.o -lgomp -lCLHEP -o kernel

testAVX.o: testAVX.cpp testAVX.h
main.o: main.cpp testAVX.h

clean:
	-rm -rf *.o kernel
