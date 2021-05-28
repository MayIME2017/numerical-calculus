FC = gfortran
FFLAGS = -march=native -mtune=native -O3
FDEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow,denormal
IDIR = -I/usr/local/include/fgsl/
LIBS = -lgsl -lfgsl #-lopenblas -lm

TARGET1 = trabalho4_fortran
TARGET2 = trabalho4_value

OBJECTS1 = system_parameters.o inter_approx.o main.o
OBJECTS2 = system_parameters.o inter_approx.o main_value.o

all: TARGET1 TARGET2

TARGET1: $(OBJECTS1)
	$(FC) $(FFLAGS) -o $(TARGET1) $(IDIR) $^ $(LIBS)

TARGET2: $(OBJECTS2)
	$(FC) $(FFLAGS) -o $(TARGET2) $(IDIR) $^ $(LIBS)

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $(IDIR) $<
