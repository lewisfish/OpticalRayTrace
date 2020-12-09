FCOMP = gfortran
FCFLAGS = -O2 -cpp -freal-4-real-8 -march=native -flto -mavx
FCDEBUG = -g -fbacktrace -fcheck=all -fbounds-check -ffpe-trap=invalid,overflow,underflow,denormal
FCBUILD = -Wall -Wextra -pedantic -std=f2008
FCOPENMP = -fopenmp

PROGRAM =  raytrace

SRCS =      constants.f90 \
			utils.f90 \
			vector_class.f90 \
			stackMod.f90 \
			random_mod.f90 \
			sourceMod.f90 \
			lens.f90 \
			main.f90 
			
OBJECTS = $(SRCS:.f90=.o)

all:    $(PROGRAM)
mp:     FCFLAGS += $(FCOPENMP)
mp:     $(PROGRAM)
debug:  FCFLAGS += $(FCDEBUG)
debug:  $(PROGRAM)
build:  FCFLAGS += $(FCBUILD)
build:  $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ 
	
%.o:  %.f90
	$(FCOMP)  $(FCFLAGS) -c $<

.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD mcgrid
