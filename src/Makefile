SRCS	:= m_particle_swarm.f90 m_io.f90 m_config.f90
OBJS	:= $(SRCS:%.f90=%.o)
INCDIRS	:= ../particle_core
LIBDIRS := ../particle_core
LIBS	:= particle_core
FC	:= gfortran
FFLAGS	:= -Wall -O2 -std=f2008 -fopenmp -Wno-unused-dummy-argument

ifeq ($(DEBUG), 1)
	FFLAGS += -fcheck=all -g -pg -ffpe-trap=invalid,zero,overflow \
	-pedantic -finit-real=nan
endif

ifeq ($(PROF), gprof)
	FFLAGS += -pg
else ifeq ($(PROF), gperftools)
	LIBS += profiler
endif

vpath %.f90 config_fortran

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

../%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: 	../particle_swarm

clean:
	$(RM) -f *.o *.mod ../particle_swarm

# Dependency information
LIBPC := ../particle_core/libparticle_core.a
../particle_swarm: $(OBJS)
$(OBJS): $(LIBPC)
m_particle_swarm.o: m_config.o
particle_swarm.o: m_particle_swarm.o m_io.o m_config.o
