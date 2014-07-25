FC 	:= gfortran
FFLAGS	:= -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow -g -O3 -fopenmp
OBJS	:= m_particle_swarm.o m_io.o

INCDIRS	:= ../../fosito ../../particle_core
LIBDIRS := ../../fosito ../../particle_core
LIBS	:= particle_core fosito

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

../%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: ../particle_swarm

clean:
	$(RM) -f *.o *.mod particle_swarm

# Dependency information
m_particle_swarm.o:	m_io.o
../particle_swarm:		$(OBJS)
