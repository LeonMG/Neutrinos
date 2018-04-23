#Variables

IDIR = /usr/local/include/aspic
FC = gfortran
FFLAGS = -I$(IDIR)

DEPS = /home/aterno/aspic-0.3.2/include/aspic.h

LDIR = /usr/local/lib
LIBS = -laspic
FLFLAGS = -L$(LDIR)

OBJS = wrapper_NI_aspic.o

# Reglas
main:  $(OBJS) 
	gfortran $(FLFLAGS) $(OBJS) -o xWT $(LIBS)
	@ mkdir -p data
	@ mv xWT data

.PHONY: clean

clean:
	-rm -f *.o data/xWT 

%.o: %.f90 $(DEPS)
	$(FC) -c -o $@ $< $(FFLAGS)



