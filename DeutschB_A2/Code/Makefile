
OBJS1 = numtype.o BUharmo.o 

PROG1 = harmo


F95 = gfortran

F95FLAGS = -O3 -funroll-loops  -fexternal-blas

LIBS = -framework Accelerate

LDFLAGS = $(LIBS)

all: $(PROG1) 

$(PROG1): $(OBJS1)
	$(F95) $(LDFLAGS) -o $@ $(OBJS1) 

clean:
	rm -f $(PROG1) *.{o,mod}

.SUFFIXES: $(SUFFIXES) .f95

.f95.o:
	$(F95) $(F95FLAGS) -c $<

