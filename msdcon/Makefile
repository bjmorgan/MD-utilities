OBJS=msdcon.o fileio.o
EXECUTABLE=msdcon

all: app

debug: FFLAGS = -check all
debug: app

opt: FFLAGS = -fast -O2
opt: app

app: $(OBJS)
	ifort $(FFLAGS) $(OBJS) -o $(EXECUTABLE)

msdcon.o: msdcon.f90 fileio.mod
	ifort $(FFLAGS) -c msdcon.f90 

fileio.o: fileio.f90
	ifort $(FFLAGS) -c fileio.f90 

fileio.mod: fileio.f90
	ifort $(FFLAGS) -c fileio.f90

clean:
	rm -f *.o $(EXECUTABLE)
