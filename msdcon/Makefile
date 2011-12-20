OBJS=msdcon.o fileio.o
EXECUTABLE=msdcon
CC=gfortran

all: app

debug: FFLAGS = -check all
debug: app

opt: FFLAGS = -fast -O2
opt: app

app: $(OBJS)
	$(CC) $(FFLAGS) $(OBJS) -o $(EXECUTABLE)

msdcon.o: msdcon.f90 fileio.mod
	$(CC) $(FFLAGS) -c msdcon.f90 

fileio.o: fileio.f90
	$(CC) $(FFLAGS) -c fileio.f90 

fileio.mod: fileio.f90
	$(CC) $(FFLAGS) -c fileio.f90

clean:
	rm -f *.o $(EXECUTABLE)
