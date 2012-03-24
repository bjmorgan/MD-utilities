OBJS=defect.o fileio.o sites.o
EXECUTABLE=defect_new

all: app

debug: FFLAGS = -check all
debug: app

opt: FFLAGS = -fast -O2
opt: app

app: $(OBJS)
	ifort $(FFLAGS) $(OBJS) -o $(EXECUTABLE)

defect.o: defect.f90 class_file.mod class_tetrahedron.mod class_octahedron.mod
	ifort $(FFLAGS) -c defect.f90 

fileio.o: fileio.f90
	ifort $(FFLAGS) -c fileio.f90 

class_file.mod: fileio.f90
	ifort $(FFLAGS) -c fileio.f90

sites.o: sites.f90
	ifort $(FFLAGS) -c sites.f90

class_tetrahedron.mod: sites.f90
	ifort $(FFLAGS) -c sites.f90

class_octahedron.mod: sites.f90
	ifort $(FFLAGS) -c sites.f90

clean:
	rm -f *.o $(EXECUTABLE)
