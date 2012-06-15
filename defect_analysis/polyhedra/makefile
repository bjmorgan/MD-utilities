OBJS=main.o cell.o atoms.o polyhedra.o
EXECUTABLE=polyhedra

all: app

debug: FFLAGS = -check all
debug: app

opt: FFLAGS = -fast -O2
opt: app

app: $(OBJS)
	ifort $(FFLAGS) $(OBJS) -o $(EXECUTABLE)

main.o: main.f90 atoms.mod cell.mod polyhedra.mod
	ifort $(FFLAGS) -c main.f90 

cell.o: cell.f90
	ifort $(FFLAGS) -c cell.f90

cell.mod: cell.f90
	ifort $(FFLAGS) -c cell.f90

atoms.o: atoms.f90
	ifort $(FFLAGS) -c atoms.f90

atoms.mod: atoms.f90
	ifort $(FFLAGS) -c atoms.f90

polyhedra.o: polyhedra.f90
	ifort $(FFLAGS) -c polyhedra.f90

polyhedra.mod: polyhedra.f90
	ifort $(FFLAGS) -c polyhedra.f90

clean:
	rm -f *.o $(EXECUTABLE)

