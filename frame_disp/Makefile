OBJS=frame_disp.o fileio.o
EXECUTABLE=frame_disp

all: app

app: $(OBJS)
	ifort $(OBJS) -o $(EXECUTABLE)

frame_disp.o: frame_disp.f90 fileio.mod
	ifort -c frame_disp.f90 

fileio.o: fileio.f90
	ifort -c fileio.f90 

fileio.mod: fileio.f90
	ifort -c fileio.f90

clean:
	rm -f *.o $(EXECUTABLE)
