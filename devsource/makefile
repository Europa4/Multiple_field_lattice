LIBS	= -I/usr/local/include/ -L/usr/lib64 -lgsl -lgslcblas -lboost_filesystem -lboost_system -lm -std=c++0x -O4
CCOMP 	= mpic++
LINKER 	= -o 
CFLAGS 	= -c

program     = evolve_exec
objects     = main.o initialise.o c_phi.o matrix_multiplication.o thimble_lattice.o

$(program): $(objects)
	mpic++ $(objects) -o $(program) $(LIBS)
	
main.o: main.cpp Prot.h
	mpic++ -c main.cpp $(LIBS)
	
initialise.o: initialise.cpp Prot.h
	mpic++ -c initialise.cpp $(LIBS)
	
c_phi.o: c_phi.cpp Prot.h
	mpic++ -c c_phi.cpp $(LIBS)
	
matrix_multiplication.o: matrix_multiplication.cpp Prot.h
	mpic++ -c matrix_multiplication.cpp $(LIBS)

thimble_lattice.o: thimble_lattice.cpp thimble_lattice.h
	mpic++ -c thimble_lattice.cpp $(LIBS)

clean:
	rm -f  core $(objects) $(program)

#the make file is of the form
#   target: source
#   #   command
#
##this means that target depends on source, and if source has changed
#then command is implemented, creating a new target. 
#command must be preceeded by a tab.


