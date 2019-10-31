LIBS	= -I/usr/local/include/ -L/usr/lib64 -lgsl -lgslcblas -lboost_filesystem -lboost_system -lm -std=c++0x -O4
CCOMP 	= mpic++
LINKER 	= -o 
CFLAGS 	= -c

program     = evolve_exec
objects     = main.o thimble_lattice.o

$(program): $(objects)
	mpic++ $(objects) -o $(program) $(LIBS)
	
main.o: main.cpp 
	mpic++ -c main.cpp $(LIBS)

thimble_lattice.o: thimble_lattice.cpp thimble_lattice.h matrix.hpp
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


