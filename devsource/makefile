LIBS	= -I/usr/local/include/ -L/usr/lib64 -lboost_filesystem -lboost_system -lquadmath -lm -std=c++0x -g
CCOMP 	= mpic++
LINKER 	= -o 
CFLAGS 	= -c

program     = evolve_exec
objects     = main.o thimble_lattice.o scalars_and_interactions.o

$(program): $(objects)
	mpic++ $(objects) -o $(program) $(LIBS)
	
main.o: main.cpp 
	mpic++ -c main.cpp $(LIBS)

thimble_lattice.o: thimble_lattice.cpp thimble_lattice.hpp matrix.hpp
	mpic++ -c thimble_lattice.cpp $(LIBS)

scalars_and_interactions.o: scalars_and_interactions.cpp thimble_lattice.hpp
	mpic++ -c scalars_and_interactions.cpp $(LIBS)

clean:
	rm -f  core $(objects) $(program)

#the make file is of the form
#   target: source
#   #   command
#
##this means that target depends on source, and if source has changed
#then command is implemented, creating a new target. 
#command must be preceeded by a tab.


