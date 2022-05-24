# files: source files:  main_chamber.cpp cell_chamber.cpp
#        header files:  cell_chamber.h MersenneTwister.h
# executable file: chamber
#

# first define target file : chamber 
# dependencies are the object files that build the program

chamber: cell_chamber.o main_chamber.o
	g++ -o chamber cell_chamber.o main_chamber.o -fopenmp

# now define how each object file is a target and list dependencies and how
# to build that object file if any dependencies change

cell_chamber.o: cell_chamber.cpp cell_chamber.h MersenneTwister.h
	g++ -c cell_chamber.cpp -fopenmp

main_chamber.o: main_chamber.cpp cell_chamber.h MersenneTwister.h
	g++ -c main_chamber.cpp -fopenmp

clean:
	rm chamber cell_chamber.o main_chamber.o
