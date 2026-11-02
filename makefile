hamiltonian: graph.o setup.o solver.o main.o
	g++ -g -o hamiltonian graph.o setup.o solver.o main.o

graph.o: graph.cpp
	g++ -g -c graph.cpp

setup.o: setup.cpp
	g++ -g -c setup.cpp

solver.o: solver.cpp
	g++ -g -c solver.cpp

main.o: main.cpp
	g++ -g -c main.cpp