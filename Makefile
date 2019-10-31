all: tools
	g++ --std=c++11 --openmp Equation.o DoubleMass.o RungeKutta.o BifurcationTree.o main.cpp -o main

eqs:
	g++ --std=c++11 Equation.cpp -c -o Equation.o
	g++ --std=c++11 DoubleMass.cpp -c -o DoubleMass.o

method: eqs
	g++ --std=c++11 RungeKutta.cpp -c -o RungeKutta.o

tools: method
	g++ --std=c++11 BifurcationTree.cpp -c -o BifurcationTree.o

clean:
	rm -rf *.o main