#MoCObench-MPI-metaheuristic

optimisation algorithm for mUBQP-like problems

##Prerequisites
* mpich
* g++
* gnuplot-x11
* gnuplot

##Compilation
mpic++ test_MUBQPEval-par.cpp

##Execution
mpirun -np [number of processes] ./a.out instance.dat
