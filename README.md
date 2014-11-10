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

##TODOs
- [ ] Send back solutions in blocks
- [ ] Filter whatever M (remove #defines)
- [ ] Implement compare function to refactor filter_solutions()
- [ ] Avoid copies
- [ ] Avoid returning seed but filter with it
- [ ] Avoid returning and neighbour solutions: return index of original solution, and number of flipped bit (requires a vector of sorted non-optimal-but-being-explored solutions)
- [ ] Compare with best results on MoCObench website
- [ ] Implement clean sequential version and compare performances
- [ ] Send pools of solutions (seeds) (parametrized size)
