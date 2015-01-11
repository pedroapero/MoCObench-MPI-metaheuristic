#MoCObench-MPI-metaheuristic
optimisation algorithm for mUBQP-like problems.

##Prerequisites
* mpich
* g++
* gnuplot-x11
* gnuplot

##Compilation
mpic++ MUBQPEval-heuristic-par.cpp

##Execution
mpirun -np [number of processes] ./a.out instance.dat

##Interruption
To interrupt the process and display efficiency informations, send a SIGUSR1 to mpirun's process:
kill -s [mpirun's PID]

##Sequential version
The sequential version takes the MUBQP instance as only parameter:
./MUBQPEval-heuristic-seq.cpp

The ctrl+c interruption will display the efficiency informations.

##TODOs
- [x] Send back solutions in blocks
- [x] Filter whatever M (remove #defines)
- [x] Implement compare function to refactor filter_solutions()
- [x] Avoid copies
- [x] Avoid returning seeds: return index of original solution, and number of flipped bit (requires a vector of sorted non-optimal-but-being-explored solutions)
- [x] Compare with best results on MoCObench website
- [x] Send pools of solutions (seeds) (parametrized size)
- [x] Implement clean sequential version and compare performances
- [-] Evaluate computation time and communication time
- [ ] Grid5000
