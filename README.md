#MoCObench-MPI-metaheuristic
Optimisation algorithm for MUBQP-like problems.

##Prerequisites
* mpich
* g++
* gnuplot-x11
* gnuplot

##Compilation
```
mpic++ MUBQPEval-heuristic-par.cpp
```

##Execution
```
mpirun -np [number of processes] ./a.out instance.dat [solutions pool size]
```

The larger the number of solutions to send to the workers, the later the filtering. The smaller the
number of solutions, the more important the losses in communication time.

##Interruption
To interrupt the process and display efficiency informations, send a SIGUSR1 to mpirun's process:```kill -s [mpirun's PID]```

##Sequential version
The sequential version takes the MUBQP instance as only parameter:
```
./MUBQPEval-heuristic-seq instance.dat
```

The ctrl+c interruption will display the efficiency informations.

##Deployment
* Local deployment:
```
./launch_instances.sh [number of instances] [instance duration (seconds)] [number of processes]
[instance file] [pool size]
```

* Passive Grid5000 deployment:
```
oarsub -l /nodes=[number of nodes],walltime=xx:xx:xx "./launch_passive_instances.sh [number of
instances] [instance duration (seconds)] [instance file] [pool size]"
```

Objective vectors will be stored in files named: [instance file]_[duration]_[number of
processes]_[pool size]_[instance number].txt

##TODOs
- [x] Send back solutions in blocks
- [x] Filter whatever M (remove #defines)
- [x] Implement compare function to refactor filter_solutions()
- [x] Avoid copies
- [x] Avoid returning seeds: return index of original solution, and number of flipped bit (requires
	a vector of sorted non-optimal-but-being-explored solutions)
- [x] Compare with best results on MoCObench website
- [x] Send pools of solutions (seeds) (parametrized size)
- [x] Implement clean sequential version and compare performances
- [ ] Evaluate computation time and communication time
- [ ] Dynamic pool size: 25% ; 50% ; 75% or 100% of (number of solutions to explore) / (number of
	processes)
- [x] Grid5000
