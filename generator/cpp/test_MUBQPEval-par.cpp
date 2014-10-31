#include <iostream>
#include <vector>
#include <string>
#include <climits>

#include "mubqpEval.h"
#include "mubqpEval-C-version.h"
#include "mubqpEval-C-version.c"
#include <unistd.h>
#include "mpi.h"
#include "signal.h"
#include "gnuplot_i.hpp"

using namespace std;

/**************************************************************
 * MPI tools
 **************************************************************/
int self;			/* my rank among processes */
int procs;			/* number of processes */
int rnbr, lnbr;			/* right neighbout and left neighbour */
int p; /* process number */
MPI_Comm com;
MPI_Status status;

#define PROC_NULL 0

void interruption_signal_handler(int sig) {
	fflush(stdout);
	MPI_Abort(com,EXIT_SUCCESS);
	MPI_Finalize();
	exit(sig);
}

/***************************************************************/
/*                 Main                                        */
/***************************************************************/

// TODO: use real instance values instead of this dirty workaround
#define SOLUTION_LENGTH 8 // N
#define RESULT_DIMENSION 2 // M

// the tuple solution - vector
struct result_t {
	// TODO: remove hard coded values!! (input is size N, below output is M, but impossible to use them to fix array size)
	unsigned int input[SOLUTION_LENGTH]; // tested matrix 
	int output[RESULT_DIMENSION]; // objvect
};

int main(int argc, char *argv[]) {
  //-----------------------------------------------
	// MPI initialisation
	com = MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(com, &procs);
	MPI_Comm_rank(com, &self);


  if (argc != 2) {
    std::cerr << "Invalid number of parameters. Use the command line: ./test_MUBQPEval instance.dat" << std::endl;
    exit(1);
  }

	if (procs < 2) {
		printf("Two processes minimum, it's master/slaves schema.");
		exit(1);
	}

	/* we will use the C version instead
  MUBQPEval mubqp(argv[1]);
  // size of the bit string
  const unsigned int N = mubqp.getN();
  // objective space dimension
  const unsigned int M = mubqp.getM();
	*/

	Instance instance;
	loadInstance(argv[1], &instance);
  unsigned int N = instance.N;
  unsigned int M = instance.M;

	unsigned int solution[SOLUTION_LENGTH]; // one solution
	unsigned int solutions[procs-1][SOLUTION_LENGTH]; // solutions sending buffer
	int objVec[M]; // one result
	int computedObjVecs[procs-1][RESULT_DIMENSION]; // solutions receiving buffer

	// Master
	if(self == PROC_NULL) {
		//-----------------------------------------------
		// SIGTERM manual interruption handling
		struct sigaction action;
		sigaction(SIGTERM, &action, NULL);

		std::vector<result_t> best_solutions;

		Gnuplot gnuplot("test");
		gnuplot.set_grid();


		while(true) {
			//-----------------------------------------------
			// build random solutions
			srand(time(NULL));
			for(unsigned i=0 ; i<(procs-1) ; i++)
				for(unsigned j=0 ; j<N ; j++)
					solutions[i][j] = (rand() / (double) RAND_MAX) < 0.5 ? 0 : 1;

			//-----------------------------------------------
			// process communications (distribute the tasks)
			// one array to be sent to each worker (procs-1 arrays)
			MPI_Scatter(solutions, N, MPI_INT, solution, N, MPI_INT, PROC_NULL, com); // TODO: avoid blocking collective communication primitives

			// retrieve computed objective vectors
			MPI_Gather(objVec, M, MPI_INT, computedObjVecs, M, MPI_INT, PROC_NULL, com);

			//-----------------------------------------------
			// filter solutions, keep the best ones // TODO: do this between scatter and gather (during the evaluation process) -> need two solutions and objVect arrays
			unsigned process_iterator=0;
			result_t new_solution;
			// the first time, keep the first solution
			if(best_solutions.size() == 0) {
				for(unsigned i=0; i<N; i++)
					new_solution.input[i] = solutions[0][i]; // can't use std::copy here...)
				for(unsigned i=0; i<M; i++)
					new_solution.output[i] = computedObjVecs[0][i];
				best_solutions.push_back(new_solution);
				process_iterator++;
			}

			for( ; process_iterator<(procs-1) ; process_iterator++) {
				bool keep_it = true;
				unsigned i = 0;
				// find out if the solution is worth being kept
				while(i < best_solutions.size() && keep_it) {
					if(best_solutions.at(i).output[0] >= computedObjVecs[process_iterator][0]) // TODO: remove hard coded values
						if(best_solutions.at(i).output[1] >= computedObjVecs[process_iterator][1])
							keep_it = false;
					i++;
				}
				// erase outdated values
				if(keep_it) {
					i = 0;
					while(i < best_solutions.size()) {
						if(best_solutions.at(i).output[0] <= computedObjVecs[process_iterator][0])
							if(best_solutions.at(i).output[1] <= computedObjVecs[process_iterator][1])
								best_solutions.erase(best_solutions.begin()+i);
						i++;
					}
					// save new optimal solution
					for(unsigned i=0; i<N; i++)
						new_solution.input[i] = solutions[process_iterator][i]; // can't use std::copy here...
					for(unsigned i=0; i<M; i++)
						new_solution.output[i] = computedObjVecs[process_iterator][i];
					best_solutions.push_back(new_solution);
				}
			}

			//-----------------------------------------------
			// plot the results // TODO: do this between scatter and gather (during the evaluation process) or in another process
			std::vector<int> x,y;
			int x_min = INT_MAX, x_max = INT_MIN, y_min = INT_MAX, y_max = INT_MIN;
			for(unsigned int i=0; i<best_solutions.size(); i++) {
				x.push_back(best_solutions.at(i).output[0]); // TODO: clean that (faster way to transmit points to gnuplot?)
				y.push_back(best_solutions.at(i).output[1]);
			}
			gnuplot.reset_plot();
			gnuplot.remove_tmpfiles();
			gnuplot.plot_xy(x, y, "best results");
			
			//-----------------------------------------------
			// ouput : print the solution and its objective vector (TEMPORARY)
			/*
			for(unsigned int i=0; i<procs-1 ; i++) {
				printf("result for process %d: ", i);
				for(unsigned int j=0; j<M ; j++)
					printf("%d ", computedObjVecs[i][j]);
				printf("%d\ncorresponding solution: ", N);
				for(unsigned int j=0; j<N; j++)
					printf("%d", solutions[i][j]);
				printf("\n\n");
			}
			*/
			printf("best solutions so far:\n");
			for(unsigned int i=0; i<best_solutions.size(); i++)
				printf("%d : (%d ; %d)\n",i,best_solutions.at(i).output[0],best_solutions.at(i).output[1]);
			/* sleep(1); */
		}
	}

	//-----------------------------------------------
	// slaves: evaluate your solution
	if(self != PROC_NULL) {
		while(true) {
			MPI_Scatter( solutions, N, MPI_INT, solution, N, MPI_INT, PROC_NULL, com);
			eval(&instance,(int*) solution,objVec); // can't use mubqp.eval (we need arrays instead of vectors)

			// send results
			MPI_Gather( objVec, M, MPI_INT, computedObjVecs, M, MPI_INT, PROC_NULL, com);
		}
	}

  return 0;
}
