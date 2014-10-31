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

		//-----------------------------------------------
		// build random solutions
		srand(time(NULL));
		for(unsigned i=0 ; i<(procs-1) ; i++)
			for(unsigned j=0 ; j<N ; j++)
				solutions[i][j] = (rand() / (double) RAND_MAX) < 0.5 ? 0 : 1;

		//-----------------------------------------------
		// process communications (distribute the tasks)
		// one array to be sent to each worker (procs-1 arrays)
		MPI_Scatter(solutions, N, MPI_INT, solution, N, MPI_INT, PROC_NULL, com);

		// retrieve computed objective vectors
		MPI_Gather(objVec, M, MPI_INT, computedObjVecs, M, MPI_INT, PROC_NULL, com);

		//-----------------------------------------------
		// ouput : print the solution and its objective vector (TEMPORARY)
		for(unsigned int i=0; i<procs-1 ; i++) {
			printf("result for process %d: ", i);
			for(unsigned int j=0; j<M ; j++)
				printf("%d ", computedObjVecs[i][j]);
			printf("%d\ncorresponding solution: ", N);
			for(unsigned int j=0; j<N; j++)
				printf("%d", solutions[i][j]);
			printf("\n\n");
		}
	}

	//-----------------------------------------------
	// slaves: evaluate your solution
	if(self != PROC_NULL) {
		MPI_Scatter( solutions, N, MPI_INT, solution, N, MPI_INT, PROC_NULL, com);
		eval(&instance,(int*) solution,objVec); // can't use mubqp.eval (we need arrays instead of vectors)

		// send results
		MPI_Gather( objVec, M, MPI_INT, computedObjVecs, M, MPI_INT, PROC_NULL, com);
	}

  return 0;
}
