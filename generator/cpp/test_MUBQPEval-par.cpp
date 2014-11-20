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
#include <cstddef>

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


// TODO: use real instance values instead of this dirty workaround
#define SOLUTION_LENGTH 1000 // N
#define RESULT_DIMENSION 2 // M

// the tuple solution - vector
struct result_t {
	// TODO: remove hard coded values!! (input is size N, below output is M, but impossible to use them to fix array size)
	unsigned int input[SOLUTION_LENGTH]; // tested matrix 
	int output[RESULT_DIMENSION]; // objvect
	unsigned short int done; // 1 if we have evaluated all its neighbors, 0 else (MPI_BOOL does not exist)
} result;


//-----------------------------------------------
// generate a solution
void generate_solution(unsigned int solution[SOLUTION_LENGTH], unsigned int size) {
	for(unsigned int i=0 ; i<size ; i++)
		solution[i] = (rand() / (double) RAND_MAX) < 0.5 ? 0 : 1;
}

void save_solution(std::vector<result_t> &best_solutions, unsigned int solution[SOLUTION_LENGTH], unsigned int solution_length, int objVec[RESULT_DIMENSION], unsigned int result_dimension, unsigned short done) {
	result_t new_solution;
	for(unsigned int i=0; i<solution_length; i++)
		new_solution.input[i] = solution[i];
	for(unsigned int i=0; i<result_dimension; i++)
		new_solution.output[i] = objVec[i];
	new_solution.done = done;

	best_solutions.push_back(new_solution);
}

//-----------------------------------------------
// filter non optimal solutions, save the best ones.
void filter_solutions(std::vector<result_t> & best_solutions, unsigned int solution[SOLUTION_LENGTH], unsigned int solution_length, int objVec[RESULT_DIMENSION], unsigned int result_dimension, unsigned short done) {
	// the first time, keep the first solution
	if(best_solutions.size() == 0) {
		save_solution(best_solutions, solution, solution_length, objVec, result_dimension, done);
		return;
	}
	bool keep_it = true;
	unsigned int i = 0;
	// find out if the solution is worth being kept
	while(i < best_solutions.size() && keep_it) {
		if(best_solutions.at(i).output[0] >= objVec[0]) // TODO: remove hard coded values
			if(best_solutions.at(i).output[1] >= objVec[1])
				keep_it = false;
		i++;
	}
	// erase outdated values
	if(keep_it) {
		i = 0;
		while(i < best_solutions.size()) {
			if(best_solutions.at(i).output[0] < objVec[0])
				if(best_solutions.at(i).output[1] < objVec[1])
					best_solutions.erase(best_solutions.begin()+i);
			i++;
		}
		// save new optimal solution
		save_solution(best_solutions, solution, solution_length, objVec, result_dimension, done);
	}
}

/***************************************************************/
/*                 Main                                        */
/***************************************************************/
int main(int argc, char *argv[]) {
  //-----------------------------------------------
	// MPI initialization
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
	int objVec[RESULT_DIMENSION]; // one result
	std::vector<result_t> best_solutions; // solutions-results storage structure

	//--------------------------------------------
	// MPIze the type result_t
	MPI_Datatype MPI_result_t;
	int array_of_blocklengths[3] = { SOLUTION_LENGTH , RESULT_DIMENSION , 1};
	MPI_Datatype array_of_types[3] = { MPI_UNSIGNED , MPI_INT , MPI_SHORT};
	MPI_Aint array_of_displacements[3];
	array_of_displacements[0] = offsetof(result_t, input);
	array_of_displacements[1] = offsetof(result_t, output);
	array_of_displacements[2] = offsetof(result_t, done);
	MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_result_t);
	MPI_Type_commit(&MPI_result_t); // TODO: free the type with MPI_Type_free at the end of execution


	/***************************************************************/
	/*                 Master                                       */
	/***************************************************************/
	if(self == PROC_NULL) {
		//-----------------------------------------------
		// SIGTERM manual interruption handling
		struct sigaction action;
		sigaction(SIGTERM, &action, NULL);


		Gnuplot gnuplot("test");
		gnuplot.set_grid();

		//-----------------------------------------------
		// initialization
		srand(time(NULL));
		while(best_solutions.size() < procs-1) {
			generate_solution(solution, N);
			eval(&instance,(int*) solution,objVec); // can't use mubqp.eval (we need arrays instead of vectors)
			filter_solutions(best_solutions, solution, N, objVec, M, 1);
		}

		//-----------------------------------------------
		// process communications (distribute the tasks)
		// one random solution to be sent to each worker
		for(unsigned int p=1; p<procs; p++)
			MPI_Send(best_solutions[p-1].input, N, MPI_UNSIGNED, p, 0, com);

		while(true) {
			//-----------------------------------------------
			// output: print best objective vectors
			printf("\n\n");
			for(unsigned int i=0; i<best_solutions.size(); i++)
				printf("%d %d %d\n",best_solutions.at(i).output[0],best_solutions.at(i).output[1], best_solutions.at(i).done);
		
			
			//-----------------------------------------------
			// receive results
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, com, &status);
			int count;
			MPI_Get_count(&status, MPI_result_t, &count);

			// allocate reception buffer
			// result_t solutions_buffer[count];
			std::vector<result_t> solutions_buffer(count);

			MPI_Recv(&solutions_buffer[0], count, MPI_result_t, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // Receiving from process we just probed

			// filter
			for(unsigned int i=0; i<count; i++)
				filter_solutions(best_solutions, solutions_buffer.at(i).input, N, solutions_buffer.at(i).output, M, 0);

			//----------------------------------------------
			// send next seed
			bool sent = false;
			unsigned solutions_iterator = 0;
			while(solutions_iterator < best_solutions.size() && !sent) {
				if(best_solutions.at(solutions_iterator).done == 0) {
					//printf("sending old solution to explore\n");
					MPI_Send(best_solutions.at(solutions_iterator).input, N, MPI_UNSIGNED, status.MPI_SOURCE, 0, com); // TODO: maybe do something with the tag? (0 else...)
					best_solutions.at(solutions_iterator).done = 1;
					sent = true;
				}
				solutions_iterator++;
			}

			/*
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
			*/		

			// sleep(1);
		}
	}

	/***************************************************************/
	/*                 Slaves (browse and evaluate neighbors)      */
	/***************************************************************/
	if(self != PROC_NULL) {
		while(true) {
			best_solutions.clear();
			
			// next seed
			MPI_Recv(solution, N, MPI_UNSIGNED, PROC_NULL, MPI_ANY_TAG, com, &status);

			// iterate through neighbours (N neighbours for a solution of length N)
			unsigned int solution_tmp[SOLUTION_LENGTH];
			for(unsigned int solution_cursor=0; solution_cursor<N; solution_cursor++) {
				for(unsigned int i=0; i<N; i++)
					solution_tmp[i] = solution[i];
				solution_tmp[solution_cursor] = solution[solution_cursor] == 0 ? 1 : 0;

				eval(&instance,(int*) solution_tmp, objVec); // can't use mubqp.eval (we need arrays instead of vectors)

				filter_solutions(best_solutions, solution_tmp, N, objVec, M, 0);
			}

			// send results
			MPI_Ssend(&best_solutions[0], best_solutions.size(), MPI_result_t, PROC_NULL, 0, com);
		}
	}

  return 0;
}
