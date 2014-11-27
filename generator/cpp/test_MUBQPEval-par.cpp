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

// the tuple solution - vector
struct result_t {
	// TODO: remove hard coded values!! (input is size N, below output is M, but impossible to use them to fix array size)
	std::vector<unsigned int> input; // tested matrix
	std::vector<int> output; // objectif vector
	unsigned short int done; // 1 if we have evaluated all its neighbors, 0 else (MPI_BOOL does not exist)
} result;


//-----------------------------------------------
// generate a solution
void generate_solution(std::vector<unsigned int> &solution, unsigned int size) {
	if(solution.size() != 0) solution.clear();
	for(unsigned int i=0 ; i<size ; i++)
		solution.push_back((rand() / (double) RAND_MAX) < 0.5 ? 0 : 1);
}

void save_solution(std::vector<unsigned int> input, std::vector<int> output, unsigned short done, std::vector<result_t> &best_solutions) {
	result_t new_solution;

	new_solution.input = input;
	new_solution.output = output;
	new_solution.done = done;

	best_solutions.push_back(new_solution);
}

// returns 1 if objVec1 dominates objVec2 ; -1 if objVec2 dominates objVec1 ; 2 in case of error ; 0 else (not comparable)
int compare_vectors(std::vector<int> objVec1, std::vector<int> objVec2) {
	if(objVec1.size() != objVec2.size()) return 2;

	for(unsigned int i=0; i < objVec1.size(); i++) {
		if(objVec1.at(i) > objVec2.at(i)) { // does 1 dominates 2?
			i++;
			while(i < objVec1.size()) {
				if(objVec1.at(i) < objVec2.at(i)) return 0;
				i++;
			}
			return 1;
		}
		if(objVec2.at(i) > objVec1.at(i)) { // does 2 dominates 1?
			i++;
			while(i < objVec1.size()) {
				if(objVec2.at(i) < objVec1.at(i)) return 0;
				i++;
			}
			return -1;
		}
	}
	return 0; // comparing same vector!
}

//-----------------------------------------------
// filter non optimal solutions, save the best ones.
void filter_solutions(std::vector<unsigned int> solution, std::vector<int> objVec, unsigned short done, std::vector<result_t> &best_solutions) {
	unsigned int i = 0;
	// find out if the solution is worth being kept
	for(unsigned int i=0; i < best_solutions.size(); i++) {
		int comparison = compare_vectors(objVec, best_solutions.at(i).output);
		if(comparison == -1) return;
		else if(comparison == 1) best_solutions.erase(best_solutions.begin() + i); // the new solution will inevitably be saved in this case
		i++;
	}
	// so objVec dominates or equals best_solutions => save it
	save_solution(solution, objVec, done, best_solutions); // will keep the solution if best_solutions is empty
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

	// using the C version
	Instance instance;
	loadInstance(argv[1], &instance);
  unsigned int N = instance.N;
  unsigned int M = instance.M;

	std::vector<unsigned int> solution(N); // one solution
	std::vector<int> objVec(M); // one result
	std::vector<result_t> best_solutions; // solutions-results storage structure
	int position; // buffer cursor


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
			eval(&instance,(int*) solution.data(), objVec.data()); // can't use mubqp.eval (we need arrays instead of vectors)
			filter_solutions(solution, objVec, 1, best_solutions);
		}

		//-----------------------------------------------
		// process communications (distribute the tasks)
		// one random solution to be sent to each worker
		for(unsigned int p=1; p<procs; p++)
			MPI_Send(best_solutions[p-1].input.data(), N, MPI_UNSIGNED, p, 0, com);

		while(true) {
			//-----------------------------------------------
			// output: print best objective vectors
			printf("\n\n");
			for(unsigned int i=0; i<best_solutions.size(); i++)
				printf("%d %d %d\n",best_solutions.at(i).output.at(0),best_solutions.at(i).output.at(1), best_solutions.at(i).done);
		
			
			//-----------------------------------------------
			// receive results
			// allocate reception buffer
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, com, &status);
			int count;
			MPI_Get_count(&status, MPI_PACKED, &count);
			int number_of_results = count / ((N + M + 1) * sizeof(int)); // one result is a triplet (input ; output ; done)
			char buffer[count];

			// receive from process we just probed
			MPI_Recv(buffer, count, MPI_PACKED, status.MPI_SOURCE, MPI_ANY_TAG, com, &status);

			// rebuild solutions and filter
			std::vector<unsigned int> input(N);
			std::vector<int> output(M);
			int done;
			position = 0;
			for(unsigned int i=0; i<number_of_results; i++) { // one iteration per triplet (input, output, done)
				MPI_Unpack(buffer, count, &position, input.data(), N, MPI_INT, com);
				MPI_Unpack(buffer, count, &position, output.data(), M, MPI_INT, com);
				MPI_Unpack(buffer, count, &position, &done, 1, MPI_INT, com);
				filter_solutions(input, output, done, best_solutions);
			}

			//----------------------------------------------
			// send next seed
			bool sent = false;
			unsigned solutions_iterator = 0;
			while(solutions_iterator < best_solutions.size() && !sent) {
				if(best_solutions.at(solutions_iterator).done == 0) {
					//printf("sending old solution to explore\n");
					MPI_Send(best_solutions.at(solutions_iterator).input.data(), N, MPI_UNSIGNED, status.MPI_SOURCE, 0, com); // TODO: maybe do something with the tag? (0 else...)
					best_solutions.at(solutions_iterator).done = 1;
					sent = true;
				}
				solutions_iterator++;
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
		}
	}

	/***************************************************************/
	/*                 Slaves (browse and evaluate neighbors)      */
	/***************************************************************/
	if(self != PROC_NULL) {
		while(true) {
			best_solutions.clear();
			
			// next seed
			MPI_Recv(solution.data(), N, MPI_UNSIGNED, PROC_NULL, MPI_ANY_TAG, com, &status);

			// iterate through neighbours (N neighbours for a solution of length N)
			for(unsigned int solution_cursor=0; solution_cursor<N; solution_cursor++) {
				solution.at(solution_cursor) = (solution.at(solution_cursor) == 0 ? 1 : 0); // flip solution_cursor'nth bit

				eval(&instance,(int*) solution.data(), objVec.data()); // can't use mubqp.eval (we need arrays instead of vectors)
				// mubqp.eval(solution, objVec);

				filter_solutions(solution, objVec, 0, best_solutions);
				solution.at(solution_cursor) = (solution.at(solution_cursor) == 0 ? 1 : 0); // reflip bit back to the original one
			}

			// send results
			unsigned int buffer_size = best_solutions.size() * (N + M + 1) ; // total number of integers to send (includes "done" flags)
			buffer_size *= sizeof(int);
			char buffer[buffer_size];
			position = 0;
			for(unsigned int i=0; i<best_solutions.size(); i++) {
				MPI_Pack(best_solutions.at(i).input.data(), N, MPI_INT, buffer, buffer_size, &position, com);
				MPI_Pack(best_solutions.at(i).output.data(), M, MPI_INT, buffer, buffer_size, &position, com);
				MPI_Pack(&best_solutions.at(i).done, 1, MPI_INT, buffer, buffer_size, &position, com);
			}
			MPI_Send(buffer, position, MPI_PACKED, PROC_NULL, 0, com);
		}
	}

  return 0;
}
