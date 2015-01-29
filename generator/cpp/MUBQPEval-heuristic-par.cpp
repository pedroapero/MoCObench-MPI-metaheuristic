#include <iostream>
#include <vector>
#include <string>
#include <climits>
#include <time.h>
#include <unistd.h>

#include "signal.h"
#include "mpi.h"

#include "mubqpEval-C-version.c"
// #include "gnuplot_i.hpp"


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


// the tuple solution - vector
struct result_t {
	std::vector<unsigned int> input; // tested matrix
	std::vector<int> output; // objectif vector
	int done; // 1 if we have evaluated all its neighbors, 0 else (MPI_BOOL does not exist) // TODO: refactor to unsigned int (or short?)
	int flipped; // used to mark input as a neighbour of another (with the flipped'nth bit flipped) // TODO: refactor to unsigned short int
	int solution_number;
};


std::vector<result_t>* _best_solutions; // needs to be global: it will be used in sigaction handler, and it's not possible to pass a parameter to a signal handler

void master_interruption_signal_handler(int sig) {
	/*
	fflush(stdout);
	std::cout << "stopping processes..." << std::endl; // TODO: display efficiency informations

	//------------------------------------
	// efficiency evaluation
	std::vector<int> best_objectives; // one result
	long double sum_best_objectives = 0;
	std::vector<int> worst_objectives; // one result
	long double sum_worst_objectives = 0;

	
	for(unsigned int i=0; i<_best_solutions->at(0).output.size(); i++) {
		best_objectives.push_back(INT_MIN);
		worst_objectives.push_back(INT_MAX);
	}
	
	for(unsigned int i=0; i<_best_solutions->size(); i++) {
		for(unsigned int j=0; j<best_objectives.size(); j++) {
			if(_best_solutions->at(i).output.at(j) > best_objectives.at(j))
				best_objectives.at(j) = _best_solutions->at(i).output.at(j);
			if(_best_solutions->at(i).output.at(j) < worst_objectives.at(j))
				worst_objectives.at(j) = _best_solutions->at(i).output.at(j);
		}
	}
	
	for(unsigned int i=0; i<best_objectives.size(); i++) {
		sum_best_objectives += best_objectives.at(i);
		sum_worst_objectives += worst_objectives.at(i);
	}


	//------------------------------------
	// display efficiency informations
	std::cout << "number of solutions found: " << _best_solutions->size() << std::endl;
	std::cout << "(best_objectives + worst_objectives) / number of solutions: " << (sum_best_objectives + sum_worst_objectives) / _best_solutions->size() << std::endl;
	*/

	//------------------------------------
	// display best objectif vectors (to be redirected to instance file)
	for(unsigned int i=0; i<_best_solutions->size(); i++)
		std::cout << _best_solutions->at(i).output.at(0) << " " << _best_solutions->at(i).output.at(1) << std::endl; // WARN: only for two dimensionnal objectif vectors

	MPI_Abort(com,EXIT_SUCCESS);
	MPI_Finalize();
	exit(sig);
}

void slaves_interruption_signal_handler(int sig) {
	// do nothing, wait for the master to kill you
	sleep(100);
}


//-----------------------------------------------
// generate a solution
void generate_solution(std::vector<unsigned int> &solution) {
	for(unsigned int i=0 ; i<solution.size() ; i++)
		solution.at(i) = ((rand() / (double) RAND_MAX) < 0.5 ? 0 : 1);
}

void save_solution(std::vector<unsigned int> input, std::vector<int> output, int flipped, int solution_number, std::vector<result_t> &best_solutions) {
	result_t new_solution;

	new_solution.input = input;
	new_solution.output = output;
	new_solution.flipped = flipped;
	new_solution.solution_number = solution_number;
	new_solution.done = 0;

	best_solutions.push_back(new_solution);
}

// returns 1 if objVec1 dominates objVec2 ; -1 if objVec2 dominates objVec1 ; 2 in case of error ; 0 else (not comparable)
int compare_vectors(std::vector<int> objVec1, std::vector<int> objVec2) {
	if(objVec1.size() != objVec2.size()) {
		printf("warning: comparing two vectors of different sizes!\n");
		return 2;
	}

	for(unsigned int i=0; i < objVec1.size(); i++) {
		if(objVec1.at(i) > objVec2.at(i)) { // does 1 dominates 2?
			for(i = i+1; i<objVec1.size(); i++)
				if(objVec1.at(i) < objVec2.at(i)) return 0; // no
			return 1; // yes
		}
		if(objVec2.at(i) > objVec1.at(i)) { // does 2 dominates 1?
			for(i = i+1; i<objVec1.size(); i++)
				if(objVec2.at(i) < objVec1.at(i)) return 0; // no
			return -1; // yes
		}
	}
	return 0; // comparing same vector!
}

//-----------------------------------------------
// filter non optimal solutions, save the best ones.
void filter_solutions(std::vector<unsigned int> solution, std::vector<int> objVec, int flipped, int solution_number, std::vector<result_t> &best_solutions) {
	unsigned int i = 0;
	for(unsigned int i=0; i < best_solutions.size(); i++) {
		int comparison = compare_vectors(objVec, best_solutions.at(i).output);
		if(comparison == -1) return; // new objVec is dominated
		else if(comparison == 1) best_solutions.erase(best_solutions.begin() + i); // the new vector dominates a best_solution, it will inevitably be saved
	}
	// so objVec dominates or equals best_solutions => save it
	save_solution(solution, objVec, flipped, solution_number, best_solutions); // will keep the solution if best_solutions is empty
}


//-----------------------------------------------
// Display best objectif vectors
void display_vectors(std::vector<result_t> best_solutions) {
	printf("\n\n");
	for(unsigned int i=0; i<best_solutions.size(); i++)
		printf("%d %d %d\n",best_solutions.at(i).output.at(0),best_solutions.at(i).output.at(1), best_solutions.at(i).done);
	printf("(%d solutions)\n", (int)best_solutions.size());
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


  if (argc != 3) {
    std::cerr << "Invalid number of parameters.\nUsage: mpirun -np [number of processes} ./MUBQPEval-heuristic-par instance.dat [pool size]" << std::endl;
    exit(1);
  }

	if (procs < 2) {
		printf("Two processes minimum, it's master/slaves schema.");
		exit(1);
	}

	// using the C version (more practical to use unsigned integers instead of booleans with MPI)
	Instance instance;
	loadInstance(argv[1], &instance);
  unsigned int N = instance.N; // length of a solution
  unsigned int M = instance.M; // dimension of objVec

	std::vector<unsigned int> solution(N); // one solution
	std::vector<int> objVec(M); // one result
	std::vector<result_t> best_solutions; // solutions-results storage structure
	_best_solutions = &best_solutions;

	int pool_size = atoi(argv[2]); // number of solutions to sent to each worker
	std::vector<std::vector<std::vector<unsigned int> > >
		sent_solutions(procs, std::vector<std::vector<unsigned int> >
				(pool_size, std::vector<unsigned int>(M))
				); // solutions whose neighbours are being evaluated (pool_size solutions per process => three dimensionnal vector)
	// note: allocates procs*pool_size + procs*pool_size*M => is it still contiguous in memory?

	int position; // buffer cursor


	/***************************************************************/
	/*                 Master                                       */
	/***************************************************************/
	if(self == PROC_NULL) {
		//-----------------------------------------------
		// SIGUSR1 manual interruption handling
		struct sigaction action;
		action.sa_handler = master_interruption_signal_handler;
		sigaction(SIGUSR1, &action, NULL); // TODO: this is overwritten by MPI's own handler

		/*
		Gnuplot gnuplot("results");
		gnuplot.set_grid();
		*/

		//-----------------------------------------------
		// initialization
		srand(time(NULL));
		while(best_solutions.size() < procs-1) {
			generate_solution(solution);
			eval(&instance,(int*) solution.data(), objVec.data());
			filter_solutions(solution, objVec, 0, 0, best_solutions);
		}

		//-----------------------------------------------
		// process communications (distribute the tasks)
		// one random solution to be sent to each worker
		for(unsigned int p=1; p<procs; p++) {
			// store solution to be sent
			sent_solutions.at(p).at(0) = best_solutions.at(p-1).input;
			best_solutions.at(p-1).done = 1;

			int buffer_size = (N * sizeof(int));
			char buffer[buffer_size];
			int position = 0;
			MPI_Pack(sent_solutions.at(p).at(0).data(), N, MPI_UNSIGNED, buffer, buffer_size, &position, com); // we will later send pools of solutions to workers in MPI_PACKED.
			MPI_Send(buffer, position, MPI_PACKED, p, 0, com);
		}

		while(true) {
			// output
			// display_vectors(best_solutions);
			
			//-----------------------------------------------
			// receive results
			// allocate reception buffer
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, com, &status);
			int count;
			MPI_Get_count(&status, MPI_PACKED, &count);
			int number_of_results = count / ((1 + 1 + M) * sizeof(int)); // one result is a triplet (solution_number ; flipped bit number ; output )
			char reception_buffer[count];

			// receive from process we just probed
			MPI_Recv(reception_buffer, count, MPI_PACKED, status.MPI_SOURCE, MPI_ANY_TAG, com, &status);

			// rebuild solutions and filter
			std::vector<unsigned int> input(N);
			std::vector<int> output(M);
			int flipped;
			int solution_number;
			position = 0;
			for(unsigned int i=0; i<number_of_results; i++) { // one iteration per couple: flipped, output (correspond to one solution)
				MPI_Unpack(reception_buffer, count, &position, &solution_number, 1, MPI_INT, com);
				MPI_Unpack(reception_buffer, count, &position, &flipped, 1, MPI_INT, com);
				MPI_Unpack(reception_buffer, count, &position, output.data(), M, MPI_INT, com);
				input = sent_solutions.at(status.MPI_SOURCE).at(solution_number);
				input.at(flipped) = (input.at(flipped) == 0 ? 1 : 0); // retrieve the input by reflipping the bit of the original solution
				filter_solutions(input, output, flipped, solution_number, best_solutions);
			}

			//----------------------------------------------
			// send next seeds
			unsigned int solutions_iterator = 0;

			// initialize the pool of solutions to send
			int buffer_size = (N * pool_size * sizeof(int));
			char sending_buffer[buffer_size];
			int position = 0;
			sent_solutions.at(status.MPI_SOURCE).clear(); // new size is reset to 0
			while(solutions_iterator < best_solutions.size() && sent_solutions.at(status.MPI_SOURCE).size()<pool_size ) {
				if(best_solutions.at(solutions_iterator).done == 0) {
					best_solutions.at(solutions_iterator).done = 1;
					sent_solutions.at(status.MPI_SOURCE).push_back(best_solutions.at(solutions_iterator).input);

					MPI_Pack(sent_solutions.at(status.MPI_SOURCE).back().data(), N, MPI_UNSIGNED, sending_buffer, buffer_size, &position, com);
				}
				solutions_iterator++;
			}

			// should only happen for small N!
			if(sent_solutions.at(status.MPI_SOURCE).size() == 0) {
				printf("no more solution to evaluate!\n");
				sleep(100);
				master_interruption_signal_handler(0);
			}

			// send the pool
			MPI_Send(sending_buffer, position, MPI_PACKED, status.MPI_SOURCE, 0, com);

			
			//-----------------------------------------------
			// plot the results
			/*
			std::vector<int> x,y;
			int x_min = INT_MAX, x_max = INT_MIN, y_min = INT_MAX, y_max = INT_MIN;
			for(unsigned int i=0; i<best_solutions.size(); i++) {
				x.push_back(best_solutions.at(i).output[0]); // TODO: clean that (faster way to transmit points to gnuplot?)
				y.push_back(best_solutions.at(i).output[1]); // only for two dimensionnal objVecs!
			}
			gnuplot.reset_plot();
			gnuplot.remove_tmpfiles();
			gnuplot.plot_xy(x, y, "best results");
			*/
		}
	}

	/***************************************************************/
	/*                 Slaves (browse and evaluate neighbors)      */
	/***************************************************************/
	if(self != PROC_NULL) {
		clock_t computation_beginning_timestamp;
		clock_t solution_generation_beginning_timestamp = clock();

		//-----------------------------------------------
		// SIGUSR1 manual interruption handling
		struct sigaction action;
		action.sa_handler = slaves_interruption_signal_handler;
		sigaction(SIGUSR1, &action, NULL); // not possible to overwrite mpirun's own SIGINT handler

		while(true) {
			best_solutions.clear();
			
			MPI_Probe(PROC_NULL, MPI_ANY_TAG, com, &status);
			int count;
			MPI_Get_count(&status, MPI_PACKED, &count);
			int number_of_seeds = count / (N * sizeof(int));
			char reception_buffer[count];

			//-----------------------------------------------
			// receive next seeds
			MPI_Recv(reception_buffer, count, MPI_PACKED, status.MPI_SOURCE, MPI_ANY_TAG, com, &status);
			// printf("%d have waited for new seeds for %lfs\n", self, (double) (clock() - solution_generation_beginning_timestamp) / CLOCKS_PER_SEC);
			computation_beginning_timestamp = clock();

			// iterate over seeds
			position = 0;
			for(unsigned int i=0; i<number_of_seeds; i++) {
				MPI_Unpack(reception_buffer, count, &position, solution.data(), N, MPI_UNSIGNED, com);
				
				//-----------------------------------------------
				// iterate through neighbours (N neighbours for a solution of length N)
				for(unsigned int solution_cursor=0; solution_cursor<N; solution_cursor++) {
					solution.at(solution_cursor) = (solution.at(solution_cursor) == 0 ? 1 : 0); // flip solution_cursor'nth bit

					eval(&instance,(int*) solution.data(), objVec.data());
					filter_solutions(solution, objVec, solution_cursor, i, best_solutions); // i is the solution_number
					solution.at(solution_cursor) = (solution.at(solution_cursor) == 0 ? 1 : 0); // reflip back to the original bit
				}
			}

			// send results
			unsigned int buffer_size = best_solutions.size() * (1 + 1 + M) ; // total number of integers to send: flipped number, solution number, objVec
			buffer_size *= sizeof(int);
			char sending_buffer[buffer_size];
			position = 0;
			for(unsigned int i=0; i<best_solutions.size(); i++) {
				MPI_Pack(&best_solutions.at(i).solution_number, 1, MPI_INT, sending_buffer, buffer_size, &position, com);
				MPI_Pack(&best_solutions.at(i).flipped, 1, MPI_INT, sending_buffer, buffer_size, &position, com);
				MPI_Pack(best_solutions.at(i).output.data(), M, MPI_INT, sending_buffer, buffer_size, &position, com);
			}
			// printf("%d: total execution time (minus Probe and reception buffer allocations): %.2fs\n", self, (double) (clock() - computation_beginning_timestamp) / CLOCKS_PER_SEC);
			MPI_Send(sending_buffer, position, MPI_PACKED, PROC_NULL, 0, com);
			solution_generation_beginning_timestamp = clock();
		}
	}

  return 0;
}
