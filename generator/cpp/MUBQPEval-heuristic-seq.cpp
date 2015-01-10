#include <iostream>
#include <vector>
#include <string>
#include <time.h>
#include <signal.h>
#include <climits>

#include "mubqpEval-C-version.c"
#include "gnuplot_i.hpp"

void interruption_signal_handler(int sig) {
	std::cout << "Ending process..." << std::endl; // TODO: display efficiency informations
	exit(sig);
}

// the tuple solution - vector
struct result_t {
	std::vector<unsigned int> input; // tested matrix
	std::vector<int> output; // objectif vector
	int done; // 1 if we have evaluated all its neighbors, 0 else (MPI_BOOL does not exist) // TODO: refactor to unsigned int (or short?)
};

void save_solution(std::vector<unsigned int> input, std::vector<int> output, std::vector<result_t> &best_solutions) {
	result_t new_solution;

	new_solution.input = input;
	new_solution.output = output;
	new_solution.done = false;

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
void filter_solutions(std::vector<unsigned int> solution, std::vector<int> objVec, std::vector<result_t> &best_solutions) {
	unsigned int i = 0;
	for(unsigned int i=0; i < best_solutions.size(); i++) {
		int comparison = compare_vectors(objVec, best_solutions.at(i).output);
		if(comparison == -1) return; // new objVec is dominated
		else if(comparison == 1) best_solutions.erase(best_solutions.begin() + i); // the new vector dominates a best_solution, it will inevitably be saved
	}
	// so objVec dominates or equals best_solutions => save it
	save_solution(solution, objVec, best_solutions); // will keep the solution if best_solutions is empty
}


/***************************************************************/
/*                 Main                                        */
/***************************************************************/

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Invalid number of parameters. \nUsage: ./test_MUBQPEval instance.dat" << std::endl;
    exit(1);
  }

	struct sigaction action;
	action.sa_handler = interruption_signal_handler;
	sigaction(SIGTERM, &action, NULL);

	// using the C version for compatibility with MPI and parallel version
	Instance instance;
	loadInstance(argv[1], &instance);
  unsigned int N = instance.N; // length of a solution
  unsigned int M = instance.M; // dimension of objVec

	std::vector<unsigned int> solution(N); // one solution
	std::vector<int> objVec(M); // one result
	std::vector<result_t> best_solutions; // solutions-results storage structure

	Gnuplot gnuplot("test");
	gnuplot.set_grid();

  //-----------------------------------------------
  // initializing the seed
  srand(time(NULL));
  for(unsigned int i = 0; i < N; i++)
    solution.at(i) = (rand() / (double) RAND_MAX) < 0.5 ? 0 : 1;
  
  //-----------------------------------------------
  // evaluate and store the first solution
	eval(&instance,(int*) solution.data(), objVec.data());
	filter_solutions(solution, objVec, best_solutions);

  //-----------------------------------------------
  // main loop: iterate through best_solutions
	int solution_iterator = 0;
	while(solution_iterator < best_solutions.size()) {
		if(best_solutions.at(solution_iterator).done == 0) {
			best_solutions.at(solution_iterator).done = 1;
			solution = best_solutions.at(solution_iterator).input;

			//-----------------------------------------------
			// iterate through neighbours (N neighbours for a solution of length N)
			for(unsigned int solution_cursor=0; solution_cursor<N; solution_cursor++) {
					solution.at(solution_cursor) = (solution.at(solution_cursor) == 0 ? 1 : 0); // flip solution_cursor'nth bit

					eval(&instance,(int*) solution.data(), objVec.data());
					filter_solutions(solution, objVec, best_solutions); // i is the solution_number

					solution.at(solution_cursor) = (solution.at(solution_cursor) == 0 ? 1 : 0); // reflip back to the original bit
			}

			solution_iterator = 0;
		}
		else
			solution_iterator++;


		//-----------------------------------------------
		// output: print best objective vectors
		printf("\n\n");
		for(unsigned int i=0; i<best_solutions.size(); i++)
			printf("%d %d %d\n",best_solutions.at(i).output.at(0),best_solutions.at(i).output.at(1), best_solutions.at(i).done);
		printf("(%d solutions)\n", best_solutions.size());


		//-----------------------------------------------
		// plot the results
		std::vector<int> x,y;
		int x_min = INT_MAX, x_max = INT_MIN, y_min = INT_MAX, y_max = INT_MIN;
		for(unsigned int i=0; i<best_solutions.size(); i++) {
			x.push_back(best_solutions.at(i).output[0]); // TODO: clean that (faster way to transmit points to gnuplot?)
			y.push_back(best_solutions.at(i).output[1]); // only for two dimensionnal objVecs!
		}
		gnuplot.reset_plot();
		gnuplot.remove_tmpfiles();
		gnuplot.plot_xy(x, y, "best results");
	}
	

  return 0;
}
