/*
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; version 3 
    of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Contact: http://mocobench.sourceforge.net

Authors:
    Arnaud Liefooghe <arnaud.liefooghe@lifl.fr>
    Sebastien Verel  <sebastien.verel@inria.fr>
*/

/*
 * To compile this soure code with gcc:
 *   g++ -o test_MUBQPEval test_MUBQPEval.cpp -I.
 *
 * Execution:
 *   ./test_MUBQPEval instance.dat 
 * where instance.dat is the file name of the mUBQP problem instance
 *
 * This simple test performs:
 *  - read a mUBQP problem instance, 
 *  - initialize a random solution, 
 *  - evaluate the solution
 *  - print the result
 */

#include <iostream>
#include <vector>
#include <string>
#include <time.h>

#include <mubqpEval.h>

/***************************************************************/
/*                 Main                                        */
/***************************************************************/

int main(int argc, char *argv[]) {
  //-----------------------------------------------
  // evaluation function of the rMNK-landscapes

  if (argc != 2) {
    std::cerr << "Invalid number of parameters. Use the command line: ./test_MUBQPEval instance.dat" << std::endl;
    exit(1);
  }

  MUBQPEval mubqp(argv[1]);

  // size of the bit string
  unsigned N = mubqp.getN();
  // objective space dimension
  unsigned M = mubqp.getM();

  //-----------------------------------------------
  // initialization of the solution

  srand(time(NULL));

  std::vector<bool> solution(N);

  for(unsigned i = 0; i < N; i++)
    solution[i] = (rand() / (double) RAND_MAX) < 0.5 ? true : false;
  
  //-----------------------------------------------
  // eval the solution

  std::vector<double> objVec;

  mubqp.eval(solution, objVec);

  //-----------------------------------------------
  // ouput : print the solution and its objective vector

  for(unsigned i = 0; i < M ; i++)
    std::cout << objVec[i] << " ";
  
  std::cout << solution.size() << " " ;
  
  for(unsigned i = 0; i < N; i++)
    std::cout << solution[i] ;
  
  std::cout << std::endl;

  return 0;
}

