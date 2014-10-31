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

#ifndef __mubqpEval
#define __mubqpEval

/*
 *  Fitness function of the multiobjective unconstrained binary quadratic programming problem (mUBQP problem)
 *  reading a m UBQP problem instance file
 *  in c++ style
 *
 *  UBQP problem instances can be generated with the mubqpGenerator.R 
 *
 */

#include <fstream>
#include <vector>
#include <string>

using namespace std;

class MUBQPEval
{
public:
  /*
   *  Constructor 
   *
   * @param _fileName file name instance of the mUBQP
   */
  MUBQPEval(const char * _fileName) { 
    load(_fileName);
  }

  /*
   *  Destructor
   */
  ~MUBQPEval() {
    if (Qs != NULL) {
      for(int n = 0; n < M; n++) {
	for(int i = 0; i < N; i++) 
	  delete[] (Qs[n][i]);
	delete[] Qs[n];
      }

      delete Qs;
      Qs = NULL;
    }
  }

  /*
   * Compute the fitness function
   *
   * @param _solution the solution to evaluate
   * @param _objVec   the objective vector of the corresponding solution
   */
  void eval(std::vector<bool> & _solution, std::vector<double> & _objVec) {
    _objVec.resize(M);
      
    for(unsigned int n = 0; n < M; n++) 
      _objVec[n] = evalUBQP(n, _solution);
      
  } 

  /*
   * to get objective space dimension
   *
   * @return dimension of the objective space
   */
  unsigned getM()   { return M ; }

  /*
   * to get bitstring size
   *
   * @return dimension of the bitstring
   */
  unsigned getN()   { return N ; }

  /*
   * to get density parameter
   *
   * @return density parameter
   */
  double getDensity()   { return density ; }

  /*
   * to get the correlation between each tuple of contributions 
   *
   * @return parameter rho
   */
  double getRho() { return rho ; }

  /*
   * to get one matrix Q
   *
   * @param dimension 
   * @return matrix Q
   */
  int** getQ(int dimension) {
    return Qs[dimension];
  }


  /*
   * to get all the matrices Q
   *
   * @return matrices Q
   */
  int*** getQs() {
    return Qs;
  }


  void print() {
    std::cout << rho << " " << M << " " << N << " " << density << std::endl;

    for(unsigned int n = 0; n < M; n++) {
      for(unsigned int i = 0; i < N; i++) {
	for(unsigned int j = 0; j < N; j++) {
	  std::cout << Qs[n][i][j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
    }
  }

protected:
  // correlation between contributions
  double rho;

  // number of objective functions
  unsigned M;

  // size of the bit string
  unsigned N;

  // density of the matrices
  double density;

  // the M matrices 
  int *** Qs;

  /***********************************************
   *
   * load the file with data for an mUBQP problem
   *
   * @param fileName file name instance of the mUBQP problem
   *
   ***********************************************/
  virtual void load(const char * fileName) {
    fstream file;
    file.open(fileName, ios::in);

    if (file.is_open()) {
      string s;
	    
      // read the commentaries
      string line;

      file >> s;
      while (s[0] == 'c') {
	getline(file, line, '\n');
	file >> s;
      }
	    
      // read the parameters
      if (s.compare("p") != 0) 
	cerr <<"Error MUBQP.load: expected line beging by \"p\" at " << s << " in " << fileName << std::endl;
	    
      file >> s;

      if (s.compare("MUBQP") != 0) 
	cerr <<"Error MUBQP.load: type MUBQP expected at " << s << " in "  << fileName << std::endl;

      // effective read of the parameters
      file >> rho >> M >> N >> density;

      init();

      // read of matrices
      file >> s;
      if (s.compare("p") != 0) 
	cerr <<"Error MUBQP.load: expected line beging by \"p\" at " << s << " in " << fileName << std::endl;

      file >> s;

      if (s.compare("matrices") == 0) 
	loadQs(file);
      else 
	cerr <<"Error MUBQP.load: line with \"matrices\" expected at " << s << " in " << fileName << std::endl;

      file.close();
    } else 
      cerr << "Error MUBQP.load: impossible to open file " << fileName << endl;
  };

  /***********************************************
   *
   * Initialization of the differents matrices
   *
   ***********************************************/
  void init() {
    Qs = new int**[M];
	
    for(unsigned n = 0; n < M; n++) {
      Qs[n] = new int*[N];
	    
      for(unsigned i = 0; i < N; i++) 
	Qs[n][i] = new int[N];
    }
  }

  /***********************************************
   *
   * load the matrices from file variable 
   *
   * @param file open file of the matrices
   *
   ***********************************************/
  void loadQs(fstream & file) {
    unsigned n, i, j;

    for(j = 0; j < N; j++)
      for(i = 0; i < N; i++)
	for(n = 0; n < M; n++)
	  file >> Qs[n][i][j];

    // put the matrix in lower triangular form
    for(n = 0; n < M; n++)
      for(unsigned i = 1; i < N; i++)
	for(unsigned int j = 0; j < i; j++) {
	  Qs[n][i][j] = Qs[n][i][j] + Qs[n][j][i];
	  Qs[n][j][i] = 0;
	}
  }

  /***********************************************
   *
   *  fitness function of a single objective UBQP problem
   *
   * @param _numObj the objective fonction to consider
   * @param _sol the solution to evaluate
   *
   ***********************************************/
  int evalUBQP(unsigned _numObj, std::vector<bool> & _sol){
    int fit = 0;
    unsigned int j;

    for(unsigned i = 0; i < N; i++)
      if (_sol[i] == 1) 
	for(j = 0; j <= i; j++)
	  if (_sol[j] == 1) 
	    fit += Qs[_numObj][i][j];
	
    return fit ;
  }

}; 

#endif

