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
 *  Fitness function of the multiobjective UBQP problem
 *  reading a mUBQP problem instance file
 *  in c style
 *
 *  mUBQP problem instances can be generated with the mubqpGenerator.R 
 *
 *
 *  version 1.0: 09/27/2012
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mubqpEval-C-version.h"

/***********************************************
 *
 * Load the file of a mubqp problem instance
 *
 * @param _fileName file name instance of the mubqp problem
 * @param _instance the structure of the instance
 *
 ***********************************************/
void loadInstance(const char * _fileName, struct Instance * _instance) {
  FILE* file = NULL;
  unsigned int MAX_LINE = 512; // maximum number of characters per line in the file
  
  file = fopen(_fileName, "r");
  
  if (file != NULL) {
    // read the commentaries
    char line[MAX_LINE];
    char initialChar ;

    initialChar = fgetc(file);

    if (initialChar == EOF) {
      printf("Error MUBQP.load: no line to read in %s\n", _fileName);
      exit(1);
    }

    while (initialChar == 'c') {
      if (fgets(line, MAX_LINE, file) == NULL) {
	printf("Error MUBQP.load: expected data information in %s\n", _fileName);
	exit(1);
      }
      initialChar = fgetc(file);
      if (initialChar == EOF) {
	printf("Error MUBQP.load: expected data information in %s after comments:\nc %s\n", _fileName, line);
	exit(1);
      }
    }
	    
    // read the parameters
    if (initialChar != 'p') {
      printf("Error MUBQP.load: expected line beging by \"p\" in %s after comments:\nc %s\n", _fileName, line);
      exit(1);
    }
	    
    char s[MAX_LINE];

    fscanf(file, "%s ", s);

    if (strcmp(s, "MUBQP") != 0) {
      printf("Error MUBQP.load: \"MUBQP\" expected in %s:\np %s\n", _fileName, s);
      exit(1);
    }

    // effective read of the parameters
    if (fscanf(file, "%lf %d %d %lf", &(_instance->rho), &(_instance->M), &(_instance->N), &(_instance->d)) != 4) {
      printf("Error MUBQP.load: error when reading the parameters after \"p MUBQP\" in %s:\n", _fileName);
      exit(1);
    }

    initInstance(_instance);

    // read the matrices
    if (fscanf(file, "\n%c %s\n", &initialChar, s) != 2) {
      printf("Error MUBQP.load: expected \"p matrices\" after the parameters in %s\n", _fileName);
      exit(1);
    }

    if (initialChar != 'p') {
      printf("Error MUBQP.load: expected line beging by \"p\"  after the parameters in %s\n", _fileName);
      exit(1);
    }

    if (strcmp(s, "matrices") == 0) 
      loadQs(file, _instance);
    else {
      printf("Error MUBQP.load: expected \"matrices\"  after \"p\" in %s\n", _fileName);
      exit(1);
    }

    fclose(file);
  } else { 
    printf("Error MUBQP.load: impossible to open file %s\n", _fileName);
    exit(1);
  }
}

/***********************************************
 *
 * Initialization of the differents tables and epistasis links
 *
 * @param _instance the structure of the instance
 *
 ***********************************************/
void initInstance(struct Instance * _instance) {
  _instance->Qs = (int***) malloc(_instance->M * sizeof(int**));

  unsigned int n , i;
  for(n = 0; n < _instance->M; n++) {
    _instance->Qs[n] = (int**) malloc(_instance->N * sizeof(int*));

    for(i = 0; i < _instance->N; i++) {
      _instance->Qs[n][i] = (int*) malloc(_instance->N * sizeof(int));
    }
  }
}

/***********************************************
 *
 * Load the matrices
 *
 * @param _file open file of the instance
 * @param _instance the structure of the instance
 *
 ***********************************************/
void loadQs(FILE* _file, struct Instance * _instance) {
  unsigned m, i;
  int j;

  for(i = 0; i < _instance->N; i++)
    for(j = 0; j < _instance->N; j++)
      for(m = 0; m < _instance->M; m++)
	if (fscanf(_file, "%d", &(_instance->Qs[m][i][j])) != 1) {
	  printf("Error MUBQPEval.load: when reading the matrices, a number should be missing.\n");
	  exit(1);
	}

  // put the matrix in lower triangular form
  for(m = 0; m < _instance->M; m++)
    for(i = 1; i < _instance->N; i++)
      for(j = 0; j < i; j++) {
	_instance->Qs[m][i][j] = _instance->Qs[m][i][j] + _instance->Qs[m][j][i];
	_instance->Qs[m][j][i] = 0;
      }
}

/***********************************************
 * Free the variable of a MUBQP problem instance
 *
 * @param _instance the structure of the instance
 *
 ***********************************************/
void freeInstance(struct Instance * _instance) {
  if (_instance != NULL) {
    unsigned int n, i;
    if (_instance->Qs != NULL) {
      for(n = 0; n < _instance->M; n++) {
	for(i = 0; i < _instance->N; i++) 
	  free(_instance->Qs[n][i]);
	free(_instance->Qs[n]);
      }

      free(_instance->Qs);
      _instance->Qs = NULL;
    }
  }
}

/***********************************************
 * Compute the fitness function
 *
 * @param _instance the structure of the instance
 * @param _solution the solution to evaluate
 * @param _objVec   the objective vector of the corresponding solution
 ***********************************************/
void eval(struct Instance * _instance, int * _solution, int * _objVec) {
  unsigned int n;
  for(n = 0; n < _instance->M; n++) 
    _objVec[n] = evalUBQP(_instance, n, _solution);
} 


/***********************************************
 *
 * Fitness function of a single-objective UBQP problem
 *
 * @param _instance the structure of the instance
 * @param _numObj the objective fonction to consider
 * @param _sol the solution to evaluate
 *
 ***********************************************/
int evalUBQP(struct Instance * _instance, unsigned _numObj, int * _sol) {
  int fit = 0;
  unsigned int i, j;

  for(i = 0; i < _instance->N; i++)
    if (_sol[i] == 1) 
      for(j = 0; j <= i; j++)
	if (_sol[j] == 1) 
	  fit += _instance->Qs[_numObj][i][j];
	
  return fit ;
}

/***********************************************
 *
 * print the MUBQP instance
 *
 * @param _instance the structure of the instance
 *
 ***********************************************/
void printInstance(struct Instance * _instance) {
  // parameters 
  printf("%lf %d %d %f\n", _instance->rho, _instance->M, _instance->N, _instance->d);

  // matrices
  int i, k, n;

  for(i = 0; i < _instance->N; i++) {
    for(k = 0; k < _instance->N; k++) {
      for(n = 0; n < _instance->M; n++) {
	printf("%d ", _instance->Qs[n][i][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
}




