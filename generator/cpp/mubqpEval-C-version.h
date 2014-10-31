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

#ifndef __mubqpEvalC_H
#define __mubqpEvalC_H

#include <stdlib.h>
#include <stdio.h>

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

/*
 * Structure of an mUBQP problem instance
 */
struct Instance {
  // correlation between contributions
  double rho;

  // number of objective functions
  unsigned M;

  // size of the bit string
  unsigned N;

  // density of the matrices (frequency of non-zero numbers)
  double d;

  // the M matrices of contributions
  int *** Qs;

};


/***********************************************
 *
 * Load the file of a mubqp problem instance
 *
 * @param fileName file name instance of the mubqp problem
 * @param _instance the structure of the instance
 *
 ***********************************************/
void loadInstance(const char * , struct Instance * ) ;

/***********************************************
 *
 * Initialization of the differents tables and epistasis links
 *
 * @param _instance the structure of the instance
 *
 ***********************************************/
void initInstance(struct Instance * ) ;

/***********************************************
 *
 * Load the matrices
 *
 * @param _file open file of the instance
 * @param _instance the structure of the instance
 *
 ***********************************************/
void loadQs(FILE* , struct Instance * ) ;

/***********************************************
 * Free the variable of an mUBQP problem instance
 *
 * @param _instance the structure of the instance
 *
 ***********************************************/
void freeInstance(struct Instance * ) ;

/***********************************************
 * Compute the fitness function
 *
 * @param _instance the structure of the instance
 * @param _solution the solution to evaluate
 * @param _objVec   the objective vector of the corresponding solution
 ***********************************************/
void eval(struct Instance * , int * , int * ) ;

/***********************************************
 *
 * Fitness function of a single-objective UBQP problem
 *
 * @param _instance the structure of the instance
 * @param _numObj the objective fonction to consider
 * @param _sol the solution to evaluate
 *
 ***********************************************/
int evalUBQP(struct Instance * , unsigned , int * ) ;
 
/***********************************************
 *
 * print the mUBQP problem instance
 *
 * @param _instance the structure of the instance
 *
 ***********************************************/
void printInstance(struct Instance * _instance) ;

#endif

