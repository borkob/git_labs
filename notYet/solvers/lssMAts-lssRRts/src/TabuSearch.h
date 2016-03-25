/*
 * TabuSearch.h
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
 */

#ifndef TABU_SEARCH_H_
#define TABU_SEARCH_H_

#include "types.h"
#include "random.h"
#include "time.h"

Fitness TabuSearch(Gene *initialSeq, int N, Random rnd, clock_t startTime, unsigned long long int * NFEs);
void setTrace();
void setWalk();
void setFile(FILE* aFile);
bool getTrace();
bool getWalk();
int getStep();
int getChain();
void setValueTarget(const int valueTarget);
#endif /* TABU_SEARCH_H_ */
