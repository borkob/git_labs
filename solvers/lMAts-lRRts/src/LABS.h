/*
 * LABS.h
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/

#ifndef LABS_H_
#define LABS_H_

#include <stdio.h>

#include "types.h"
#include "random.h"


Fitness autocorrelation(Gene *seq, int N);

double meritFactor(Fitness f, int N);

double seqMeritFactor(Gene *seq, int N);

void LABSRandomSeq(Gene *seq, int N, Random rnd);

void readLABSOpts(char* fileName);

Fitness knownOptimum(int N);

void setOpts(const int size, const int target);

void fprintSeq(FILE *fp, Gene *seq, int N);

#endif /* LABS_H_ */
