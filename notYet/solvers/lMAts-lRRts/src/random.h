/*
 * random.h
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/

#ifndef RANDOM_H_
#define RANDOM_H_

#include <stdbool.h>
#include <stdint.h>

#include "libs/dcmt0.6.1/include/dc.h"

typedef mt_struct *Random;
typedef uint32_t uint32;

Random newRandom(int id, uint32 seed);


void freeRandom(Random *rnd);

// random integer between 0 and UINT32_MAX inclusive
uint32 RandomNextInt(Random rnd);

// random double in [0.0,1.0) (0.0 inclusive, 1.0 exclusive)
double RandomNextDouble(Random rnd);

int RandomNextIntUntil(Random rnd, uint32 n);

bool RandomNextBool(Random rnd);

#define withProb(rnd,p)	(RandomNextDouble((rnd))<(p))

#endif /* RANDOM_H_ */
