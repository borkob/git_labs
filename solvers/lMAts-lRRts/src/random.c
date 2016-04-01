/*
 * random.c
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/

#include<stdio.h>
#include<stdlib.h>

#include "random.h"

// Use different ids for different generators 
// Use different seeds for starting at different positions in sequence
Random newRandom(int id, uint32 seed) {
    Random mt = get_mt_parameter_id_st(32,521,id,seed);  
    if (mt == NULL) {
		fprintf(stderr,"\nnewRandom: get_mt_parameter_id_st failed\n");
		exit(1);	
    }
	sgenrand_mt(seed, mt);
	return mt;
}	


void freeRandom(Random *rnd) {
	free_mt_struct(*rnd);
	*rnd = NULL;
}

// random integer between 0 and UINT32_MAX inclusive
uint32 RandomNextInt(Random rnd) {
	return genrand_mt(rnd);
}

// random double in [0.0,1.0) (0.0 inclusive, 1.0 exclusive)
double RandomNextDouble(Random rnd) {
	return ((double)genrand_mt(rnd) / ((double)(UINT32_MAX) + 1.0));
}

int RandomNextIntUntil(Random rnd, uint32 n) {
	return genrand_mt(rnd) % n;
}

bool RandomNextBool(Random rnd) {
	return genrand_mt(rnd) % 2;
}
