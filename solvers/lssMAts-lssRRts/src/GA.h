/*
 * GA.h
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
 */

#ifndef GA_H_
#define GA_H_

#include <pthread.h>

#include "types.h"
#include "random.h"
//#include "timer.h"
#include "time.h"


typedef struct {
    Fitness fitness;
    Chromosome chromosome;
} Individual;

typedef Individual *IndividualP;


typedef struct {
    int seed;
    Random rnd; // random generator for this GA

    int GAId;

    int chromosomeLgth; 		// Assumed to be same for all individuals in population
    int chromosomeLgthPlus1Div2;

    int populationSz;         	// number of individuals in population
    IndividualP *population; 	// Array with individuals
    IndividualP child;     		// Use for reproduction

    int iter;

    double mutationProb;
    double crossoverProb;

    int tournamentSz; // number of tournaments per selected individual

    IndividualP best; // best for whole execution
    IndividualP aux;

    double maxExecTime;

    // TimerP timerP;
    clock_t startTime;
    bool stop;
    bool ic;
    int tr;

    unsigned long long int NFEs; // Borko

} GA;

typedef GA *GAP;

#endif /* GA_H_ */
