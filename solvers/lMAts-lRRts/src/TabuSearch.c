/*
 * TabuSearch.c
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

#include "TabuSearch.h"

Fitness TabuSearch(Gene *seq, int N, Random rnd, Timer *timer, unsigned long long int* NFEs) {
    #include "EfficientFlip.c"

    Gene currentSeq[N];

    int bestMovesBits[N];
    int numBestMoves;

    Fitness thisMoveFitness, bestMoveFitness;

    // Maximum number of iterations
    int maxIter = N/2 + RandomNextIntUntil(rnd,N);

    // Range of tabu iterations
    int MinTabu = maxIter/10;
    int MaxTabu = maxIter/50;

    // Initialize tabu table
    int tabuList[N];
    memset(tabuList,0,N*sizeof(*tabuList));

    Fitness initialFitness = initTables(seq);
    (*NFEs) ++; // Borko

    memcpy(currentSeq, seq, N*sizeof(Gene));
    Fitness currentFitness = initialFitness;

    Fitness bestFitness = initialFitness;

    for (int iter=1; iter<=maxIter; iter++) {

        bestMoveFitness = INT_MAX; // infinite
        numBestMoves = 0;

        // Locate best moves
        for(int bit=0; bit<N; bit++) {

            // Fitness if bit is mutated
            Fitness newFitness = fitnessIfMutated(bit);
            (*NFEs)++; // Borko

            // Movement allowed if not tabu or if aspiration criterion holds
            if( (tabuList[bit]<iter) || ((newFitness<bestFitness)) )  {
                if(newFitness<bestMoveFitness)  {
                    bestMovesBits[0] = bit;
                    numBestMoves = 1;
                    bestMoveFitness = newFitness;
                } else if (newFitness==bestMoveFitness) {
                    bestMovesBits[numBestMoves] = bit;
                    numBestMoves++;
                }
            }
        }

        // All movements may be in tabu state
        if(numBestMoves>0) {
            // Choose randomnly one best movement and accept it
            int bestMoveBit = bestMovesBits[RandomNextIntUntil(rnd,numBestMoves)];

            currentSeq[bestMoveBit] = -currentSeq[bestMoveBit];
            currentFitness = bestMoveFitness;
            updateTables(bestMoveBit);

            // make it tabu
            tabuList[bestMoveBit] = iter + MinTabu + (MaxTabu > 0 ? RandomNextIntUntil(rnd,MaxTabu) : 0);

            // Update best solution found, if needed
            if(bestMoveFitness<bestFitness) {
                memcpy(seq,currentSeq,N*sizeof(Gene));
                bestFitness = bestMoveFitness;
            }
        }
    }

    return bestFitness;
}
