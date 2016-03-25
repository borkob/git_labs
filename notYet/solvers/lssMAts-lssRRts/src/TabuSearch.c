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
#include "LABS.h"

bool trace = false;
bool walk = false;
int step = 0;
int chain = 0;
FILE * file;
int valueTarget = 0;

void setValueTarget(const int aValueTarget){ valueTarget = aValueTarget; }

void setTrace(){ trace = true; }
void setWalk(){ walk = true; }
void setFile(FILE* aFile){ file = aFile; }
bool getTrace() { return trace; }
bool getWalk() { return walk; }
int getStep() { return step; }
int getChain() { return chain; }

Fitness TabuSearch(Gene *seq, int N, Random rnd, clock_t startTime, unsigned long long int* NFEs) { // Borko (unsigned int* NFEs)
    #include "EfficientFlip.c"

    Gene currentSeq[NPlus1Div2];

    int bestMovesBits[NPlus1Div2];
    int numBestMoves;

    Fitness thisMoveFitness, bestMoveFitness;

    // Maximum number of iterations
    int maxIter = (N/2 + RandomNextIntUntil(rnd,N));

    // Range of tabu iterations
    int MinTabu = maxIter/10;
    int MaxTabu = maxIter/50;

    // Initialize tabu table
    int tabuList[NPlus1Div2];
    memset(tabuList, 0, NPlus1Div2*sizeof(*tabuList));

    Fitness initialFitness = initTables(seq);

    (*NFEs) ++;
    #ifndef NDEBUG // Borko (start)
    if(trace){
        printf("walkLength=0 cntRestart=%d cntProbe=%llu pivotInit=",chain,*NFEs);
        fprintSeq2(stdout, seq, N);
        printf("\n");
    }
    if(walk){
        printf("%d\t%d\t",step,chain);
        fprintSeq3(stdout, seq, N);
        if(initialFitness <= valueTarget) printf("# targetReached 1\n");
        else printf("#\n");
        fprintf(file,"%d\t%d\t",step,chain);
        fprintSeq3(file, seq, N);
        fprintf(file,"# Init\n");
        if(initialFitness <= valueTarget) fprintf(file,"# targetReached 1\n");
        else fprintf(file,"#\n");
    }
    #endif     // Borko (end)

    memcpy(currentSeq, seq, NPlus1Div2*sizeof(Gene));
    Fitness currentFitness = initialFitness;

    Fitness bestFitness = initialFitness;

    for (int iter=1; iter<=maxIter; iter++) {

        bestMoveFitness = INT_MAX; // infinite
        numBestMoves = 0;

        // Locate best moves
        for(int bit=0; bit<NPlus1Div2; bit++) {

            // Fitness if bit is mutated
            Fitness newFitness = fitnessIfMutated(bit);

            (*NFEs)++; // Borko (start)
            #ifndef NDEBUG
            if(trace){
                printf("\tprobe=%llu ",*NFEs);
                printf("flipped bit = %d, ",(bit+1));
                currentSeq[bit] = -currentSeq[bit];
                fprintSeq2(stdout, currentSeq, N);
                printf("\n");
                currentSeq[bit] = -currentSeq[bit];
            }
            #endif // Borko (end)

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

        #ifndef NDEBUG // Borko (start)
        if(walk) step ++;
        #endif // Borko (end)

        // All movements may be in tabu state
        if(numBestMoves>0) {
            // Choose randomnly one best movement and accept it
            int bestMoveBit = bestMovesBits[RandomNextIntUntil(rnd, numBestMoves)];

            currentSeq[bestMoveBit] = -currentSeq[bestMoveBit];
            currentFitness = bestMoveFitness;
            updateTables(bestMoveBit);

            #ifndef NDEBUG // Borko (start)
            if(trace){
                printf("walkLength=%d cntRestart=%d cntProbe=%llu ",iter,chain,*NFEs);
                printf("bit=%d pivot=",(bestMoveBit+1));
                fprintSeq2(stdout, currentSeq, N);
                printf("\n");
            }
            if(walk){
                printf("%d\t%d\t",step,chain);
                fprintSeq3(stdout, currentSeq, N);
                if(currentFitness <= valueTarget)
                    printf("# targetReached 1\n");
                else
                    printf("#\n");
                fprintf(file,"%d\t%d\t",step,chain);
                fprintSeq3(file, currentSeq, N);
                if(currentFitness <= valueTarget)
                    fprintf(file,"# targetReached 1\n");
                else
                    fprintf(file,"#\n");
            }
            #endif // Borko (end)

            // make it tabu
            tabuList[bestMoveBit] = iter + MinTabu + (MaxTabu > 0 ? RandomNextIntUntil(rnd,MaxTabu) : 0);

            // Update best solution found, if needed
            if(bestMoveFitness<bestFitness) {
                memcpy(seq, currentSeq, NPlus1Div2*sizeof(Gene));
                bestFitness = bestMoveFitness;
            }
        }
    }

    #ifndef NDEBUG // Borko (start)
    chain ++;
    #endif // Borko (end)
    return bestFitness;
}
