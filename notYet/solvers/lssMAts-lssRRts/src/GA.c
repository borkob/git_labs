/*
 * GA.c
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
 */

#include "GA.h"
#include "LABS.h"
#include "random.h"
#include "dynamicMem.h"
#include "timer.h"
#include "TabuSearch.h"
#include "LABS.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>

#define _isBetterThan_(x,y) ((x) < (y))
#define _hasBetterFitnessThan_(x,y) (_isBetterThan_((x)->fitness, (y)->fitness))


// Borko (start)
static char solverName[255];
static char coordInit[512];
int size;
double maxExecTime;
int seed;
int threads;
 FILE * file;

void setCoordInit(const char aCoordInit[]){
    memcpy(coordInit,aCoordInit,strlen(aCoordInit)*sizeof(char));
}
// Borko (end)

double individualMeritFactor(IndividualP indP, GAP gaP) {
    return meritFactor(indP->fitness, gaP->chromosomeLgth) ;
}

void initIndividual(IndividualP indP, GAP gaP) {
    LABSRandomSeq(indP->chromosome, gaP->chromosomeLgth, gaP->rnd);
}

void initIndividual2(IndividualP indP, GAP gaP, char init[]) {
    const int N = gaP->chromosomeLgth;
    int NPlus1Div2 = (N+1) / 2;
    for(int i=0; i<NPlus1Div2; i++) {
        if(init[i] == '1') indP->chromosome[i] = -1;
        else if(init[i] == '0') indP->chromosome[i] = 1;
        else {
            fprintf(stderr,"Error: wrong init coordinate:%s",init);
            exit(1);
        }
    }
}

void setCoordInit2(IndividualP indP, GAP gaP, char init[]){
    const int N = gaP->chromosomeLgth;
    int NPlus1Div2 = (N+1) / 2;
    for(int i=0; i<NPlus1Div2; i++) {
        if(indP->chromosome[i] == -1) init[i] = '1';
        else if(indP->chromosome[i] == 1) init[i] = '0';
        else {
            fprintf(stderr,"Error: wrong init seq:%s",init);
            exit(1);
        }
    }
    init[NPlus1Div2] = '\0';
}

Fitness evalIndividual(IndividualP indP, GAP gaP) {
    return autocorrelation(indP->chromosome, gaP->chromosomeLgth);
}

void mutateIndividual(IndividualP indP, GAP gaP) {
    for(int i=0; i<gaP->chromosomeLgthPlus1Div2; i++) {
        if(withProb(gaP->rnd, gaP->mutationProb))
            indP->chromosome[i] = - indP->chromosome[i];
    }
}

IndividualP newIndividual(GAP gaP) {
    IndividualP indP = allocate(sizeof(Individual));
    indP->chromosome = allocate(gaP->chromosomeLgthPlus1Div2*sizeof(Gene));
    return indP;
}

void freeIndividual(IndividualP indP) {
    deallocate(indP->chromosome);
    deallocate(indP);
}

void copyIndividualFromTo(IndividualP from, IndividualP to, GAP gaP) {
    to->fitness = from->fitness;
    memcpy(to->chromosome,from->chromosome,gaP->chromosomeLgthPlus1Div2*sizeof(Gene));
}

void crossover(IndividualP parent1P, IndividualP parent2P, IndividualP childP, GAP gaP) {
    for(int i=0; i<gaP->chromosomeLgthPlus1Div2; i++)
        childP->chromosome[i] = (RandomNextBool(gaP->rnd) ? parent1P : parent2P)->chromosome[i];
}

bool isOptimal(IndividualP indP, GAP gaP) {
    gaP->tr = indP->fitness == knownOptimum(gaP->chromosomeLgth);
    if(indP->fitness < knownOptimum(gaP->chromosomeLgth)) gaP->tr = 2;
    return indP->fitness <= knownOptimum(gaP->chromosomeLgth);
}

GAP newGA(int id, int populationSz, int chromosomeLgth) {

    GAP gaP = allocate(sizeof(GA));
    gaP->GAId = id;
    gaP->populationSz = populationSz;
    gaP->chromosomeLgth = chromosomeLgth;
    gaP->chromosomeLgthPlus1Div2 = (chromosomeLgth+1) / 2;

    gaP->population = allocate(populationSz*sizeof(IndividualP));
    for(int i=0; i<populationSz; i++) {
        gaP->population[i] = newIndividual(gaP);
    }

    gaP->child = newIndividual(gaP);
    gaP->aux = newIndividual(gaP);

    gaP->best = newIndividual(gaP);
    gaP->best->fitness = INT_MAX;

    return gaP;
}

void freeGA(GAP gaP) {
    for(int i=0; i<gaP->populationSz; i++) {
        freeIndividual(gaP->population[i]);
    }
    freeIndividual(gaP->child);
    freeIndividual(gaP->best);
    freeIndividual(gaP->aux);

    deallocate(gaP);
}

// Tournament selection. Returns index into population
int tournament(GAP gaP, int tournamentSz) {
    int bestIndividual = RandomNextIntUntil(gaP->rnd,gaP->populationSz);
    for(int i=0; i<tournamentSz; i++) {
        int bestIndividual2 = RandomNextIntUntil(gaP->rnd,gaP->populationSz);
        if(_hasBetterFitnessThan_(gaP->population[bestIndividual2],gaP->population[bestIndividual]))
            bestIndividual = bestIndividual2;
    }
    return bestIndividual;
}

// selects an individual. Returns index into population
int randomIndividual(GAP gaP) {
    return RandomNextIntUntil(gaP->rnd, gaP->populationSz);
}

bool areEqIndividuals(GAP gaP, IndividualP indP1, IndividualP indP2) {
    if(indP1->fitness != indP2->fitness)
        return false;
    else
        return memcmp(indP1->chromosome, indP2->chromosome, gaP->chromosomeLgthPlus1Div2) == 0;
}

bool isAlreadyInPopulation(GAP gaP, IndividualP indP) {
    for(int i=0; i<gaP->populationSz; i++) {
        if(areEqIndividuals(gaP, gaP->population[i], indP))
            return true;
    }
    return false;
}

void updateBestIndividual(IndividualP ind, GAP gaP) {
    if(_hasBetterFitnessThan_(ind, gaP->best)) {

        copyIndividualFromTo(ind, gaP->best, gaP);

    double elapsedTime =  ((double)(clock() - gaP->startTime))/CLOCKS_PER_SEC;
    if(!getWalk()) // Borko (start)
        printf("# GA(%2d - %2d) iter: %12d  %4d  %.2f  %.1f(secs) %e(cntProbe) %e(speed) %d(targetReached) %d(isCensored)\n", // Borko
                gaP->GAId,
                gaP->seed,
                gaP->iter,
                gaP->best->fitness,
                individualMeritFactor(gaP->best, gaP),
                elapsedTime,
                (double)gaP->NFEs,
                (double)gaP->NFEs/elapsedTime,
                gaP->tr,
                gaP->ic
                ); // Borko (end)
        if(isOptimal(gaP->best, gaP)){
            gaP->stop = true;
        }
    }
}

// Initializes all individuals in population
void initPopulation(GAP gaP) {
    int bestIndividualIdx = 0;

    for(int i=0; i<gaP->populationSz; i++) {
        initIndividual(gaP->population[i], gaP);
        gaP->population[i]->fitness = evalIndividual(gaP->population[i], gaP);

        if(_hasBetterFitnessThan_(gaP->population[i], gaP->population[bestIndividualIdx]))
            bestIndividualIdx = i;
    }
    gaP->NFEs = gaP->populationSz; // Borko
    updateBestIndividual(gaP->population[bestIndividualIdx], gaP);
}

void print(){ // Borko (start)
    char filename[500];
    if(size < 10)
        sprintf(filename,"n=00%d-B-0_%s_-%s-walk0.txt",size,solverName,coordInit);
    else if(size < 100)
        sprintf(filename,"n=0%d-B-0_%s_-%s-walk0.txt",size,solverName,coordInit);
    else
        sprintf(filename,"n=%d-B-0_%s_-%s-walk0.txt",size,solverName,coordInit);


    #ifndef NDEBUG
    file = fopen(filename,"w");
    if(getWalk()) setFile(file);
    #endif

    printf("# %21s %s\n"," ",filename); // Borko (start)
    printf("# -->>>------------------------------------------------------------------\n");
    printf("# %-21s labs\n","labs");
    printf("# %-21s %d\n","instanceDef",size);
    printf("# %-21s %d\n","laevus",0);
    printf("# %-21s %d\n","isSkew",1);
    printf("# %-21s %s\n","solverName",solverName);
    printf("# %-21s %s\n","coordinateType","B");
    printf("# %-21s %d\n","nDim",(size+1)/2);
    printf("# %-21s %d\n","valueTarget",knownOptimum(size));
    printf("# %-21s %0.2f\n","meritTarget",meritFactor(knownOptimum(size), size));
    printf("# %-21s %0.2f\n","runtimeLmt",maxExecTime);
    printf("# LABS %d (opt=%d,%.2f), seed=%d, maxExecTime=%.1f(secs), threads=%d\n", size, knownOptimum(size), meritFactor(knownOptimum(size), size),seed, maxExecTime, threads);
    printf("# -->>>------------------------------------------------------------------\n");

    #ifndef NDEBUG
    if(getWalk()){
        fprintf(file,"# %21s %s\n"," ",filename);
        fprintf(file,"# -->>>------------------------------------------------------------------\n");
        fprintf(file,"# %-21s labs\n","labs");
        fprintf(file,"# %-21s %d\n","instanceDef",size);
        fprintf(file,"# %-21s %d\n","laevus",0);
        fprintf(file,"# %-21s %d\n","isSkew",1);
        fprintf(file,"# %-21s %s\n","solverName",solverName);
        fprintf(file,"# %-21s %s\n","coordinateType","B");
        fprintf(file,"# %-21s %d\n","nDim",(size+1)/2);
        fprintf(file,"# %-21s %d\n","valueTarget",knownOptimum(size));
        fprintf(file,"# %-21s %0.2f\n","meritTarget",meritFactor(knownOptimum(size), size));
        fprintf(file,"# %-21s %0.2f\n","runtimeLmt",maxExecTime);
        fprintf(file,"# LABS %d (opt=%d,%.2f), seed=%d, maxExecTime=%.1f(secs), threads=%d\n", size, knownOptimum(size), meritFactor(knownOptimum(size), size),seed, maxExecTime, threads);
        fprintf(file,"# -->>>------------------------------------------------------------------\n");
        printf("%s\t%s\t%s\t%s\t%s\n","step","chain","coord","value","comment");
        fprintf(file,"%s\t%s\t%s\t%s\t%s\n","step","chain","coord","value","comment");
    }
    #endif
} // Borko (end)

void *RunGA(void *ptr) {
    GAP gaP = (GAP) ptr;

    gaP->iter = 0;
    initPopulation(gaP);

    unsigned int i = 0;
    if(isOptimal(gaP->best, gaP)){

        setCoordInit2(gaP->best, gaP, coordInit);
        print();
    }
    if(!isOptimal(gaP->best, gaP)){
      while(!(gaP->stop)) {

    #if defined lssRRts || defined lssRRts0 // Borko (start)
        #ifdef lssRRts
            initIndividual(gaP->child, gaP);
        #else
            copyIndividualFromTo(gaP->population[i], gaP->child, gaP);
            i ++;
            if(i >= gaP->populationSz) i = 0;
        #endif
    #else
        // reproduction
        if(withProb(gaP->rnd, gaP->crossoverProb)) {
            int parentIdx1 = tournament(gaP, gaP->tournamentSz);
            int parentIdx2 = tournament(gaP, gaP->tournamentSz);

            crossover(gaP->population[parentIdx1], gaP->population[parentIdx2], gaP->child, gaP);
        } else
            copyIndividualFromTo(gaP->population[randomIndividual(gaP)], gaP->child, gaP);

        // mutation
        mutateIndividual(gaP->child, gaP);
    #endif // Borko (end)
        //gaP->child->fitness = evalIndividual(gaP->child, gaP);
        // Local Search

        if(gaP->iter==0){
            if(coordInit[0] != '\n') initIndividual2(gaP->child, gaP,coordInit);
            else setCoordInit2(gaP->child, gaP, coordInit);
            print();
        }

        gaP->child->fitness = TabuSearch(gaP->child->chromosome, gaP->chromosomeLgth, gaP->rnd, gaP->startTime,&gaP->NFEs); // Borko &gaP->NFEs

        // replacement
        if(!isAlreadyInPopulation(gaP, gaP->child)) {
            int forDeathIdx = randomIndividual(gaP);
            IndividualP aux = gaP->child;
            gaP->child = gaP->population[forDeathIdx];
            gaP->population[forDeathIdx] = aux;
            //updateBestIndividual(gaP->child, gaP);                 // Borko
            updateBestIndividual(gaP->population[forDeathIdx], gaP); // Borko
        }
        (gaP->iter)++;

    double elapsedTime = ((double)(clock() - gaP->startTime))/CLOCKS_PER_SEC;
        if(elapsedTime > gaP->maxExecTime){
            gaP->stop = true;
            gaP->ic = true;
        }
    }
    }
    else{
    double elapsedTime = ((double)(clock() - gaP->startTime))/CLOCKS_PER_SEC;
    if(!getWalk()) // Borko (start)
            printf("GA(%2d - %2d) iter: %12d  %4d  %.2f  %.1f(secs) %e(cntProbe) %e(speed) %d(targetReached) %d(isCensored)\n", // Borko
                gaP->GAId,
                gaP->seed,
                gaP->iter,
                gaP->best->fitness,
                individualMeritFactor(gaP->best, gaP),
                elapsedTime,
                (double)gaP->NFEs,
                (double)gaP->NFEs/elapsedTime,
                gaP->tr,
                gaP->ic
                ); // Borko (end)
    }
}

int main(int argc, char *args[]) {

    setbuf(stdout, NULL);
    installBackTracer();

    if(argc < 5 || argc > 9) {
    #ifndef NDEBUG     // Borko (start)
        printf("usage: %s <binary seq length> <random seed> <max time (secs)> <no. of threads> <valueTarget> (<-trace> | <-walk> | <-coordInit coord>)\n\n",args[0]);
    #else
        printf("usage: %s <binary seq length> <random seed> <max time (secs)> <no. of threads> <valueTarget>\n\n",args[0]);
    #endif     // Borko (end)
        printf("note1: The last argument is optional. If no <valueTarget> is specified,\n");
        printf("\tthe \"best known value\", stored internally, will be accessed.\n");
        printf("note2: This program stops either\n");
        printf("\tif runtime   >= <max time (secs)> or\n");
        printf("\tif valueBest <= <valueTarget>\n\n");
        printf("Copyright 2012\n");
        printf("*  José E. Gallardo, Carlos Cotta, and Antonio J. Fernández\n");
        printf("*  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.\n");
        printf("*  Applied Soft Computing. 9(4): 1252-1262 (2009).\n");
        printf("*  Modified by Borko.\n");
        exit(1);
    }

    size = atoi(args[1]);
    seed = atoi(args[2]);
    maxExecTime = atoi(args[3]);
    threads = atoi(args[4]);
    file = NULL;

    if(size%2==0)  {
        printf("<seq length> should be odd for skew-symmetric sequences\n",args[0]);
        exit(1);
    }

    // Borko (start)
    #ifdef lssMAts
    memcpy(solverName,"lssMAts",8);
    #endif
    #ifdef lssRRts
    memcpy(solverName,"lssRRts",8);
    #endif
    #ifdef lssRRts0
    memcpy(solverName,"lssRRts0",9);
    #endif

    coordInit[0] = '\n';
    if(argc >= 6 && atoi(args[5]) != 0){
        setOpts(size,atoi(args[5]));
    }
    setValueTarget(knownOptimum(size));
    for(unsigned int i=6; i<argc; i++){
        if(!strcmp("-trace",args[i])) setTrace();
        else if(!strcmp("-walk",args[i])) setWalk();
        else if(!strcmp("-coordInit",args[i])){
            i++;
            setCoordInit(args[i]);
        }
    }

    // Create random generators
    Random rnds[threads];
    for(int i=0; i<threads; i++) {
        rnds[i] = newRandom(i, seed)    ;
    }

   /* TimerP timerP = newTimer();
    StartTimer(timerP);*/
   clock_t startTime = clock() - 0.001;

    GAP gasP[threads];
    for(int i=0; i<threads; i++) {
        #if defined lssRRts
        gasP[i] = newGA(i, 1, size);
        #else
        gasP[i] = newGA(i, 100, size);
        #endif
        gasP[i]->seed = seed;
        gasP[i]->rnd = rnds[i];
        gasP[i]->startTime = startTime;
        gasP[i]->crossoverProb = 0.9;
        gasP[i]->mutationProb = 1.0/(double)gasP[i]->chromosomeLgthPlus1Div2;
        gasP[i]->tournamentSz = 2; // binary tournament
        gasP[i]->maxExecTime = maxExecTime;
        gasP[i]->stop = false;
        gasP[i]->ic = false;
    }

    int threadId[threads];
    pthread_t workerThreads[threads];

    for(int i=0; i<threads; i++) {
        threadId[i] = pthread_create(&(workerThreads[i]), NULL, RunGA, (void*) gasP[i]);
    }

    for(int i=0; i<threads; i++) {
        pthread_join(workerThreads[i], NULL);
    }

    /* StopTimer(timerP);*/
    for(int i=0; i<threads; i++) {
        int tr = 0;
        if(gasP[i]->best->fitness == knownOptimum(size)) tr = 1;
        if(gasP[i]->best->fitness < knownOptimum(size)) tr = 2;

    double elapsedTime =  ((double)(clock() - gasP[i]->startTime))/CLOCKS_PER_SEC;
        if(getWalk()){ // Borko (start)
            if(getChain() == 0) printf("%d\t%d\t",getStep()+1, getChain());
            else printf("%d\t%d\t",getStep()+1, getChain()-1);
            fprintSeq3(stdout, gasP[i]->best->chromosome, gasP[i]->chromosomeLgth);
            printf("# PIVOT (TARGET) **END**\n");
            if(getChain() == 0) fprintf(file,"%d\t%d\t",getStep()+1, getChain());
            else fprintf(file,"%d\t%d\t",getStep()+1, getChain()-2);
            fprintSeq3(file, gasP[i]->best->chromosome, gasP[i]->chromosomeLgth);
            fprintf(file,"# PIVOT (TARGET) **END**\n");
        }
        if(getTrace()){
            if(getChain() == 0) printf("walkLength=0 cntRestart=%d cntProbe=%llu pivot=",getChain(),gasP[i]->NFEs);
            else printf("walkLength=0 cntRestart=%d cntProbe=%llu pivot=",getChain()-1,gasP[i]->NFEs);
            fprintSeq2(stdout, gasP[i]->best->chromosome, gasP[i]->chromosomeLgth);
            printf("\t PIVOT (TARGET) **END**\n");
        }
        if(!getWalk()){
            printf("# GA(%2d - %2d) iter: %12d  %4d  %.2f  %.1f(secs) %e(cntProbe) %e(speed) %d(targetReached) %d(isCensored)\n# BEST:", // Borko
                gasP[i]->GAId,
                gasP[i]->seed,
                gasP[i]->iter,
                gasP[i]->best->fitness,
                individualMeritFactor(gasP[i]->best, gasP[i]),
                elapsedTime,
                (double)gasP[i]->NFEs,
                (double)gasP[i]->NFEs/elapsedTime,
                gasP[i]->tr,
                gasP[i]->ic
                );
        	fprintSeq(stdout, gasP[i]->best->chromosome, gasP[i]->chromosomeLgth);

    		printf("\n%-21s %d\n","instanceDef",size); // Borko - results
		printf("%-21s %d\n","nDim",(size+1)/2);
		printf("%-21s %s\n","progName",solverName);
		printf("%-21s %s\n","progVesrsion","0.1");
    		printf("%-21s %0.4f\n","meritTarget",meritFactor(knownOptimum(size), size));
    		printf("%-21s %d\n","valueTarget",knownOptimum(size));
    		printf("%-21s %0.4f\n","meritBest",individualMeritFactor(gasP[i]->best, gasP[i]));
    		printf("%-21s %d\n","valueBest", gasP[i]->best->fitness);
    		printf("%-21s %d\n","targetReached", gasP[i]->tr);
    		printf("%-21s %d\n","isCensored",gasP[i]->ic);
    		printf("%-21s %d\n","runtimeLmt",(unsigned int)maxExecTime);
    		printf("%-21s %.2f\n","runtime",elapsedTime);
    		printf("%-21s %llu\n","cntProbe",(unsigned long long int)gasP[i]->NFEs);
    		printf("%-21s %d\n","speed",(unsigned int)((double)gasP[i]->NFEs/elapsedTime));
    		printf("%-21s %d\n","seedFirst",seed);
    		printf("%-21s %s\n","coordInit",coordInit);
    		printf("%-21s ","coordBest");
        	fprintSeq4(stdout, gasP[i]->best->chromosome, gasP[i]->chromosomeLgth);
    		printf("%-21s ","coordBestFull");
        	fprintSeq5(stdout, gasP[i]->best->chromosome, gasP[i]->chromosomeLgth);
        }
    }
    if(!getWalk()) printf("\n"); // Borko (end)

    for(int i=0; i<threads; i++) {
        freeGA(gasP[i]);
        freeRandom(&rnds[i]);
    }
    // freeTimer(&timerP);

    #ifndef NDEBUG // Borko (start)
    if(getWalk() && file != NULL) fclose(file);
    #endif // Borko (end)
}

#if 0

void scatterCrossover(IndividualP parent1P, IndividualP parent2P, IndividualP childP, GAP gaP) {
    int N = gaP->chromosomeLgth;

    #include "EfficientFlip.c"

    Gene currentSeq[N], LSSeq[N];

    memcpy(currentSeq, parent1P->chromosome, N*sizeof(Gene));

    Fitness initialFitness = initTables(currentSeq);

    Fitness bestFitness = INT_MAX;

    for(int bit = 0; bit<N; bit++) {
        if(currentSeq[bit] != parent2P->chromosome[bit]) {

            currentSeq[bit] = -currentSeq[bit];
            Fitness newFitness = fitnessIfMutated(bit);
            updateTables(bit);

            if(newFitness < bestFitness) {
                bestFitness = newFitness;
                memcpy(childP->chromosome, currentSeq, N*sizeof(Gene));
            }
        }

    }
    childP->fitness = bestFitness;
}

void scatterCrossover2(IndividualP parent1P, IndividualP parent2P, IndividualP childP, GAP gaP) {
    int N = gaP->chromosomeLgth;

    Gene currentSeq[N], LSSeq[N];

    memcpy(currentSeq, parent1P->chromosome, N*sizeof(Gene));

    Fitness bestFitness = INT_MAX;

    for(int bit = 0; bit<N; bit++) {
        if(currentSeq[bit] != parent2P->chromosome[bit]) {
        if(withProb(gaP->rnd,0.5)) {
            currentSeq[bit] = -currentSeq[bit];
            memcpy(LSSeq, currentSeq, N*sizeof(Gene));
            Fitness LSFitness = TabuSearch(LSSeq, N, gaP->rnd, NULL);
            if(LSFitness < bestFitness) {
                bestFitness = LSFitness;
                memcpy(childP->chromosome, LSSeq, N*sizeof(Gene));
            }
        }
        }
    }
    childP->fitness = bestFitness;
}
#endif
