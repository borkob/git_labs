/*
 * LABS.c
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/


#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "LABS.h"
#include "random.h"


Fitness autocorrelation(Gene *seq, int N) {
    Fitness s, f = 0;

    for(int i=0; i<N-1; i++) {
        s = 0;
        for(int j=0; j<N-1-i; j ++)
            s +=  seq[j]*seq[i+1+j];
        f += (SQR(s));
    }
    return f;
}


double meritFactor(Fitness f, int N) {
    return SQR(N) / (2.0 * f);
}

double seqMeritFactor(Gene *seq, int N) {
    return meritFactor(autocorrelation(seq,N), N);
}

void LABSRandomSeq(Gene *seq, int N, Random rnd) {
    for(int i=0; i<N; i++) {
        seq[i] = RandomNextBool(rnd) ? 1 : (-1);
    }
}


void fprintSeq(FILE *fp, Gene *seq, int N) {
    for(int i=0; i<N; i++)
        fprintf(fp, "%c", seq[i] == 1 ? '+' : '-');
    fprintf(fp, " %d %.2f", autocorrelation(seq, N), seqMeritFactor(seq, N));
}

#define MAX_OPTS  304
Fitness labsOpts[1+MAX_OPTS] = {
       0,    0,    1,    1,    2,    2,    7,    3,    8,   12,   13,    5,
      10,    6,   19,   15,   24,   32,   25,   29,   26,   26,   39,   47,
      36,   36,   45,   37,   50,   62,   59,   67,   64,   64,   65,   73,
      82,   86,   87,   99,  108,  108,  101,  109,  122,  118,  131,  135,
     140,  136,  153,  153,  166,  170,  175,  171,  192,  188,  197,  205,
     218,  226,  235,  207,  208,  240,  257,  241,  250,  274,  295,  275,
     300,  308,  341,  329,  334,  358,  347,  339,  352,  372,  377,  377,
     430,  414,  439,  431,  448,  432,  453,  477,  498,  486,  499,  479,
     520,  536,  545,  577,  578,  578,  567,  555,  612,  620,  701,  677,
     702,  662,  723,  687,  788,  752,  817,  745,  814,  786,  847,  835,
     872,  844,  885,  893,  922,  846,  875,  887,  932,  920,  945,  913,
    1014, 1010, 1063, 1027, 1076, 1052, 1117, 1133, 1178, 1126, 1235, 1191,
    1248, 1208, 1273, 1265, 1298, 1218, 1279, 1275, 1388, 1340, 1429, 1437,
    1438, 1366, 1467, 1439, 1504, 1512, 1585, 1529, 1594, 1474, 1591, 1563,
    1620, 1532, 1693, 1677, 1674, 1606, 1699, 1719, 1780, 1808, 1929, 1897,
    1898, 1898, 2015, 1995, 2040, 2028, 2069, 1973, 1970, 1966, 2123, 2191,
    2304, 2272, 2337, 2281, 2374, 2218, 2343, 2275, 2412, 2460, 2541, 2421,
    2542, 2662, 2723, 2695, 2720, 2664, 2761, 2801, 2878, 2698, 2799, 2831,
    2968, 3036, 3173, 3189, 3322, 3206, 3211, 3215, 3392, 3416, 3569, 3409,
    3566, 3474, 3687, 3587, 3752, 3692, 3821, 3757, 3674, 3590, 3651, 3711,
    3948, 3992, 4073, 4073, 4150, 4098, 4223, 4291,    0, 4280, 4341, 4165,
    4386, 4382, 4587, 4463, 4584, 4472, 4705, 4705,    0, 4790, 4887, 4803,
    4892, 4948,    0, 5037, 5130, 4950, 5203, 5243,    0,    0,    0,    0,
       0,    0,    0,    0,    0, 5564,    0, 5697,    0, 5790,    0,    0,
       0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
       0,    0,    0, 6455, 6568
};

void readLABSOpts(char* fileName) {
    FILE* f;

    f = fopen(fileName, "rt");
    if(!f) {
        printf("\nCan't open %s file for reading", fileName);
        exit(1);
    }

    int size, opt;
    while(fscanf(f, "%d %d", &size, &opt) != EOF) {
        labsOpts[size] = opt;
    }

    fclose(f);
}

void setOpts(const int size, const int target){ // Borko
    labsOpts[size] = target;
}

Fitness knownOptimum(int N) {
    return ((N <= MAX_OPTS) ? labsOpts[N] : 0);
}


