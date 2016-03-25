/*
 * EfficientFlip.h
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
 */

int NDiv2 = N / 2;
int NPlus1Div2 = (N+1) / 2;
int NMinus1 = N-1;

#undef seqPos
#define seqPos(i) ((i) < NPlus1Div2 ? seq[i] : ( (NMinus1-(i))%2==(NDiv2%2) ? seq[NMinus1-(i)] : -seq[NMinus1-(i)] ) )

int tableC[N][N];
int vectorC[N];
int vectorCAux[N];

Fitness initTables(Gene *seq) {
    for(int i=0; i<NDiv2; i++){
        for(int j=0; j<NDiv2-i; j++){
            tableC[i][j] = 4*seq[j]*seqPos(j+(2*(i+1)));
        }
    }

    Fitness initialFitness = 0;
    for(int i=0; i<NDiv2; i++) {
        vectorC[i] = tableC[i][NDiv2-i-1]/4;
        for(int j=0; j<NDiv2-i-1; j++)
            vectorC[i] += (tableC[i][j]/2);
        initialFitness += SQR(vectorC[i]);
    }

    return initialFitness;
}

Fitness fitnessIfMutated(int bit) {
    int i, *ptrTableC, *ptrI, *ptrIEnd;
    Fitness newFitness;

    memcpy(vectorCAux,vectorC,NPlus1Div2*sizeof(*vectorC));

    // Mutar el bit bit-ésimo y recalcular fitness en f
    for(ptrI=&vectorCAux[0], ptrTableC=&tableC[0][bit-2*(0+1)], ptrIEnd=ptrI+bit/2; ptrI<ptrIEnd; ptrTableC+=(N-2),ptrI++) {
        *ptrI -= (*ptrTableC);
    }

    int simBit = NMinus1-bit;
    if(simBit != bit) {
        for(i=NDiv2-bit, ptrI=&vectorCAux[i], ptrIEnd=ptrI+(bit/2),ptrTableC=&tableC[i][bit-2*(0+1)]; ptrI<ptrIEnd; ptrTableC+=(N-2),ptrI++) {
            *ptrI -= (*ptrTableC);
        }
    }

    for(ptrI=&vectorCAux[0], ptrIEnd=ptrI+NDiv2-(bit+1), ptrTableC=&tableC[0][bit]; ptrI<ptrIEnd; ptrTableC+=N,ptrI++) {
        *ptrI -= (*ptrTableC);
    }

    newFitness = 0;
    for(ptrI=&vectorCAux[0],ptrIEnd=ptrI+NDiv2; ptrI<ptrIEnd; ptrI++) {
        newFitness += SQR(*ptrI);
    }

    return newFitness;
}

void updateTables(int bit) {
    int i, val, *ptrTableC, *ptrI, *ptrIEnd;

    for(ptrI=&vectorC[0], ptrTableC=&tableC[0][bit-2*(0+1)], ptrIEnd=ptrI+bit/2; ptrI<ptrIEnd; ptrTableC+=(N-2),ptrI++) {
        val = *ptrTableC;
        *ptrTableC = -val;
        *ptrI -= val;
    }

    int simBit = NMinus1-bit;
    if(simBit != bit) {
        for(i=NDiv2-bit, ptrI=&vectorC[i], ptrIEnd=ptrI+(bit/2),ptrTableC=&tableC[i][bit-2*(0+1)]; ptrI<ptrIEnd; ptrTableC+=(N-2),ptrI++) {
            val = *ptrTableC;
            *ptrTableC = -val;
            *ptrI -= val;
        }
    }

    for(ptrI=&vectorC[0], ptrIEnd=ptrI+NDiv2-(bit+1), ptrTableC=&tableC[0][bit]; ptrI<ptrIEnd; ptrTableC+=N,ptrI++) {
        val = *ptrTableC;
        *ptrTableC = -val;
        *ptrI -= val;
    }
}

