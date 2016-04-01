/*
 * EfficientFlip.h
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/

int tableC[N][N];
int vectorC[N];	

Fitness initTables(Gene *seq) {
	for(int i=0; i<N-1; i++)
		for(int j=0; j<N-1-i; j++)
			tableC[j][i] = seq[j]*seq[i+1+j];

	Fitness initialFitness = 0;
	for(int i=0; i<N-1; i++) {
		vectorC[i] = 0;
		for(int j=0; j<N-1-i; j++)
			vectorC[i] += tableC[j][i];
		initialFitness += (vectorC[i]*vectorC[i]);
	}
	return initialFitness;
}


Fitness fitnessIfMutated(int bit) {
	int *ptrTableC_bit_i;
	int *ptrTableC_bit_i_plus1_i;
	int *ptrC, C;

	ptrTableC_bit_i_plus1_i = &tableC[bit-(0+1)][0];
	ptrTableC_bit_i = &tableC[bit][0];

	Fitness newFitness = 0;
	ptrC = &vectorC[0];

	for(int i=0,topePrimInd=N-1;i<N-1;i++,topePrimInd--,ptrTableC_bit_i++,ptrC++) {
		C = *ptrC;

		// Est?dentro del rango de los primeros ?dices?
		if(bit<topePrimInd) 
			C -=  ((*ptrTableC_bit_i) << 1);

		// Est?dentro del rango de los segundos ?dices?
		if(bit>i) 
			C -=  ((*ptrTableC_bit_i_plus1_i) << 1);

		ptrTableC_bit_i_plus1_i -= (N-1);

		newFitness += (C*C);
	}
	return newFitness;
}

void updateTables(int bit) {
	int *ptrTableC_bit_i;
	int *ptrTableC_bit_i_plus1_i;
	int *ptrC;	

	ptrTableC_bit_i = &tableC[bit][0];
	ptrTableC_bit_i_plus1_i = &tableC[bit-(0+1)][0];
	ptrC = &vectorC[0];


	for(int i=0,topePrimInd=N-1;i<N-1;i++,topePrimInd--,ptrTableC_bit_i++,ptrC++) {

		// Est?dentro del rango de los primeros ?dices?
		if(bit<topePrimInd) {
			(*ptrTableC_bit_i) *= -1;
			(*ptrC) += ((*ptrTableC_bit_i) << 1);
		}

		// Est?dentro del rango de los segundos ?dices?
		if(bit>i) {
			(*ptrTableC_bit_i_plus1_i) *= -1;
			(*ptrC) += ((*ptrTableC_bit_i_plus1_i) << 1) ;
		}
		ptrTableC_bit_i_plus1_i -= (N-1);
	}
}	

