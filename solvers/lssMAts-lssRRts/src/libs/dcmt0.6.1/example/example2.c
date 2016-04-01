/*
   This is a second simple example for Dynamic Creator library.
*/

#include <stdio.h>
#include "dc.h"

int main(void)
{
    int i;
    mt_struct *mts0, *mts1, *mts2;

    init_dc(4172);

    /* This trys to find three independent small Mersenne Twisters
       with period 2^521-1. */
	printf ("get first MT\n");
    mts0 = get_mt_parameter_id(32,521,0);  /* id=0 */
    if (mts0 == NULL) {
	printf ("error\n"); return 0; /* if did not find */
    }
	printf ("get second MT\n");
    mts1 = get_mt_parameter_id(32,521,1);  /* id=1 */
    if (mts1 == NULL) {
	printf ("error\n"); return 0;
    }
	printf ("get third MT\n");
    mts2 = get_mt_parameter_id(32,521,999);
	/* id may be any=16bit integers, e.g. id=999 */
    if (mts2 == NULL) {
	printf ("error\n"); return 0;
    }
	sgenrand_mt(1234, mts0); /* initialize mts0 with seed 1234 */
	sgenrand_mt(4567, mts1);
	sgenrand_mt(8901, mts2);
	for (i=0; i<10; i++) {
		printf("%8"PRIx32" ", genrand_mt(mts0));
		printf("%8"PRIx32" ", genrand_mt(mts1));
		printf("%8"PRIx32" ", genrand_mt(mts2));
	    printf("\n");
    /* print output of mts0, mts1, mts2, ten times */
	}
	free_mt_struct(mts0);
	free_mt_struct(mts1);
	free_mt_struct(mts2);
    return 0;
}

