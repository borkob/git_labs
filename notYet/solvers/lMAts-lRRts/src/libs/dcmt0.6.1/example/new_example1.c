/*
    This is a simple example for Dynamic Creator library.
*/

#include <stdio.h>
#include "dc.h"

int main(void)
{
    int i,j;
    mt_struct *mts;

    /* This trys to find a small Mersenne Twister with period 2^521-1. */
    mts = get_mt_parameter_st(32,521,4172);
    if (mts == NULL) {
        printf ("error\n");
    }
    else {
        sgenrand_mt(3241, mts);
        for (i=0; i<100; i++) {
            for (j=0; j<5; j++)
                printf("%8"PRIx32" ", genrand_mt(mts));
            printf("\n");
        }
        free_mt_struct(mts);
    }

    return 0;
}

