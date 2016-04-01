/*
   This is a third simple example for Dynamic Creator library.
*/

#include <stdio.h>
#include "dc.h"

int main(void)
{
    mt_struct **mtss;
    int count;
    int seed[] = {1234, 4567, 8901};
    int i, j;

    init_dc(4172);

    /* This trys to find three independent small Mersenne Twisters
       with period 2^521-1. */
    printf ("get MT parameters.\n");
    /* start_id = 3, max_id = 5 */
    mtss = get_mt_parameters(32,521,2,&count);
    if (mtss == NULL) {
	printf ("error\n");
	return 1;
    }
    for (i = 0; i < count; i++) {
	sgenrand_mt(seed[i], mtss[i]);
    }
    for (i=0; i<10; i++) {
	for (j = 0; j < count; j++) {
	    printf("%8"PRIx32" ", genrand_mt(mtss[j]));
	}
	printf("\n");
    }
    free_mt_struct_array(mtss, count);
    return 0;
}

