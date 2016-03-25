/* genmtrand.c */

/* Copyright (C) 2001-2009 Makoto Matsumoto and Takuji Nishimura.  */
/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */ 
/* 02111-1307  USA                                                 */

#include <stdio.h>
#include <stdlib.h>
#include "dci.h"

#define SHIFT1 18

void sgenrand_mt(uint32_t seed, mt_struct *mts) 
{
    int i;

    for (i=0; i<mts->nn; i++) {
	mts->state[i] = seed;
        seed = (UINT32_C(1812433253) * (seed  ^ (seed >> 30))) + i + 1;
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
    }
    mts->i = mts->nn;

    for (i=0; i<mts->nn; i++) 
	mts->state[i] &= mts->wmask;
}

uint32_t genrand_mt(mt_struct *mts) 
{
    uint32_t *st, uuu, lll, aa, x;
    int k,n,m,lim;

    if ( mts->i >= mts->nn ) {
	n = mts->nn; m = mts->mm;
	aa = mts->aaa;
	st = mts->state;
	uuu = mts->umask; lll = mts->lmask;

	lim = n - m;
	for (k=0; k<lim; k++) {
	    x = (st[k]&uuu)|(st[k+1]&lll);
	    st[k] = st[k+m] ^ (x>>1) ^ (x&1U ? aa : 0U);
	}
	lim = n - 1;
	for (; k<lim; k++) {
	    x = (st[k]&uuu)|(st[k+1]&lll);
	    st[k] = st[k+m-n] ^ (x>>1) ^ (x&1U ? aa : 0U);
	}
	x = (st[n-1]&uuu)|(st[0]&lll);
	st[n-1] = st[m-1] ^ (x>>1) ^ (x&1U ? aa : 0U);
	mts->i=0;
    }
		
    x = mts->state[mts->i];
    mts->i += 1;
    x ^= x >> mts->shift0;
    x ^= (x << mts->shiftB) & mts->maskB;
    x ^= (x << mts->shiftC) & mts->maskC;
    x ^= x >> mts->shift1;

    return x;
}

