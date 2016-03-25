/*
 * types.h
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/

#ifndef TYPES_H_
#define TYPES_H_
#include <stdbool.h>

typedef void *voidP;

typedef signed char SignedByte;
typedef unsigned char Byte;
typedef Byte *ByteP;


typedef SignedByte Gene;

typedef Gene *Chromosome;

typedef int Fitness;


#define MIN(x,y) ((x)<=(y) ? (x) : (y))
#define MAX(x,y) ((x)>=(y) ? (x) : (y))
#define SQR(x) 	 ((x)*(x))
#define ABS(x)   ((x)>=0 ? (x) : (-(x)))

#endif /* TYPES_H_ */
