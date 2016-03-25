/*
 * dynamicMem.h
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/

#ifndef DYNAMICMEM_H_
#define DYNAMICMEM_H_

#include <stddef.h>
#include "types.h"

#define allocate(size) allocateFileLine((size), __FILE__, __LINE__)
#define deallocate(ptr) deallocateFileLine(((voidP)&(ptr)), __FILE__, __LINE__)
#define deallocate2(ptr) deallocateFileLine(((voidP)(ptr)), __FILE__, __LINE__)

#define notNull(ptr) notNullFileLine(((voidP)(ptr)), __FILE__, __LINE__)

void *allocateFileLine(size_t size, char *file, int line);
void deallocateFileLine(void **ptr, char *file, int line);
void notNullFileLine(void *ptr, char *file, int line);
void backtraceHandler(int sig);
void installBackTracer();

#endif /* DYNAMICMEM_H_ */
