/*
 * dynamicMem.c
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/

#include <stdlib.h>
#include <stdio.h>

#include "types.h"


#ifdef __RUN_IN_LINUX__
#include <signal.h>
#include <execinfo.h>
void backtraceHandler(int sig) {
  void *array[100];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 100);

  // print out all the frames to stderr
  fprintf(stderr, "\nError: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, 2);
  exit(1);
}

void installBackTracer() {
	signal(SIGSEGV, backtraceHandler);
}	

#else
#define backtraceHandler(x)	exit(x) 

void installBackTracer() {
	;
}	
#endif


void notNullFileLine(void *ptr, char *file, int line) {
	if(ptr==NULL){
		fprintf(stderr,"\nERROR null pointer. File: %s line: %d\n",file,line);
		backtraceHandler(1);
	}
}

void *allocateFileLine(size_t size, char* file, int line) {
	void *ptr = malloc(size);
	if(ptr==NULL) {
		fprintf(stderr,"\nERROR on allocate. File: %s line: %d\n",file,line);
		backtraceHandler(1);
	} else
		return ptr;
}

void deallocateFileLine(ByteP *ptr, char *file, int line) {
	if(*ptr==NULL) {
		fprintf(stderr,"\nERROR on deallocate. File: %s line: %d\n",file,line);
		backtraceHandler(1);
	} else {
		free(*ptr);
		(*ptr) = NULL;
	}
}
