/*
 * timer.h
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.  
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <sys/time.h>
#include <unistd.h>

typedef struct timeval Clock;

typedef enum {
	STOPPED = 0,
	TIMER_RUNNING = 1,
	PAUSED  = 2
} TimerStatus;

typedef struct {
	Clock beginTime;	 	/* time at which the timer was started */
	double elapsedTime; 	/* time elapsed since last start */
	TimerStatus status;	 	/* status of the timer */
} Timer;

typedef Timer *TimerP;


TimerP newTimer();

void freeTimer(TimerP *t);

void ResetTimer (TimerP t);

void StartTimer (TimerP t);

void PauseTimer (TimerP t);

void StopTimer (TimerP t);

double GetElapsedTime (Timer *t);

#endif /* TIMER_H_ */
