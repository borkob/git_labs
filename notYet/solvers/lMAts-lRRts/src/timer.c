/*
 * timer.c
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/

#include "timer.h"
#include "dynamicMem.h"

#define GetTimeFunction(t)   {  gettimeofday(&t, NULL); }

TimerP newTimer() {
	TimerP t = (TimerP) allocate(sizeof(Timer));
	ResetTimer(t);
	return t;
}


void freeTimer(TimerP *t) {
	deallocate2(t);
}
 

double TimeDifference(Clock end, Clock start)
{ long seconds  = end.tv_sec  - start.tv_sec;
  long useconds = end.tv_usec - start.tv_usec;
  return  ((seconds) + (double)useconds/1E6) ;
}

/*
 *  Resets a timer, i.e., stops it if it was running and
 * sets the elapsed time to 0.
 *
 */
void ResetTimer (Timer* t)
{
	t->status = STOPPED;
	t->elapsedTime = 0.0;
}

/*
 *  Starts the timer. If it was stopped, the elapsed time
 * is reset to 0. Otherwise (it was paused) the elapsed time
 * is retained.
 *
 */
void StartTimer (Timer* t)
{
	if (t->status == STOPPED)
		t->elapsedTime = 0.0;
	t->status = TIMER_RUNNING;
	GetTimeFunction(t->beginTime);
}




/*
 *  Pauses a timer. Any further call to GetElapsedTime will
 * return the total time measured since it was started. The
 * timer can be further resumed by calling StartTimer.
 *
 */
void PauseTimer (Timer* t)
{
	Clock timeNow;

	GetTimeFunction(timeNow);
	t->elapsedTime += TimeDifference(timeNow, t->beginTime);
	t->status = PAUSED;
}




/*
 *  Stops a timer. This function is similar to PauseTimer,
 * with the diference that the timer cannot be resumed.
 *
 */
void StopTimer (Timer* t)
{
	Clock timeNow;

	GetTimeFunction(timeNow);
	t->elapsedTime += TimeDifference(timeNow, t->beginTime);
	t->status = STOPPED;
}



/*
 *  Gets the elapsed time. This function works in any status
 * of the timer but notice that this function requires some
 * time and, if the timer is running, the final measure will
 * include it.
 *
 */
double GetElapsedTime(Timer *t)
{
	if (t->status == TIMER_RUNNING) {
		Clock timeNow;

		GetTimeFunction(timeNow);
		return (t->elapsedTime + TimeDifference(timeNow, t->beginTime));
	}

	return t->elapsedTime;
}
