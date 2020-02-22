
/**
 * A simple timer for calculating user time.
 */

#ifndef _UNIX_TIMER_H_
#define _UNIX_TIMER_H_

#ifdef __cplusplus
extern "C" { 
#endif

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

// For execution timing
#if !defined(SERIAL) && defined(CILKVIEW_TIMING)
	#include <cilktools/cilkview.h>
#endif

typedef struct timeval timer_time;
typedef struct {
	timer_time start_time;
} timer_id;

static timer_time user_time() {
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);

	return usage.ru_utime;
}

static timer_id start_timer() {
	timer_time start_time = user_time();
	timer_id id;
	id.start_time = start_time;
	return id;
}

static timer_time elapsed_time(timer_id* timerID) {
	if (timerID == NULL) {
		return user_time();
	}

	timer_time t = user_time();
	t.tv_sec -= timerID->start_time.tv_sec;
	t.tv_usec -= timerID->start_time.tv_usec; 
	return t;
}





static inline void _startTimer(unsigned long long *start){
#if defined(SERIAL) && SERIAL
        timer_id id = start_timer();
        *start = (unsigned long long) ((id.start_time.tv_sec)*1000000ll) + id.start_time.tv_usec;
        // *start = (unsigned long long)clock();
#else
	    struct timeval t;
	    gettimeofday(&t, 0);
        *start = (unsigned long long) ((t.tv_sec)*1000000ll) + t.tv_usec;
#endif
}

static inline void _stopTimer(unsigned long long *start, float *elapsed){
#if defined(SERIAL) && SERIAL
	timer_id now = start_timer();
	unsigned long long nowll = (unsigned long long) ((now.start_time.tv_sec)*1000000ll) + now.start_time.tv_usec;
	nowll -= *start;
	*elapsed = (float) (nowll / 1000000.0);
#else
    struct timeval now;
    gettimeofday(&now, 0);
	unsigned long long nowll = (unsigned long long) ((now.tv_sec)*1000000ll) + now.tv_usec;
	// nowll -= *start;
	*elapsed = (float) ((nowll - *start) / 1000000.0);
	*start = nowll;

#endif
}

static inline void _stopTimerAddElapsed(unsigned long long *start, float *elapsed){
#if defined(SERIAL) && SERIAL
	timer_id now = start_timer();
	unsigned long long nowll = (unsigned long long) ((now.start_time.tv_sec)*1000000ll) + now.start_time.tv_usec;
	nowll -= *start;
	*elapsed += (float) (nowll / 1000000.0);
#else
    struct timeval now;
    gettimeofday(&now, 0);
	unsigned long long nowll = (unsigned long long) ((now.tv_sec)*1000000ll) + now.tv_usec;
	nowll -= *start;
	*elapsed += (float) (nowll / 1000000.0);
#endif
}








#ifdef __cplusplus
} 
#endif

#endif

