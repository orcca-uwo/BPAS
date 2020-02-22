#ifndef CPU_TIMER_H
#define CPU_TIMER_H

#define CPU_TIMER_TYPE 1

#define MAX_TIMER_DISPLAY_WIDTH 40
/**************************************/

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
//#include <>

typedef struct cpu_timer
{
  struct timeval start_time, stop_time;
  double elapsed_time;
//	int is_enabled = 0;
} cpu_timer;

/**************************************/

void
timer_record_start (cpu_timer *t);

/**************************************/

void
timer_record_stop (cpu_timer *t);
/**************************************/

// returns in mili-seconds
void
timer_get_elapsed_time (cpu_timer *t, const char* msg, int n_iterations);
/**************************************/

// returns in mili-seconds
void
timer_print_time (float elapsed_time, char* msg);

/**************************************/

// returns in mili-seconds
void
print_quantity(int q, char* msg);

/**************************************/

// returns in mili-seconds
void
timer_print_time_percentage (float elapsed_time, char* msg,
			     float total_elapsed_time);

/**************************************/
#endif
