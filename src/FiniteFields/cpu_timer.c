
#include "FiniteFields/cpu_timer.h"
/**************************************/

void
timer_record_start (cpu_timer *t)
{
  t->elapsed_time = 0;
  if (gettimeofday (&t->start_time, NULL) != 0)
    {
      printf ("Error in measuring time!\n");
    }
//  else
//    {
//      printf("Correct start!\n");
//   }
//  t.start_time = clock ();
} 

/**************************************/

void
timer_record_stop (cpu_timer *t)
{
  if (gettimeofday (&t->stop_time, NULL) != 0)
    {
      printf ("Error in measuring time!");
    }
//  else
//      {
//        printf("Correct stop!\n");
//      }

//  t.stop_time = clock ();
}

/**************************************/

// returns in mili-seconds
void
timer_get_elapsed_time (cpu_timer *t, const char* msg, int n_iterations)
{
  t->elapsed_time = (t->stop_time.tv_usec + t->stop_time.tv_sec * CLOCKS_PER_SEC)
      - (t->start_time.tv_usec + t->start_time.tv_sec * CLOCKS_PER_SEC);
  t->elapsed_time /= 1000;

//  n_iterations = 1;
  if (msg != NULL)
    {
      printf ("[%s]", msg);

      int length = strlen (msg);
      for (int i = length; i < MAX_TIMER_DISPLAY_WIDTH; i++)
	printf (".");
	printf ("[%.8f (ms)] \n", t->elapsed_time / (double) n_iterations);
    }

}

/**************************************/

// returns in mili-seconds
void
timer_print_time (float elapsed_time, char* msg)
{
  if (msg != NULL)
    {
      printf ("[%s]", msg);

      int length = strlen (msg);
      for (int i = length; i < MAX_TIMER_DISPLAY_WIDTH; i++)
	printf (".");
      printf ("[%.8f (ms)] \n", elapsed_time);
    }

//  if (msg != NULL)
//    {
//      printf ("[%s]:", msg);
//      printf ("[%.3f (ms)] \n", elapsed_time);
//    }

}

/**************************************/

// returns in mili-seconds
void
print_quantity(int q, char* msg)
{
  if (msg != NULL)
    {
      printf ("[%s]", msg);

      int length = strlen (msg);
      for (int i = length; i < MAX_TIMER_DISPLAY_WIDTH; i++)
	printf (".");
      printf ("[%d] \n", q);
    }

//  if (msg != NULL)
//    {
//      printf ("[%s]:", msg);
//      printf ("[%.3f (ms)] \n", elapsed_time);
//    }

}

/**************************************/

// returns in mili-seconds
void
timer_print_time_percentage (float elapsed_time, char* msg, float total_elapsed_time)
{
  if (msg != NULL)
    {
      int percentage=(1000.0*elapsed_time)/(1.0*total_elapsed_time);
      printf ("[%s]", msg);
      int length = strlen (msg);
      for (int i = length; i < MAX_TIMER_DISPLAY_WIDTH; i++)
	printf (".");
      printf ("[%.8f (ms)] (%.1f%%) \n", elapsed_time, (float)(1.0*percentage/10.0));
    }

//  if (msg != NULL)
//    {
//      printf ("[%s]:", msg);
//      printf ("[%.3f (ms)] \n", elapsed_time);
//    }

}
