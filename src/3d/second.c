#include <stdio.h>
#include <sys/time.h>

double second();   /* Use for IBM processors. */
double second_();  /* Use for Intel processors. */

double second(void)
{
  struct timeval timestr;
  void *Tzp=0;
  gettimeofday(&timestr, Tzp);
  return (double)timestr.tv_sec+1.0E-06*(double)timestr.tv_usec;
}

double second_(void)
{
  struct timeval timestr;
  void *Tzp=0;
  gettimeofday(&timestr, Tzp);
  return (double)timestr.tv_sec+1.0E-06*(double)timestr.tv_usec;
}
