/* $Header: /cvs/homes/llcbench/llcbench/cachebench/cachebench.c,v 1.1 2000/06/09 19:15:38 mucci Exp $ */
/* $Log: cachebench.c,v $
/* Revision 1.1  2000/06/09 19:15:38  mucci
/* Full blown CVS commit
/*
 * Revision 1.4  1998/03/31 05:52:14  mucci
 * Added RCS Keywords
 * */
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <sys/times.h>
#include <sys/types.h>
#include <fcntl.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <memory.h>
#include <string.h>
#include <malloc.h>
#include <sys/time.h>

static struct timeval tv1, tv2;

#define TIMER_START     go(1);gettimeofday(&tv1, (struct timezone*)0)
#define TIMER_STOP      gettimeofday(&tv2, (struct timezone*)0)
#define TIMER_ELAPSED_US	(double)(((tv2.tv_usec-tv1.tv_usec)+((tv2.tv_sec-tv1.tv_sec)*1000000)))

#define ARRAY(x,i,j,m) *(x+((j)*m)+(i))

#ifdef REGISTER
#undef REGISTER
#define REGISTER register
#else
#define REGISTER
#endif

#if defined(USE_DOUBLE)
#ifndef DATATYPE
#define DATATYPE double
#define DATASTRING "double"
#endif
#elif defined(USE_INT)
#ifndef DATATYPE
#define DATATYPE int
#define DATASTRING "int"
#endif
#elif defined(USE_CHAR)
#ifndef DATATYPE
#define DATATYPE char
#define DATASTRING "char"
#endif
#elif defined(USE_FLOAT)
#ifndef DATATYPE
#define DATATYPE float
#define DATASTRING "float"
#endif
#endif

#ifndef DATATYPE
#define DATATYPE double
#define DATASTRING "double"
#endif

#ifdef DEBUG
#define DBG(a) a;
#else
#define DBG(a)
#endif

#define READONLY 1<<0
#define WRITEONLY 1<<1
#define READWRITE 1<<2
#define NOCALIBRATE 1<<3
#define GUESSCACHESIZE 1<<4
#define MEMORYSET 1<<5
#define MEMORYCOPY 1<<6
#define HANDTUNED 1<<7

/* DEFAULT SETTINGS */

#define CACHE_MIN_BITS (8) /* 1<<bits is smallest cache in bytes */
#define CACHE_MAX_BITS (24) /* 1<<bits is largest cache in bytes */
#define CACHE_MIN (1<<CACHE_MIN_BITS)
#define CACHE_MAX (1<<CACHE_MAX_BITS)
#define RESOLUTION (2)
#define TIME_SLOTS (2*(CACHE_MAX_BITS-CACHE_MIN_BITS+1))
#define MEASURE_TIME 5
#define DOUBLE_NS(x) ((double)x)*1.0e9
#define THRESHOLD (0.9) /* 1st level must be faster than 2nd level */
#define REPEAT_COUNT 2

extern char *optarg;

int logmemsize = CACHE_MAX_BITS, memsize = 1<<CACHE_MAX_BITS, repeat_count = REPEAT_COUNT;
int resolution = RESOLUTION, timeslots, duration = MEASURE_TIME, type = NOCALIBRATE;

volatile int keepgoing;

void flushall(int yes)
{
  static char *buffer = NULL;
  static char val = 1;

  if (yes)
    {
      if (buffer == NULL)
	assert(buffer = (char *)malloc(memsize*sizeof(char)));
      memset(buffer,val++,memsize);
    }
  else
    free(buffer);
}

int go(int a)
{
  if (a==-1)
    return(keepgoing);
  else
    return(keepgoing = a);
}

void handler(int a)
{
  go(0);
}

void compute_stats(int i, int j, double refcnt, double overhead_per_ref,
		   double *times, double *bws, double *percents, int *sizes,
		   double tmicrosec, int datasize)
{
  ARRAY(times,i,j,repeat_count) = 
    ((tmicrosec*1000.0) - (refcnt*overhead_per_ref))/refcnt;

  ARRAY(bws,i,j,repeat_count) = 
    ((refcnt*datasize)/(1024.0*1024.0)) / 
    ((ARRAY(times,i,j,repeat_count)*refcnt)/1.0E9);

  if (i == 0)
    {
      if (j == 0)
	ARRAY(percents,i,j,repeat_count) = 1.0;
      else
	{
	  ARRAY(percents,i,j,repeat_count) = 
	    ARRAY(times,i,j,repeat_count) / 
	    ARRAY(times,i,j-1,repeat_count);
	}
    }
  else
    {
      if (j == 0)
	{
	  ARRAY(percents,i,j,repeat_count) = 
	    ARRAY(times,i,j,repeat_count) / 
	    ARRAY(times,i-1,repeat_count-1,repeat_count);
	}
      else
	{
	  ARRAY(percents,i,j,repeat_count) = 
	    ARRAY(times,i,j,repeat_count) / 
	    ARRAY(times,i,j-1,repeat_count);
	}
    }

  if (isatty(1))
    {
      printf("%-15d %-15.2f %-15.2f %-15.2f\n", 
	       sizes[i],
	       ARRAY(times,i,j,repeat_count),
	       ARRAY(bws,i,j,repeat_count),
	       ARRAY(percents,i,j,repeat_count));
    }
  else
    printf("%d %f\n", 
	   sizes[i],
	   ARRAY(bws,i,j,repeat_count));
  sleep(1);
}

void compute_cache_sizes(double *times, double *bws, double *percents, int *sizes)
{
  int i, maxb = 0, maxa = 0, cachea, cacheb;

  /* Look for 2 highest percentages */

  for (i=0; i < timeslots; i++)
    if (percents[i] > percents[maxa])
      maxa = i;

  for (i=0; i < timeslots; i++)
    {
      if ((i != maxa) && (percents[i] > percents[maxb]))
	{
	  maxb = i;
	}
    }

  printf("\n\t\t\tCache Capacity Analysis\n\n");
  
  if (1.0/percents[maxa] >= THRESHOLD)
    {
      printf("No L1 limit found, must be larger than %d bytes.\n",memsize);
      return;
    }
    
  /* Set them to the index if the entry for that cache size */
  /* Remember our percents are relative to the previous entry */

  cachea = sizes[maxa-1];
  cacheb = sizes[maxb-1];

  if (cachea > cacheb)
    printf("Level 1 Cache: %d bytes\n",cachea);
  else
    {
      printf("Level 1 Cache: %d bytes\n",cachea);
      if ((times[maxa]/times[maxb] < THRESHOLD) && (1.0/percents[maxb] < THRESHOLD))
	printf("Level 2 Cache: %d bytes\n",cacheb);
    }
}

int usage(int argc, char **argv)
{
  int c;
  int errflg = 0;

  while ((c = getopt(argc, argv, "m:d:hrwbcsptx:e:")) != -1)
    switch (c) {
    case 'm':
      if (((logmemsize = atoi(optarg)) < 0) || (logmemsize <= CACHE_MIN_BITS))
	errflg++;
      memsize = 1<<logmemsize;
      break;
    case 'd':
      if ((duration = atoi(optarg)) < 0)
	errflg++;
      break;
    case 's':
      type |= MEMORYSET;
      break; 
    case 'p':
      type |= MEMORYCOPY;
      break; 
    case 'r':
      type |= READONLY;
      break; 
    case 'h':
      errflg++;
      break; 
    case 't':
      type |= HANDTUNED;
      break;
    case 'w':
      type |= WRITEONLY;
      break; 
    case 'b':
      type |= READWRITE;
      break; 
    case 'c':
      type ^= NOCALIBRATE;
      break; 
    case 'g':
      type |= GUESSCACHESIZE;
      break; 
    case 'x':
      if ((resolution = atoi(optarg)) < 0)
	errflg++;
      resolution++; /* Correct for my usage */
      break; 
    case 'e':
      if ((repeat_count = atoi(optarg)) <= 0)
	errflg++;
      break; 
    case '?':
      errflg++;
      break; }

  if ((type & (READONLY|WRITEONLY|READWRITE|MEMORYCOPY|MEMORYSET)) == 0)
    {
      if (isatty(1))
	type |= READONLY|WRITEONLY|READWRITE|MEMORYCOPY|MEMORYSET;
      else
	type |= READWRITE;
    }
  
  if (errflg) 
    {
      fprintf(stderr, "Usage: %s -rwbtsp [-x #] [-m #] [-d #] [-e #]\n",argv[0]);
      fprintf(stderr, "\t -r Read benchmark\n");
      fprintf(stderr, "\t -w Write benchmark\n");
      fprintf(stderr, "\t -b Read/Modify/Write benchmark\n");
      fprintf(stderr, "\t -t Use hand tuned versions of the above\n");
      fprintf(stderr, "\t -s memset() benchmark\n");
      fprintf(stderr, "\t -p memcpy() benchmark\n");
      /* fprintf(stderr, "\t -c Enable calibration code\n"); 
      fprintf(stderr, "\t -g Enable cache size guessing code\n"); */
      fprintf(stderr, "\t -x Number of measurements between powers of 2.\n");
      fprintf(stderr, "\t -m Specify the log2(available physical memory)\n");
      fprintf(stderr, "\t -d Number of seconds per iteration\n");
      fprintf(stderr, "\t -e Repeat count per cache size\n\n");
      fprintf(stderr, "Datatype used is %s, %d bytes.\n",DATASTRING,(int)sizeof(DATATYPE));
      fprintf(stderr, "Defaults if  tty: -rwbsp -x%d -m%d -d%d -e%d\n",
	      RESOLUTION-1,CACHE_MAX_BITS,MEASURE_TIME,REPEAT_COUNT);
      fprintf(stderr, "Defaults if file: -b   -x%d -m%d -d%d -e1\n",
	      RESOLUTION-1,CACHE_MAX_BITS,MEASURE_TIME);
      exit(1);
    }

  timeslots = resolution*(logmemsize - CACHE_MIN_BITS) + 1;

  DBG(printf("%d %d %d %d\n",logmemsize,memsize,duration,timeslots)); 

  return(type);
}

void fake_out_optimizations(DATATYPE *x, int bytes)
{
  static int fd = -1;

  if (fd == -1)
    assert(fd=open("/dev/null",O_WRONLY));
  assert(write(fd,(void *)x,bytes)!=-1);
}

/* double calibrate_benchmark_ronly(REGISTER DATATYPE *x, REGISTER int limit, int *oloops, double *ous)
{
  REGISTER int index = 0, loops = 0;
  REGISTER DATATYPE sum = 0, foolem = 0;
  
  flushall(x);
  keepgoing = 1;
  assert(signal(SIGALRM,handler) != SIG_ERR);
  limit -= 4; foolem = (DATATYPE)limit;

  alarm(duration);
  TIMER_START;
again:
  sum += foolem + foolem+(DATATYPE)1 + foolem+(DATATYPE)2 + foolem+(DATATYPE)3;
  if (((index+=4) < limit) && (keepgoing))
    goto again;
  else if (keepgoing)
    {
      index = 0;
      loops++;
      goto again;
    }
  TIMER_STOP;
  index += 4;

  x[0] = (DATATYPE)sum;
  x[1] = (DATATYPE)index;
  fake_out_optimizations(x,2*sizeof(DATATYPE));

  *oloops = loops;
  *ous = TIMER_ELAPSED_US;
  return(((double)loops*(double)limit)+(double)index);  
}

double calibrate_benchmark_wonly(REGISTER DATATYPE *x, REGISTER int limit, int *oloops, double *ous)
{
  REGISTER int index = 0, loops = 0;
  REGISTER DATATYPE sum1 = (DATATYPE)0, sum2 = (DATATYPE)0, sum3 = (DATATYPE)0, sum4 = (DATATYPE)0;

  flushall(x);
  keepgoing = 1;
  assert(signal(SIGALRM,handler) != SIG_ERR);
  limit -= 4; 

  alarm(duration);
  TIMER_START;
again:
  sum1++;
  sum2++;
  sum3++;
  sum4++;
  if (((index+=4) < limit) && (keepgoing))
    goto again;
  else if (keepgoing)
    {
      index = 0;
      loops++;
      goto again;
    }
  TIMER_STOP;
  index += 4;

  x[0] = (DATATYPE)sum1;
  x[1] = (DATATYPE)sum2;
  fake_out_optimizations(x,2*sizeof(DATATYPE));
    
  return(((double)loops*(double)limit)+(double)index);
} 

double calibrate_benchmark(REGISTER DATATYPE *x, REGISTER int to_do_loops, REGISTER int limit, double *ous)
{
  REGISTER int index = 0, loops = 0;
  REGISTER DATATYPE sum1 = (DATATYPE)0;

  TIMER_START;
  while (loops < to_do_loops)
    {
      for (index = 0; index < limit; index++)
	{
	  sum1++;
	}
      loops++;
    }
  TIMER_STOP;

  x[0] = (DATATYPE)sum1;
  x[1] = (DATATYPE)index;
  fake_out_optimizations(x,2*sizeof(DATATYPE));
    
  *ous = TIMER_ELAPSED_US;
  {
    double refcnt = ((double)loops*(double)limit)+(double)index;
    DBG(fprintf(stderr,"C: %d loops at limit %d took %f us, %f refs.\n",loops,limit,*ous,refcnt));
    return(refcnt);
  }
} */

double benchmark_cache_ronly(REGISTER DATATYPE *x, REGISTER int limit, int *oloops, double *ous)
{
  REGISTER int index = 0, loops = 0;
  REGISTER DATATYPE sum = (DATATYPE)0;

  flushall(1);
  assert(signal(SIGALRM,handler) != SIG_ERR);

  alarm(duration);
  TIMER_START;

#ifdef SOLARIS
  while (go(-1))
#else
  while (keepgoing)
#endif  
    {
      for  (index = 0; index < limit; index++)
	{
	  sum += x[index];
	}
      loops++;
    }

  TIMER_STOP;

  x[0] = (DATATYPE)sum;
  x[1] = (DATATYPE)index;
  fake_out_optimizations(x,2*sizeof(DATATYPE));

  *oloops = loops;
  *ous = TIMER_ELAPSED_US;
  {
    double refcnt = ((double)loops*(double)limit)+(double)index;
    DBG(fprintf(stderr,"T: %d loops at limit %d took %f us, %f refs.\n",loops,limit,*ous,refcnt));
    return(refcnt);
  }
}

double hand_benchmark_cache_ronly(REGISTER DATATYPE *x, REGISTER int limit, int *oloops, double *ous)
{
  REGISTER int index = 0, loops = 0;
  REGISTER DATATYPE sum = (DATATYPE)0;
  REGISTER DATATYPE sum2 = (DATATYPE)0;

  flushall(1);
  keepgoing = 1;
  assert(signal(SIGALRM,handler) != SIG_ERR);
  limit -= 8;

  alarm(duration);
  TIMER_START;

again:
  sum += x[index] + x[index+1] + x[index+2] + x[index+3];
  sum2 += x[index+4] + x[index+5] + x[index+6] + x[index+7];
  if ((index += 8) < limit)
    goto again;
  else if (keepgoing)
    {
      index = 0;
      loops++;
      goto again;
    }

  TIMER_STOP;
  index += 8;

  x[0] = (DATATYPE)sum;
  x[1] = (DATATYPE)index;
  fake_out_optimizations(x,2*sizeof(DATATYPE));

  *oloops = loops;
  *ous = TIMER_ELAPSED_US;
  {
    double refcnt = ((double)loops*(double)limit)+(double)index;
    DBG(fprintf(stderr,"T: %d loops at limit %d took %f us, %f refs.\n",loops,limit,*ous,refcnt));
    return(refcnt);
  }
}

double benchmark_cache_wonly(REGISTER DATATYPE *x, REGISTER int limit, int *oloops, double *ous)
{
  REGISTER int index = 0, loops = 0;
  REGISTER DATATYPE wval = (DATATYPE)0xf;

  flushall(1);
  assert(signal(SIGALRM,handler) != SIG_ERR);
  wval = (DATATYPE)limit;

  alarm(duration);
  TIMER_START;

  while (keepgoing)
    {
      for  (index = 0; index < limit; index++)
	x[index] = wval;
      loops++;
    }

  TIMER_STOP;

  fake_out_optimizations(x,limit*sizeof(DATATYPE));

  *oloops = loops;
  *ous = TIMER_ELAPSED_US;
  {
    double refcnt = ((double)loops*(double)limit)+(double)index;
    DBG(fprintf(stderr,"T: %d loops at limit %d took %f us, %f refs.\n",loops,limit,*ous,refcnt));
    return(refcnt);
  }
}

double hand_benchmark_cache_wonly(REGISTER DATATYPE *x, REGISTER int limit, int *oloops, double *ous)
{
  REGISTER int index = 0, loops = 0;
  REGISTER DATATYPE wval = (DATATYPE)0xf;

  flushall(1);
  keepgoing = 1;
  assert(signal(SIGALRM,handler) != SIG_ERR);
  limit -= 8; wval = (DATATYPE)limit;

  alarm(duration);
  TIMER_START;

again:
  x[index] = wval;
  x[index+1] = wval;
  x[index+2] = wval;
  x[index+3] = wval;
  x[index+4] = wval;
  x[index+5] = wval;
  x[index+6] = wval;
  x[index+7] = wval;
  if ((index+=8) < limit)
    goto again;
  else if (keepgoing)
    {
      index = 0;
      loops++;
      goto again;
    }

  TIMER_STOP;
  index += 8;

  fake_out_optimizations(x,limit*sizeof(DATATYPE));

  *oloops = loops;
  *ous = TIMER_ELAPSED_US;
  {
    double refcnt = ((double)loops*(double)limit)+(double)index;
    DBG(fprintf(stderr,"T: %d loops at limit %d took %f us, %f refs.\n",loops,limit,*ous,refcnt));
    return(refcnt);
  }
}

double benchmark_cache(REGISTER DATATYPE *x, REGISTER int limit, int *oloops, double *ous)
{
  REGISTER int index = 0, loops = 0;
  
  flushall(1);
  assert(signal(SIGALRM,handler) != SIG_ERR);

  alarm(duration);
  TIMER_START;

  while (keepgoing)
    {
      for (index = 0; index < limit; index++)
	x[index]++;
      loops++;
    }

  TIMER_STOP;

  fake_out_optimizations(x,limit*sizeof(DATATYPE));

  *oloops = loops;
  *ous = TIMER_ELAPSED_US;
  {
    double refcnt = ((double)loops*(double)limit)+(double)index;
    DBG(fprintf(stderr,"T: %d loops at limit %d took %f us, %f refs.\n",loops,limit,*ous,refcnt));
    return(refcnt);
  }
}

double hand_benchmark_cache(REGISTER DATATYPE *x, REGISTER int limit, int *oloops, double *ous)
{
  REGISTER int index = 0, loops = 0;
  
  flushall(1);
  keepgoing = 1;
  assert(signal(SIGALRM,handler) != SIG_ERR);

  alarm(duration);
  TIMER_START;

again:
  x[index]++;
  x[index+1]++;
  x[index+2]++;
  x[index+3]++;
  x[index+4]++;
  x[index+5]++;
  x[index+6]++;
  x[index+7]++;
  if ((index+=8) <= limit-8)
    goto again;
  else if (keepgoing)
    {
      index = 0;
      loops++;
      goto again;
    }

  TIMER_STOP;

  fake_out_optimizations(x,limit*sizeof(DATATYPE));

  *oloops = loops;
  *ous = TIMER_ELAPSED_US;
  {
    double refcnt = ((double)loops*(double)limit)+(double)index;
    DBG(fprintf(stderr,"T: %d loops at limit %d took %f us, %f refs.\n",loops,limit,*ous,refcnt));
    return(refcnt);
  }
}


double benchmark_cache_memory_copy(REGISTER void *x, REGISTER void *y, REGISTER int bytes, int *oloops, double *ous)
{
  REGISTER int loops = 0;
  
  flushall(1);
  assert(signal(SIGALRM,handler) != SIG_ERR);

  alarm(duration);
  TIMER_START;

  while (keepgoing)
    {
      memcpy(x,y,bytes);
      loops++;
    }

  TIMER_STOP;

  fake_out_optimizations(x,bytes);
  fake_out_optimizations(y,bytes);

  *ous = TIMER_ELAPSED_US;
  *oloops = loops;
  return((double)loops*(double)bytes);
}

double benchmark_cache_memory_set(REGISTER void *x, REGISTER int bytes, int *oloops, double *ous)
{
  REGISTER int loops = 0;
  
  flushall(1);
  assert(signal(SIGALRM,handler) != SIG_ERR);

  alarm(duration);
  TIMER_START;

  while (keepgoing)
    {
      memset(x,0xf0,bytes);
      loops++;
    }

  TIMER_STOP;

  fake_out_optimizations(x,bytes);

  *ous = TIMER_ELAPSED_US;
  *oloops = loops;
  return((double)loops*(double)bytes);
}

void initialize_sizes(int *sizes)
{
  int i,j;
  
  for (j=0; j<timeslots; j+=resolution)
    {
      sizes[j] = 1<<(CACHE_MIN_BITS+j/resolution);
      DBG(printf("POW: %d %d\n",j,sizes[j]));

      for (i=1;i<resolution;i++)
	{
	  if (j+i < timeslots)
	    {
	      sizes[j+i] = sizes[j] + i*(sizes[j]/resolution);
	      sizes[j+i] = sizes[j+i] - sizes[j+i]%(int)sizeof(DATATYPE);
	      DBG(printf("SUB: %d %d\n",j+i,sizes[j+i]));
	    }
	}
    }
}

void do_memory_copy(int *sizes, void *x, double *times, double *bws, double *percents)
{
  int limit, j, i, tloops;
  double refcnt, overhead_per_ref = 0.0, tmicrosec;
  /* double nullcnt = 0.0, cmicrosec = 0.0; */
  void *y;

  assert(y = (void *)malloc(memsize));
  memset(y,0x0f,memsize);

  if (isatty(1))
    {
      printf("\n\t\t%sMemory Copy Library Cache Test\n\n",
	     (!(type & NOCALIBRATE)) ? "Calibrated " : "");
      printf("C Size\t\tNanosec\t\tMB/sec\t\t%% Chnge\n");
      printf("-------\t\t-------\t\t-------\t\t-------\n");
    }

  for (i = 0; i < timeslots; i++) 
    {
      limit = sizes[i];
      for (j = 0; j < repeat_count; j++)
	{
	  refcnt = benchmark_cache_memory_copy(x, y, limit, &tloops, &tmicrosec) * 2.0;
	  
	  /* if (type & NOCALIBRATE)
	    {
	      nullcnt = 0.0;
	      overhead_per_ref = 0.0;
	      cmicrosec = 0.0;
	    }
	  else
	    {
	      nullcnt = calibrate_benchmark_cache_memory_copy(x, tloops, limit, &cmicrosec);
	      overhead_per_ref = (cmicrosec*1000.0) / nullcnt;
	      DBG(fprintf(stderr,"C: %f ns per ref.\n",overhead_per_ref));
	    } */

	  compute_stats(i,j,refcnt,overhead_per_ref,times,bws,percents,sizes,tmicrosec,1);
	}
    }
  
  /* if ((isatty(1))&&(type&GUESSCACHESIZE))
    compute_cache_sizes(times,bws,percents,sizes); */

  free(y);
}

void do_memory_set(int *sizes, DATATYPE *x, double *times, double *bws, double *percents)
{
  int limit, j, i, tloops;
  double refcnt, overhead_per_ref = 0.0, tmicrosec;
  /* double nullcnt = 0.0, cmicrosec = 0.0; */

  if (isatty(1))
    {
      printf("\n\t\t%sMemory Set Library Cache Test\n\n",
	     ((!(type & NOCALIBRATE)) ? "Calibrated " : ""));

      printf("C Size\t\tNanosec\t\tMB/sec\t\t%% Chnge\n");
      printf("-------\t\t-------\t\t-------\t\t-------\n");
    }

  for (i = 0; i < timeslots; i++) 
    {
      limit = sizes[i];
      for (j = 0; j < repeat_count; j++)
	{
	  refcnt = benchmark_cache_memory_set(x, limit, &tloops, &tmicrosec);
	  
	  /* if (type & NOCALIBRATE)
	    {
	      nullcnt = 0.0;
	      overhead_per_ref = 0.0;
	      cmicrosec = 0.0;
	    }
	  else
	    {
	      nullcnt = calibrate_benchmark_cache_memory_copy(x, tloops, limit, &cmicrosec);
	      overhead_per_ref = (cmicrosec*1000.0) / nullcnt;
	      DBG(fprintf(stderr,"C: %f ns per ref.\n",overhead_per_ref));
	    } */

	  compute_stats(i,j,refcnt,overhead_per_ref,times,bws,percents,sizes,tmicrosec,1);
	}
    }
  
  /* if ((isatty(1))&&(type&GUESSCACHESIZE))
    compute_cache_sizes(times,bws,percents,sizes); */
}

void do_read_only(int *sizes, DATATYPE *x, double *times, double *bws, double *percents)
{
  int limit, j, i, tloops;
  double refcnt, overhead_per_ref = 0.0, tmicrosec;
  /* double nullcnt = 0.0, cmicrosec = 0.0; */

  if (isatty(1))
    {
      printf("\n\t\t%s%s%s Read Cache Test\n\n",
	     ((type & HANDTUNED) ? "Tuned " : ""),
	     ((!(type & NOCALIBRATE)) ? "Calibrated " : ""), DATASTRING);

      printf("C Size\t\tNanosec\t\tMB/sec\t\t%% Chnge\n");
      printf("-------\t\t-------\t\t-------\t\t-------\n");
    }

  for (i = 0; i < timeslots; i++) 
    {
      limit = sizes[i] / (int)sizeof(DATATYPE);
      for (j = 0; j < repeat_count; j++)
	{
	  if (type & HANDTUNED)
	    refcnt = hand_benchmark_cache_ronly(x, limit, &tloops, &tmicrosec);
	  else
	    refcnt = benchmark_cache_ronly(x, limit, &tloops, &tmicrosec);

	  /* if (type & NOCALIBRATE)
	    {
	      nullcnt = 0.0;
	      overhead_per_ref = 0.0;
	      cmicrosec = 0.0;
	    }
	  else
	    {
	      nullcnt = calibrate_benchmark(x, tloops, limit, &cmicrosec);
	      overhead_per_ref = (cmicrosec*1000.0) / nullcnt;
	      DBG(fprintf(stderr,"C: %f ns per ref.\n",overhead_per_ref));
	    } */

	  compute_stats(i,j,refcnt,overhead_per_ref,times,bws,percents,sizes,tmicrosec,sizeof(DATATYPE));
	}
    }
      
  /* if ((isatty(1))&&(type&GUESSCACHESIZE))
    compute_cache_sizes(times,bws,percents,sizes); */
}

void do_write_only(int *sizes, DATATYPE *x, double *times, double *bws, double *percents)
{
  int limit, j, i, tloops;
  double refcnt, overhead_per_ref = 0.0, tmicrosec;
  /* double nullcnt = 0.0, cmicrosec = 0.0; */

  if (isatty(1))
    {
      printf("\n\t\t%s%s%s Write Cache Test\n\n",
	     ((type & HANDTUNED) ? "Tuned " : ""),
	     ((!(type & NOCALIBRATE)) ? "Calibrated " : ""), DATASTRING);

      printf("C Size\t\tNanosec\t\tMB/sec\t\t%% Chnge\n");
      printf("-------\t\t-------\t\t-------\t\t-------\n");
    }

  for (i = 0; i < timeslots; i++) 
    {
      limit = sizes[i] / (int)sizeof(DATATYPE);
      for (j = 0; j < repeat_count; j++)
	{
	  if (type & HANDTUNED)
	    refcnt = hand_benchmark_cache_wonly(x, limit, &tloops, &tmicrosec);
	  else
	    refcnt = benchmark_cache_wonly(x, limit, &tloops, &tmicrosec);

	  /* if (type & NOCALIBRATE)
	    {
	      nullcnt = 0.0;
	      overhead_per_ref = 0.0;
	      cmicrosec = 0.0;
	    }
	  else
	    {
	      nullcnt = calibrate_benchmark(x, tloops, limit, &cmicrosec);
	      overhead_per_ref = (cmicrosec*1000.0) / nullcnt;
	      DBG(fprintf(stderr,"C: %f ns per ref.\n",overhead_per_ref));
	    } */
	  
	  compute_stats(i,j,refcnt,overhead_per_ref,times,bws,percents,sizes,tmicrosec,sizeof(DATATYPE));
	}
    }
  
  /* if ((isatty(1))&&(type&GUESSCACHESIZE))
    compute_cache_sizes(times,bws,percents,sizes); */
}

void do_read_write(int *sizes, DATATYPE *x, double *times, double *bws, double *percents)
{
  int limit, j, i, tloops;
  double refcnt, overhead_per_ref = 0.0, tmicrosec;
  /* double nullcnt = 0.0, cmicrosec = 0.0; */

  if (isatty(1))
    {
      printf("\n\t\t%s%s%s RMW Cache Test\n\n",
	     ((type & HANDTUNED) ? "Tuned " : ""),
	     ((!(type & NOCALIBRATE)) ? "Calibrated " : ""), DATASTRING);

      printf("C Size\t\tNanosec\t\tMB/sec\t\t%% Chnge\n");
      printf("-------\t\t-------\t\t-------\t\t-------\n");
    }

  for (i = 0; i < timeslots; i++) 
    {
      limit = sizes[i] / (int)sizeof(DATATYPE);
      for (j = 0; j < repeat_count; j++)
	{
	  if (type & HANDTUNED)
	    refcnt = hand_benchmark_cache(x, limit, &tloops, &tmicrosec) * 2.0;
	  else
	    refcnt = benchmark_cache(x, limit, &tloops, &tmicrosec) * 2.0;
  
	  /* if (type & NOCALIBRATE)
	    {
	      nullcnt = 0.0;
	      overhead_per_ref = 0.0;
	      cmicrosec = 0.0;
	    }
	  else
	    {
	      nullcnt = calibrate_benchmark(x, tloops, limit, &cmicrosec);
	      nullcnt *= 2; 
	      overhead_per_ref = (cmicrosec*1000.0) / nullcnt;
	      DBG(fprintf(stderr,"C: %f ns per ref.\n",overhead_per_ref));
	    } */
	  
	  compute_stats(i,j,refcnt,overhead_per_ref,times,bws,percents,sizes,tmicrosec,sizeof(DATATYPE));
	}
    }
  
  /* if ((isatty(1))&&(type&GUESSCACHESIZE))
    compute_cache_sizes(times,bws,percents,sizes); */
}

int main(int argc, char **argv) 
{
  DATATYPE *x;
  int *sizes;
  double *times, *bws, *percents;

  type = usage(argc, argv);

  assert(sizes = (int *)malloc(timeslots*sizeof(int)));
  memset(sizes,0x00,(timeslots*sizeof(int)));
  assert(times = (double *)malloc(timeslots*repeat_count*sizeof(double)));
  memset(times,0x00,(timeslots*repeat_count*sizeof(double)));
  assert(bws = (double *)malloc(timeslots*repeat_count*sizeof(double)));
  memset(bws,0x00,(timeslots*repeat_count*sizeof(double)));
  assert(percents = (double *)malloc(timeslots*repeat_count*sizeof(double)));
  memset(percents,0x00,(timeslots*repeat_count*sizeof(double)));
  assert(x = (DATATYPE *)malloc(memsize));
  memset((void *)x,0x00,memsize);

  initialize_sizes(sizes);

  /* Measure cache */

  if (type & MEMORYSET)
    {
      do_memory_set(sizes,x,times,bws,percents);
    }

  if (type & MEMORYCOPY)
    {
      do_memory_copy(sizes,x,times,bws,percents);
    }

  if (type & READONLY)
    {
      do_read_only(sizes,x,times,bws,percents);
    }

  if (type & WRITEONLY)
    {
      do_write_only(sizes,x,times,bws,percents);
    }

  if (type & READWRITE)
    {
      do_read_write(sizes,x,times,bws,percents);
    }

  flushall(0);
  exit(0);
}

