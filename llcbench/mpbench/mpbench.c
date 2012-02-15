#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <sys/time.h>
#include <sys/signal.h>

#include "mpi.h"

/* Macros */

#ifdef INLINE
#undef INLINE
#define INLINE inline
#else
#define INLINE
#endif

#ifdef DEBUG
#define DBG(a) a;
#else
#define DBG(a)
#endif

/* DEFAULT SETTINGS */

#define LOG_MAX_MESSAGE_SIZE 16 
#define LOG_MIN_MESSAGE_SIZE 2 
#define MAX_MESSAGE_SIZE (1 << LOG_MAX_MESSAGE_SIZE)
#define RESOLUTION 1      
#define TIMESLOTS (RESOLUTION * (LOG_MAX_MESSAGE_SIZE - LOG_MIN_MESSAGE_SIZE) +1)
#define REPEAT_COUNT 1
#define ITERATIONS 1

/* Definitions */ 

#define TEST_EMPTY      0
#define TEST_LATENCY	1<<0
#define TEST_ROUNDTRIP	1<<1
#define TEST_BANDWIDTH	1<<2
#define TEST_ALLTOALL	1<<3
#define TEST_BROADCAST	1<<4
#define TEST_REDUCE	1<<5
#define TEST_ALLREDUCE  1<<6
#define TEST_COLLECTIVE 1<<7
#define TEST_BIBANDWIDTH 1<<8

#define TIMER_CLEAR     (tv1.tv_sec = tv1.tv_usec = tv2.tv_sec = tv2.tv_usec = 0)
#define TIMER_START     gettimeofday(&tv1, (struct timezone*)0)
#define TIMER_ELAPSED	((tv2.tv_usec-tv1.tv_usec)+((tv2.tv_sec-tv1.tv_sec)*1000000))
#define TIMER_STOP      gettimeofday(&tv2, (struct timezone*)0)

#define BROADCASTRESPONSE_TAG    3

/* Globals */

extern char *optarg;

int logmemsize = LOG_MAX_MESSAGE_SIZE, memsize = MAX_MESSAGE_SIZE, repeat_count = REPEAT_COUNT;
int resolution = RESOLUTION, timeslots = TIMESLOTS;
int *sizes;

char *sendbuf = NULL;
char *destbuf = NULL;
struct timeval tv1, tv2;
int type = TEST_EMPTY;
int max_rank = -1;
int my_rank = -1;
int source_rank = -1;
int dest_rank = -1;

/*************************/
/*  Prototypes           */
/*************************/

void flushall(int);

/*************************/
/*  Macros               */
/*************************/

INLINE int mp_broadcast(void *data,int bytes)
{
  return(MPI_Bcast(data,bytes,MPI_BYTE,0,MPI_COMM_WORLD));
}
INLINE int mp_reduce(void *data,int bytes)
{
  return(MPI_Reduce(data, destbuf, (bytes/(int)sizeof(int)), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD));
}
INLINE int mp_allreduce(void *data,int bytes)
{
  return(MPI_Allreduce(data, destbuf, (bytes/(int)sizeof(int)), MPI_INT, MPI_SUM, MPI_COMM_WORLD));
}
INLINE int mp_alltoall(void *sendbuffer, void *receivebuffer, int bytes)
{
  return(MPI_Alltoall(sendbuffer, (bytes/max_rank), MPI_BYTE, receivebuffer, (bytes/max_rank), MPI_BYTE, MPI_COMM_WORLD));
}
INLINE int mp_send(int dest, int tag, void *data, int bytes)
{
  return(MPI_Send(data, bytes, MPI_BYTE, dest, tag, MPI_COMM_WORLD));
}
INLINE int mp_recv(int dest,int tag, void *data,int bytes)
{
  MPI_Status stat; 
  return(MPI_Recv(data, bytes, MPI_BYTE, dest, tag, MPI_COMM_WORLD, &stat));
}
INLINE int mp_any_recv(int tag,void *data, int bytes)
{
  return(mp_recv(MPI_ANY_SOURCE, tag, data, bytes));
}
INLINE int mp_isend(int dest, int tag, void *data, int bytes, MPI_Request *request)
{
  return(MPI_Isend(data, bytes, MPI_BYTE, dest, tag, MPI_COMM_WORLD, request));
}
INLINE int mp_irecv(int dest, int tag, void *data, int bytes, MPI_Request *request)
{
  return(MPI_Irecv(data, bytes, MPI_BYTE, dest, tag, MPI_COMM_WORLD, request));
}
INLINE int mp_exit(void)
{
  return(MPI_Finalize());
}
INLINE int mp_barrier(void)
{
  return(MPI_Barrier(MPI_COMM_WORLD));
}
#if 0
INLINE int mp_wait(void)
{
}
#endif

/*************************/
/*  Utility functions    */
/*************************/

int am_i_the_master(void)
{
  return(my_rank == 0);
}

int am_i_the_slave(void)
{
  return(my_rank == max_rank-1);
}

void my_two_printf(char *str,unsigned int a, float b)
{
  if (am_i_the_master())
    printf(str,a,b);
}


void initialize_sizes(void)
{
  int i,j;

  assert(sizes = (int *)malloc(timeslots*sizeof(int)));
  memset(sizes,0x00,(timeslots*sizeof(int)));
 
  
  for (j=0; j < timeslots; j+=resolution)
    {
      sizes[j] = 1 << (LOG_MIN_MESSAGE_SIZE + j/resolution);
      DBG(printf("POW: %d %d\n",j,sizes[j]));
      for (i=1; i<resolution; i++)
	{
	  if (j+i < timeslots)
	    {
	      sizes[j+i] = sizes[j] + i*(sizes[j]/resolution);
	      sizes[j+i] = sizes[j+i] - sizes[j+i]%(int)sizeof(int);
	      DBG(printf("SUB: %d %d\n",j+i,sizes[j+i]));
	    }
	}
    }
}

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
  else{
    free(buffer);
    buffer = NULL;
  }
}

/*************************/
/*  Benchmarks           */
/*************************/

/* ALso known as the gap time by Berkeley, 
   the time to launch a message in the network's buffers. */

float latency(int cnt, int bytes)
{
  int i;
  float total = 0.0;

  if (am_i_the_master())
    {
      TIMER_START;
      for (i=0; i<cnt; i++){
	mp_send(dest_rank, 1, sendbuf, bytes);
      }
      TIMER_STOP;
      mp_recv(dest_rank, 2, destbuf, 4);
      total += (float)TIMER_ELAPSED;
      return(total/(float)cnt);   
    }
  else if (am_i_the_slave())
    {
      for (i=0; i<cnt; i++){
	mp_recv(source_rank, 1, destbuf, bytes);
      }
      mp_send(source_rank, 2, sendbuf, 4); 
      return(0.0);
    }
  else
    return(0.0);
}

/* This might be your more familiar definition of latency...*/

float roundtrip(int cnt, int bytes)
{
  int i;
  float total = 0.0;

  if (am_i_the_master())
    {
      TIMER_START;
      for (i=0; i<cnt; i++)
	{
	  mp_send(dest_rank, 1, sendbuf, bytes);
	  mp_recv(dest_rank, 2, destbuf, bytes);
	}
      TIMER_STOP;
      total += (float)TIMER_ELAPSED;
      return((float)cnt / (total*1.0E-6)); /* Transactions/sec */
    }
  else if (am_i_the_slave())
    {
      for (i=0; i<cnt; i++)
	{
	  mp_recv(source_rank, 1, destbuf, bytes);
	  mp_send(source_rank, 2, sendbuf, bytes);
	}
      return(0.0);
    }
  else
    return(0.0);
}

float bandwidth(int cnt, int bytes)
{
  int i;
  float total = 0.0;

  if (am_i_the_master())
    {
      TIMER_START;
      for (i=0; i<cnt; i++)
	mp_send(dest_rank, 1, sendbuf, bytes);
      mp_recv(dest_rank, 2, destbuf, 4);
      TIMER_STOP;
      total += (float)TIMER_ELAPSED;
      return(((float)(4+(bytes*cnt))/1024.0) / (total*1.0E-6)); /* KB/sec */
    }
  else if (am_i_the_slave())
    {
      for (i=0; i<cnt; i++)
	mp_recv(source_rank, 1, destbuf, bytes);
      mp_send(source_rank, 2, sendbuf, 4);
      return(0.0);
    }
  else
    return(0.0);
}

/*  New bidirectional bandwidth test */

float bibandwidth(int cnt, int bytes)
{
  int i;
  float total = 0.0;
  MPI_Request requestarray[2];
  MPI_Status  statusarray[2];

  if (am_i_the_master()) 
    {
      TIMER_START;
      for (i=0; i<cnt; i++)
	{
	  mp_irecv(dest_rank, 2, destbuf, bytes, &requestarray[1]);
	  mp_isend(dest_rank, 1, sendbuf, bytes, &requestarray[0]);
	  MPI_Waitall(2, requestarray, statusarray);
	} 
      TIMER_STOP;
      total += (float)TIMER_ELAPSED;
      return(((float)((2.0*(float)bytes*(float)cnt))/1024.0) / (total*1.0E-6)); /* KB/sec */
    }
  else if (am_i_the_slave())
    {
      /* This is the last process */
      for (i=0; i<cnt; i++)
	{
	  mp_irecv(source_rank, 1, destbuf, bytes, &requestarray[0]);
	  mp_isend(source_rank, 2, sendbuf, bytes, &requestarray[1]);
	  MPI_Waitall(2, requestarray, statusarray);
	}
      return(0.0);
    }
  else
    return(0.0);
}

float alltoall(int cnt, int bytes)
{
  int i;
  float total = 0.0;
  
  TIMER_START;
  for (i=0;i < cnt;i++)
    {
      mp_alltoall(sendbuf, destbuf, bytes);
    }
  TIMER_STOP;
  total += (float)TIMER_ELAPSED;
  return(((float)(cnt*bytes)/1024.0) / (total*1.0E-6)); /* KB/sec */
}

float reduce(int cnt, int bytes)
{
  int i;
  float total = 0.0;

  TIMER_START;
  for (i=0; i<cnt; i++)
    {
      mp_reduce(sendbuf, bytes);
    }
  TIMER_STOP;
  total += (float)TIMER_ELAPSED;
  return(((float)(cnt*bytes)/1024.0) / (total*1.0E-6)); /* KB/sec */
}

float allreduce(int cnt, int bytes)
{
  int i;
  float total = 0.0;

  TIMER_START;
  for (i=0; i<cnt; i++)
    {
      mp_allreduce(sendbuf, bytes);
    }
  TIMER_STOP;
  total += (float)TIMER_ELAPSED;
  return(((float)(cnt*bytes)/1024.0) / (total*1.0E-6)); /* KB/sec */
}

float broadcast(int cnt, int bytes)
{
  int i;
  float total = 0.0;

  if (am_i_the_master())
    {
      TIMER_START;
      for (i=0; i<cnt; i++)
	mp_broadcast(sendbuf, bytes);
      for (i=1; i<max_rank; i++)
	mp_any_recv(BROADCASTRESPONSE_TAG, destbuf, 4);
      TIMER_STOP;
      total += (float)TIMER_ELAPSED;
      return((((float)(cnt*bytes)+(float)(max_rank*4))/1024.0) / (total*1.0E-6)); /* KB/sec */
    }
  else
    {  
      for (i=0; i<cnt; i++)
	mp_broadcast(destbuf, bytes);
      mp_send(source_rank, BROADCASTRESPONSE_TAG, destbuf, 4);
      return(0.0);
    }
}

void loop(int iterations, float (*test)(int, int))
{
  int i, j;

  /* Allow routing/cache setup ahead of time */
  test(1,memsize);
  
  for( i=0; i < timeslots; i++){
    for( j=1; j <= repeat_count; j++){
      flushall(1);
      mp_barrier();
      my_two_printf("%u %f\n",sizes[i],test(iterations, sizes[i]));
      mp_barrier();
    }
  }
  flushall(0);
}

/*******************************/
/* new --- added by JT 7/98 */
/*******************************/

void usage(int argc, char **argv, int *iterations)
{
  int c;
  int errflg = 0;
  int i;  /* counter */

  /* I put earlier to get rid of extraneous command line stuff from MPI */
   assert(MPI_Init(&argc, &argv) == 0);
   MPI_Comm_size(MPI_COMM_WORLD, &max_rank);
   if (max_rank < 2){
      fprintf(stderr,"\nMust have at least 2 tasks.\n");
      MPI_Finalize();
      exit(1);
   }

   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   source_rank = 0;
   dest_rank = max_rank-1;

   DBG(
     for (i=0; i < argc; i++)
       printf("argc: %d   rank: %d    argv[%d]: %s\n", argc, rank, i, argv[i]);
   )

 
  while ((c = getopt(argc, argv, "i:m:x:e:bdlracyzh")) != -1)
    switch (c) {   
    case 'i':
      if ((*iterations = atoi(optarg)) < 0)
        errflg++;
      break;
    case 'm':
      if ((logmemsize = atoi(optarg)) < 0)
	errflg++;
      memsize = 1<<logmemsize;
      break;
    case 'x':
      if ((resolution = atoi(optarg)) < 0)
	errflg++;
      resolution++; /* Correct for my usage, since zero measurements between powers of 2 is */
      break;        /* a resolution of 1 */
    case 'e':
      if ((repeat_count = atoi(optarg)) <= 0)
	errflg++;
      break;
    case 'b':
      type = TEST_BANDWIDTH;
      break;
    case 'd':
      type = TEST_BIBANDWIDTH;
      break;
    case 'l':
      type = TEST_LATENCY;
      break;
    case 'r':
      type = TEST_ROUNDTRIP;
      break;
    case 'a':
      type = TEST_ALLTOALL;
      break;
      /* there aren't easy to remember */
    case 'c':
      type = TEST_BROADCAST;
      break;
    case 'y':
      type = TEST_REDUCE;
      break;
    case 'z':
      type = TEST_ALLREDUCE;
      break;
    case 'h': /* prints usage */
      errflg++;
      break;
    case '?': /* error condition */
      errflg++;
      break; }
  
  if (type == TEST_EMPTY)
    errflg++;

  if (*iterations < 1)
    errflg++;

  if (errflg) {
    if(am_i_the_master()){
      fprintf(stderr, "Usage: %s -blracyz [-i #] [-x #] [-m #] [-d #] [-e #]\n",argv[0]);
      fprintf(stderr, "\t -b Do bandwidth benchmark\n");
      fprintf(stderr, "\t -d Do bidirectional bandwidth benchmark\n");
      fprintf(stderr, "\t -l Do latency benchmark\n");
      fprintf(stderr, "\t -r Do roundtrip benchmark\n"); 
      fprintf(stderr, "\t -a Do all-to-all benchmark\n"); 
      fprintf(stderr, "\t -c Do broadcast benchmark\n"); 
      fprintf(stderr, "\t -y Do reduce benchmark\n"); 
      fprintf(stderr, "\t -z Do allreduce benchmark\n"); 
      fprintf(stderr, "\t -i Specify the iterations over which to average. \n");
      fprintf(stderr, "\t -x Specify the number of measurements between powers of 2.  \n");
      fprintf(stderr, "\t -m Specify the log2(available physical memory) to be used \n\t    as the maximum message size.\n");
      fprintf(stderr, "\t -e Specify the repeat count per message size. \n");
    }  
    mp_exit();
    exit(1);
  }

  timeslots = resolution * (logmemsize - LOG_MIN_MESSAGE_SIZE) + 1;

  /* Allocate space for send and receive buffers. */
  assert((sendbuf = (char *)malloc(memsize)));
  assert((destbuf = (char *)malloc(memsize)));
  for (i=0; i < memsize; i++){
    sendbuf[i] = destbuf[i] = (char)(i&0xff);
  }
 
  DBG(printf("logmemsize = %d, memsize =  %d, timeslots =  %d\n",logmemsize,memsize,timeslots)); 
}

/********************/
/*  main           **/
/********************/

int main(int argc, char **argv)
{
  int iterations = ITERATIONS;

  /* set up the run */
  setbuf(stdout,NULL);
  usage(argc, argv, &iterations);
  initialize_sizes();

  /* run appropriate test */
  if (iterations) {
      switch (type)
	{
	case TEST_LATENCY:
	  loop(iterations,latency);
	  break;
	case TEST_ROUNDTRIP:
	  loop(iterations,roundtrip);
	  break;
	case TEST_BANDWIDTH:
	  loop(iterations,bandwidth);
	  break;
	case TEST_BROADCAST:
	  loop(iterations, broadcast);
	  break;
	case TEST_REDUCE:
	  loop(iterations, reduce);
	  break;
	case TEST_ALLREDUCE:
	  loop(iterations, allreduce);
	  break;
	case TEST_ALLTOALL:
	  loop(iterations, alltoall);
	  break;
	case TEST_BIBANDWIDTH:
	  loop(iterations, bibandwidth);
	  break;
	default:
	  abort();
	}
  }
  free(sendbuf);
  free(destbuf);
  mp_exit();
  exit(0);
}
