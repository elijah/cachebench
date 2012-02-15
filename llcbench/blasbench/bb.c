#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <memory.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#ifdef T3E
#include <fortran.h>
#define sgemv_ HGEMV
#define dgemv_ SGEMV
#define saxpy_ HAXPY
#define daxpy_ SAXPY
#define sgemm_ HGEMM
#define dgemm_ SGEMM
#define foo foofn
#else
#define foo &fooc
#ifdef NO_UNDERSCORE
#define sgemv_ sgemv
#define dgemv_ dgemv
#define saxpy_ saxpy
#define daxpy_ daxpy
#define sgemm_ sgemm
#define dgemm_ dgemm
#endif
#endif

/* DEFAULT SETTINGS */

#ifdef DEBUG
#define CACHE_MIN_BITS  (8) /* 1<<bits is lower bound */
#define CACHE_MAX_BITS  (8) /* 1<<bits is upper bound */
#define ITERATIONS      (1)
#define LOWEST_ITERATION_COUNT (1)
#define RESOLUTION      (1)
#define REPEAT_COUNT    (1)
#else
#define CACHE_MIN_BITS  (8) /* 1<<bits is lower bound */
#define CACHE_MAX_BITS  (24) /* 1<<bits is upper bound */
#define ITERATIONS      (100000)
#define LOWEST_ITERATION_COUNT (1)
#define RESOLUTION      (2)
#define REPEAT_COUNT    (2)
#endif
#define VECTORVECTOR    (1<<0)
#define VECTORMATRIX    (1<<1)
#define MATRIXMATRIX    (1<<2)
#define NOCALIBRATE     (1<<3)
#define SINGLEPRECISION (1<<4)
#define CONSTANTITERATIONS (1<<5)
#define OPSPSECPMREF (1<<6)
#define HOLDLDA (1<<7)
#define REPORTDIMS (1<<8)

/* MACROS */

#define TIMER_START     gettimeofday(&tv1, (struct timezone*)0)
#define TIMER_STOP      gettimeofday(&tv2, (struct timezone*)0)
#define TIMER_ELAPSED_US  (double)(((tv2.tv_usec-tv1.tv_usec)+((tv2.tv_sec-tv1.tv_sec)*1000000)))

#define ARRAY(x,i,j,m) *(x+((j)*m)+(i))

#ifdef REGISTER
#undef REGISTER
#define REGISTER register
#else
#define REGISTER
#endif

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

#define CACHE_MIN (1<<CACHE_MIN_BITS)
#define CACHE_MAX (1<<CACHE_MAX_BITS)
#define DOUBLE_NS(x) ((double)x)*1.0e9

/* EXTERNALS */

extern char *optarg;

/* INTERNALS */

static float falpha = 1.0, fbeta = 1.0;
static double dalpha = 1.0, dbeta = 1.0;
static struct timeval tv1, tv2;
static int type = 0;
static char fooc = 'n';
static int stride = 1;
static int *sizes, timeslots, logmemsize = CACHE_MAX_BITS, memsize = 1<<CACHE_MAX_BITS, repeat_count = -1, resolution = RESOLUTION;
static double *t, *bws, *percents, *opssmref;

#ifdef T3E
static _fcd foofn;
#endif

/* Entry points to fool fortran linkers...*/
#ifdef __linux__
int MAIN__()
#endif
#if defined(__hppa) || defined(_HPUX_SOURCE)
int __main()
#endif
#if defined(__linux__) || defined(__hppa)
{
#ifdef __linux__
    /* Subroutine */ int s_stop();
    s_stop("", 0L);
#endif
    return(0);
}
#endif

void compute_stats(int i, int j, double refcnt, double datum_size, double opcnt,
		   double tmicrosec, int iterations, int size, int maxsize, int dim, int maxdim)
{
  double bytecnt = refcnt*datum_size;

  /* Compute NS per memory reference. */

  ARRAY(t,i,j,repeat_count) = (tmicrosec*1000.0)/refcnt;

  /* Compute MFlops/sec */
  
  ARRAY(opssmref,i,j,repeat_count) = opcnt / tmicrosec;

  /* Compute MB/sec */
  
  ARRAY(bws,i,j,repeat_count) = (((bytecnt)/(1024.0*1024.0)) / tmicrosec) * 1.0E6;
    
  if (i == 0)
    {
      if (j == 0)
	ARRAY(percents,i,j,repeat_count) = 100.0;
      else
	{
	  ARRAY(percents,i,j,repeat_count) = 100.0 *
	    ARRAY(t,i,j,repeat_count) / 
	    ARRAY(t,i,j-1,repeat_count);
	}
    }
  else
    {
      if (j == 0)
	{
	  ARRAY(percents,i,j,repeat_count) = 100.0 *
	    ARRAY(t,i,j,repeat_count) / 
	    ARRAY(t,i-1,repeat_count-1,repeat_count);
	}
      else
	{
	  ARRAY(percents,i,j,repeat_count) = 100.0 *
	    ARRAY(t,i,j,repeat_count) / 
	    ARRAY(t,i,j-1,repeat_count);
	}
    }

  if (isatty(1))
    {
      if (type & HOLDLDA)
	{
	  if (type & REPORTDIMS)
	    printf("%d\t%d\t%2.2f\t%d\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n", 
		   maxdim,
		   dim,
		   (double)((double)dim/(double)maxdim)*100.0,
		   iterations,
		   ARRAY(t,i,j,repeat_count),
		   ARRAY(percents,i,j,repeat_count),
		   ARRAY(opssmref,i,j,repeat_count),
		   ARRAY(bws,i,j,repeat_count));
	  else
	    printf("%d\t%d\t%2.2f\t%d\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n", 
		   maxsize,
		   size,
		   (double)((double)size/(double)maxsize)*100.0,
		   iterations,
		   ARRAY(t,i,j,repeat_count),
		   ARRAY(percents,i,j,repeat_count),
		   ARRAY(opssmref,i,j,repeat_count),
		   ARRAY(bws,i,j,repeat_count));
	}
      else
	{
	  if (type & REPORTDIMS)
	    printf("%d\t%d\t%2.2f\t%d\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n", 
		   maxdim,
		   dim,
		   (double)((double)dim/(double)maxdim)*100.0,
		   iterations,
		   ARRAY(t,i,j,repeat_count),
		   ARRAY(percents,i,j,repeat_count),
		   ARRAY(opssmref,i,j,repeat_count),
		   ARRAY(bws,i,j,repeat_count));
	  else
	    printf("%d\t%d\t%2.2f\t%d\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n", 
		   sizes[i],
		   size,
		   (double)((double)size/(double)sizes[i])*100.0,
		   iterations,
		   ARRAY(t,i,j,repeat_count),
		   ARRAY(percents,i,j,repeat_count),
		   ARRAY(opssmref,i,j,repeat_count),
		   ARRAY(bws,i,j,repeat_count));
	}
    }
  else
    {
      if (type & OPSPSECPMREF)
	{
	  if (type & REPORTDIMS)
	    printf("%d %f\n", dim, ARRAY(opssmref,i,j,repeat_count));
	  else
	    printf("%d %f\n", sizes[i], ARRAY(opssmref,i,j,repeat_count));
	}
      else
	{
	  if (type & REPORTDIMS)
	    printf("%d %f\n", dim, ARRAY(bws,i,j,repeat_count));
	  else
	    printf("%d %f\n", sizes[i], ARRAY(bws,i,j,repeat_count));
	}
    }
  sleep(1);
}

int compute_axpy_dimension(int maxmem, int datasize)
{
  return(maxmem/(2*datasize));
}

int compute_gemv_dimension(int maxmem, int datasize)
{
  double rt;

  rt = (4.0*(double)maxmem/(double)datasize) + 4.0;
  rt = sqrt(rt) - 2.0;
  rt = rt / 2;

  /* DBG(printf("DIM is now %d\n",dim));
  while (((dim*dim)+dim+dim) < maxmem/datasize)
    {
      dim++;
      DBG(printf("DIM is now %d\n",dim));
    }
  while (((dim*dim)+dim+dim) > maxmem/datasize)
    {
      dim--;
      DBG(printf("DIM is now %d\n",dim));
    } 
  printf("%f %d\n",rt,(int)rt); */
  return((int)rt);
}

int compute_gemm_dimension(int maxmem, int datasize)
{
  int dim;

  dim = (int)sqrt((double)(maxmem/(datasize*3)));
  
  return(dim);
}


INLINE void do_saxpy(float *x, float *y, int iterations, int *limit)
{
  REGISTER int i = 0;
  extern int saxpy_();

  for (;i<iterations;i++)
    {
      saxpy_(limit,&falpha,x,&stride,y,&stride);
    }
}

INLINE void do_daxpy(double *x, double *y, int iterations, int *limit)
{
  REGISTER int i = 0;
  extern int daxpy_();

  for (;i<iterations;i++)
    {
      daxpy_(limit,&dalpha,x,&stride,y,&stride);
    }
}

INLINE void do_sgemv(float *a, float *x, float *y, int iterations, int *limit, int *lda)
{
  REGISTER int i = 0;
  extern int sgemv_();

  for (;i<iterations;i++)
    {
      sgemv_(foo,limit,limit,&falpha,a,lda,x,&stride,&fbeta,y,&stride);
    }
}

INLINE void do_dgemv(double *a, double *x, double *y, int iterations, int *limit, int *lda)
{
  REGISTER int i = 0;
  extern int dgemv_();

  for (;i<iterations;i++)
    {
      dgemv_(foo,limit,limit,&dalpha,a,lda,x,&stride,&dbeta,y,&stride);
    }
}

INLINE void do_sgemm(float *a, float *b, float *c, int iterations, int *limit, int *lda)
{
  REGISTER int i = 0;
  extern int sgemm_();

  for (;i<iterations;i++)
    {
      sgemm_(foo,foo,limit,limit,limit,&falpha,a,lda,b,lda,&fbeta,c,lda);
    }
}

INLINE void do_dgemm(double *a, double *b, double *c, int iterations, int *limit, int *lda)
{
  REGISTER int i = 0;
  extern int dgemm_();

  for (;i<iterations;i++)
    {
      dgemm_(foo,foo,limit,limit,limit,&dalpha,a,lda,b,lda,&dbeta,c,lda);
    }
}

void flushall(int maxmem)
{
  static char *buffer = NULL;
  static char val = 1;

  if (maxmem)
    {
      if (buffer == NULL)
	assert(buffer = (char *)malloc(maxmem*sizeof(char)));
      memset(buffer,val++,maxmem);
    }
  else
    free(buffer);
}

void initall_flt(REGISTER float *x, REGISTER int number)
{
  REGISTER int i;

  for (i=0;i<number;i++)
    {
      x[i] = 1.0;
    }
}

void initall_dbl(REGISTER double *x, REGISTER int number)
{
  REGISTER int i;

  for (i=0;i<number;i++)
    {
      x[i] = 1.0;
    }
}

int usage(int argc, char **argv, int *iterations, int *real_iterations)
{
  int c;
  int errflg = 0;

  *iterations = *real_iterations = ITERATIONS;

  while ((c = getopt(argc, argv, "hm:x:e:i:b:ovatscld")) != -1)
    switch (c) {
    case 'h':
      errflg++;
      break;
    case 'm':
      if (((logmemsize = atoi(optarg)) < 0) || (logmemsize <= CACHE_MIN_BITS))
	errflg++;
      memsize = 1<<logmemsize;
      break;
    case 'v':
      if (type & (VECTORVECTOR|VECTORMATRIX|MATRIXMATRIX))
	if ((isatty(1)) == 0)
	  errflg++;
      type |= VECTORVECTOR;
      break; 
    case 'a':
      if (type & (VECTORVECTOR|VECTORMATRIX|MATRIXMATRIX))
	if ((isatty(1)) == 0)
	  errflg++;
      type |= VECTORMATRIX;
      break; 
    case 't':
      if (type & (VECTORVECTOR|VECTORMATRIX|MATRIXMATRIX))
	if ((isatty(1)) == 0)
	  errflg++;
      type |= MATRIXMATRIX;
      break; 
    case 's':
#ifdef T3E      
      errflg++;
#else
      type |= SINGLEPRECISION;
#endif
      break; 
    case 'o':
      type |= OPSPSECPMREF;
      break; 
    case 'd':
      type |= REPORTDIMS;
      break; 
    case 'l':
      type |= HOLDLDA;
      break; 
    case 'x':
      if ((resolution = atoi(optarg)) < 0)
	errflg++;
      resolution++; /* Correct for my usage */
      break; 
    case 'e':
      if ((repeat_count = atoi(optarg)) < 0)
	errflg++;
      break; 
    case 'i':
      if ((*iterations = *real_iterations = atoi(optarg)) < 0)
	errflg++;
      break; 
    case 'c':
      type |= CONSTANTITERATIONS;
      break; 
    case '?':
      errflg++;
      break; }

  if (repeat_count == -1)
    {
      if (isatty(1))
	repeat_count = REPEAT_COUNT;
      else
	repeat_count = 1;
    }

  if ((type & (VECTORVECTOR|VECTORMATRIX|MATRIXMATRIX)) == 0)
    {
      if (isatty(1)) 
	{
	  type |= VECTORVECTOR|VECTORMATRIX|MATRIXMATRIX;
	}
      else
	type |= MATRIXMATRIX;
    }

  if (errflg) 
    {
      fprintf(stderr, "Usage: %s [-vatsco -x # -m # -e # -i #]\n",argv[0]);
      fprintf(stderr, "\t -v AXPY dot product benchmark\n");
      fprintf(stderr, "\t -a GEMV matrix-vector multiply benchmark\n");
      fprintf(stderr, "\t -t GEMM matrix-matrix multiply benchmark\n");
#ifndef T3E
      fprintf(stderr, "\t -s Use single precision floating point data\n");
#endif
      fprintf(stderr, "\t -c Use constant number of iterations\n");
      fprintf(stderr, "\t -o Report Mflops/sec instead of MB/sec\n");
      fprintf(stderr, "\t -e Repeat count per problem size\n");
      fprintf(stderr, "\t -l Hold LDA and loop over sizes of square submatrices\n");
      fprintf(stderr, "\t -i Maximum iteration count\n");
      fprintf(stderr, "\t -x Number of measurements between powers of 2.\n");
      fprintf(stderr, "\t -m Specify the log2(maximum problem size) in bytes\n");
      fprintf(stderr, "\t -d Report/use dimension statistics instead of bytes\n");

      fprintf(stderr, "Default datatype   : %s, %d bytes\n","double",(int)sizeof(double));
      fprintf(stderr, "Default datatype   : %s, %d bytes\n","float",(int)sizeof(float));

      fprintf(stderr, "Defaults if to tty : -vat -x%d -m%d -e%d -i%d\n",
	      RESOLUTION-1,CACHE_MAX_BITS,REPEAT_COUNT,ITERATIONS);
      fprintf(stderr, "Defaults if to file: -t   -x%d -m%d -e1 -i%d\n",
	      RESOLUTION-1,CACHE_MAX_BITS,ITERATIONS);
      exit(1);
    }

  timeslots = resolution*(logmemsize - CACHE_MIN_BITS) + 1;

  DBG(printf("%d %d %d\n",logmemsize,memsize,timeslots)); 

  return(type);
}

void initialize_arrays(int type)
{
  int i,j;
  
  assert(t = (double *)malloc(timeslots*repeat_count*(int)sizeof(double)));
  memset(t,0x00,(timeslots*repeat_count*(int)sizeof(double)));
  assert(bws = (double *)malloc(timeslots*repeat_count*(int)sizeof(double)));
  memset(bws,0x00,(timeslots*repeat_count*(int)sizeof(double)));
  assert(opssmref = (double *)malloc(timeslots*repeat_count*(int)sizeof(double)));
  memset(opssmref,0x00,(timeslots*repeat_count*(int)sizeof(double)));
  assert(percents = (double *)malloc(timeslots*repeat_count*(int)sizeof(double)));
  memset(percents,0x00,(timeslots*repeat_count*(int)sizeof(double)));

  assert(sizes = (int *)malloc(timeslots*(int)sizeof(int)));
  for (j=0; j<timeslots; j+=resolution)
    {
      (sizes)[j] = 1<<(CACHE_MIN_BITS+j/resolution);
      DBG(printf("POW: %d %d\n",j,(sizes)[j]));

      for (i=1;i<resolution;i++)
	{
	  if (j+i < timeslots)
	    {
	      (sizes)[j+i] = (sizes)[j] + i*((sizes)[j]/resolution);
	      DBG(printf("SUB: %d %d\n",j+i,(sizes)[j+i]));
	    }
	}
    }
}

void print_header(int type, char *string)
{
  if (isatty(1))
    {
      printf("\n\t\t%s%s %s\n\n",
	     (type & NOCALIBRATE) ? "Uncalibrated " : "\t",
	     (type & SINGLEPRECISION) ? "float" : "double",string);
      if (type & REPORTDIMS)
	printf("BDims\tPDims\t%% Fit\tIters\tNansec\t%% Chng\tMFlopsS\tMB/sec\n");
      else
	printf("BSize\tPSize\t%% Fit\tIters\tNansec\t%% Chng\tMFlopsS\tMB/sec\n");
      printf("-----\t-----\t-----\t------\t------\t------\t-------\t------\n");
    }
}

int main(int argc, char **argv) 
{
  int limit, i, j, dim, iterations, real_iterations;
  double tmicrosec = 0.0, refcnt= 0.0, prefcnt = 0.0, opcnt = 0.0;

  type = usage(argc, argv, &iterations, &real_iterations);
  iterations += (iterations % 4);
  initialize_arrays(type);
  setbuf(stdout,NULL);
#ifdef T3E
  foofn = _cptofcd(&fooc, sizeof(fooc));
#endif

  /* Measure cache */

   if (type & VECTORVECTOR)
    {

      print_header(type, "AXPY Cache Test");

      if (type & SINGLEPRECISION)
	{
	  float *sx, *sy;
	  
	  dim = compute_axpy_dimension(memsize,(int)sizeof(float));
	  DBG(fprintf(stderr,"Max dimension in bytes is %d\n",dim*(int)sizeof(float)));
	  DBG(fprintf(stderr,"Total bytes used %d out of %d\n",(dim+dim)*(int)sizeof(float),memsize));

	  assert(sx = (float *)malloc(dim*(int)sizeof(float)));
	  assert(sy = (float *)malloc(dim*(int)sizeof(float)));
	  memset(sx,0x00,dim*(int)sizeof(float));
	  memset(sy,0x00,dim*(int)sizeof(float));

	  initall_flt(sx,dim);
	  initall_flt(sy,dim);

	  for (i = 0; i < timeslots; i++) 
	    {
	      limit = compute_axpy_dimension(sizes[i],(int)sizeof(float));
	      DBG(fprintf(stderr,"Max dimension in floats is %d\n",limit));
	      DBG(fprintf(stderr,"Cache size is %d bytes, %d floats\n",
			  sizes[i],sizes[i]/(int)sizeof(float)));

	      if ((type & CONSTANTITERATIONS) == 0)
		{
		  real_iterations = (i == 0 ? iterations : (int)(prefcnt/(3.0*(double)limit)));
		  if (real_iterations < LOWEST_ITERATION_COUNT)
		    real_iterations = LOWEST_ITERATION_COUNT;
		}

	      refcnt = (double)real_iterations*3.0*(double)limit;
	      opcnt = (double)real_iterations*2.0*(double)limit;
	      DBG(printf("refcnt now %f, was %f, doing %d iterations\n",refcnt,prefcnt,real_iterations));
	      prefcnt = refcnt;

	      for (j = 0; j < repeat_count; j++)
		{
		  flushall(memsize);

		  TIMER_START;
		  do_saxpy(sx,sy,real_iterations,&limit);
		  TIMER_STOP;

		  tmicrosec = TIMER_ELAPSED_US;
		  
		  DBG(fprintf(stderr,"R: %f us. total\n",tmicrosec));	  
		  DBG(fprintf(stderr,"R: %f ns. per ref\n",(tmicrosec*1000.0) / refcnt));
		  compute_stats(i,j,refcnt,(double)sizeof(float),opcnt,
				tmicrosec,real_iterations,
				(limit+limit)*(int)sizeof(float),
				(limit+limit)*(int)sizeof(float),
				limit,limit);
		}
	    }
	  free(sx);
	  free(sy);
	}
      else
	{
	  double *dx, *dy;

	  dim = compute_axpy_dimension(memsize,(int)sizeof(double));
	  DBG(fprintf(stderr,"Max dimension in bytes is %d\n",dim*(int)sizeof(double)));
	  DBG(fprintf(stderr,"Total bytes used %d out of %d\n",(dim+dim)*(int)sizeof(double),memsize));

	  assert(dx = (double *)malloc(dim*(int)sizeof(double)));
	  assert(dy = (double *)malloc(dim*(int)sizeof(double)));
	  memset(dx,0x00,dim*(int)sizeof(double));
	  memset(dy,0x00,dim*(int)sizeof(double));

	  initall_dbl(dx,dim);
	  initall_dbl(dy,dim);

	  for (i = 0; i < timeslots; i++) 
	    {
	      limit = compute_axpy_dimension(sizes[i],(int)sizeof(double));
	      DBG(fprintf(stderr,"Max dimension in doubles is %d\n",limit));
	      DBG(fprintf(stderr,"Cache size is %d bytes, %d doubles\n",
			  sizes[i],sizes[i]/(int)sizeof(double)));

	      if ((type & CONSTANTITERATIONS) == 0)
		{
		  real_iterations = (i == 0 ? iterations : (int)(prefcnt/(3.0*(double)limit)));
		  if (real_iterations < LOWEST_ITERATION_COUNT)
		    real_iterations = LOWEST_ITERATION_COUNT;
		}

	      refcnt = (double)real_iterations*3.0*(double)limit;
	      opcnt = (double)real_iterations*2.0*(double)limit;
	      DBG(printf("refcnt now %f, was %f, doing %d iterations\n",refcnt,prefcnt,real_iterations));
	      prefcnt = refcnt;

	      for (j = 0; j < repeat_count; j++)
		{
		  flushall(memsize);
		  
		  TIMER_START;
		  do_daxpy(dx,dy,real_iterations,&limit);
		  TIMER_STOP;

		  tmicrosec = TIMER_ELAPSED_US;
		  
		  DBG(fprintf(stderr,"R: %f us. total\n",tmicrosec));	  
		  DBG(fprintf(stderr,"R: %f ns. per ref\n",(tmicrosec*1000.0) / refcnt));
		  compute_stats(i,j,refcnt,(double)sizeof(double),opcnt,
				tmicrosec,real_iterations,
				(limit+limit)*(int)sizeof(double),
				(limit+limit)*(int)sizeof(double),
				limit,limit);
		}
	    }
	  free(dx);
	  free(dy);
	}
    }

   if (type & VECTORMATRIX)
    {
      int lda;

      print_header(type, "GEMV Cache Test");

      if (type & SINGLEPRECISION)
	{
	  float *sa, *sx, *sy;
      
	  dim = compute_gemv_dimension(memsize,(int)sizeof(float));
	  DBG(fprintf(stderr,"Max dimension in bytes is %d\n",dim*(int)sizeof(float)));
	  DBG(fprintf(stderr,"Total bytes used %d out of %d\n",
		      ((dim*dim)+dim+dim)*(int)sizeof(float),memsize));

	  assert(sa = (float *)malloc(dim*dim*(int)sizeof(float)));
	  assert(sx = (float *)malloc(dim*(int)sizeof(float)));
	  assert(sy = (float *)malloc(dim*(int)sizeof(float)));
	  memset(sa,0x00,dim*dim*(int)sizeof(float));
	  memset(sx,0x00,dim*(int)sizeof(float));
	  memset(sy,0x00,dim*(int)sizeof(float));

	  initall_flt(sa,dim*dim);
	  initall_flt(sx,dim);
	  initall_flt(sy,dim);

	  for (i = 0; i < timeslots; i++) 
	    {
	      limit = compute_gemv_dimension(sizes[i],(int)sizeof(float));
	      DBG(fprintf(stderr,"Max dimension in floats is %d\n",limit));
	      DBG(fprintf(stderr,"Cache size is %d bytes, %d floats\n",
			  sizes[i],sizes[i]/(int)sizeof(float)));
	      
	      if ((type & CONSTANTITERATIONS) == 0)
		{
		  real_iterations = (i == 0 ? iterations : (int)(prefcnt/(2*limit*limit + 2*limit)));
		  if (real_iterations < LOWEST_ITERATION_COUNT)
		    real_iterations = LOWEST_ITERATION_COUNT;
		}

	      refcnt = (double)real_iterations*(2*limit*limit + 2*limit);
	      opcnt = (double)real_iterations*(2*limit*limit + 2*limit);
	      DBG(printf("refcnt now %f, was %f, doing %d iterations\n",refcnt,prefcnt,real_iterations));
	      prefcnt = refcnt;
	      lda = ((type & HOLDLDA) ? dim : limit);
	      
	      for (j = 0; j < repeat_count; j++)
		{
		  flushall(memsize);

		  TIMER_START;
		  do_sgemv(sa,sx,sy,real_iterations,&limit,&lda);
		  TIMER_STOP;

		  tmicrosec = TIMER_ELAPSED_US;
		  
		  DBG(fprintf(stderr,"R: %f us. total\n",tmicrosec));	  
		  DBG(fprintf(stderr,"R: %f ns. per ref\n",(tmicrosec*1000.0) / refcnt));
		  compute_stats(i,j,refcnt,(double)sizeof(float),opcnt,
				tmicrosec,real_iterations,
				((limit*limit)+limit+limit)*(int)sizeof(float),
				((lda*lda)+lda+lda)*(int)sizeof(float),
				limit,lda);
		}
	    }
	  free(sa);
	  free(sx);
	  free(sy);
	}
      else
	{
	  double *da, *dx, *dy;

	  dim = compute_gemv_dimension(memsize,(int)sizeof(double));
	  DBG(fprintf(stderr,"Max dimension in bytes is %d\n",dim*(int)sizeof(double)));
	  DBG(fprintf(stderr,"Total bytes used %d out of %d\n",
		      ((dim*dim)+dim+dim)*(int)sizeof(double),memsize));

	  assert(da = (double *)malloc(dim*dim*(int)sizeof(double)));
	  assert(dx = (double *)malloc(dim*(int)sizeof(double)));
	  assert(dy = (double *)malloc(dim*(int)sizeof(double)));
	  memset(da,0x00,dim*dim*(int)sizeof(double));
	  memset(dx,0x00,dim*(int)sizeof(double));
	  memset(dy,0x00,dim*(int)sizeof(double));

	  initall_dbl(da,dim*dim);
	  initall_dbl(dx,dim);
	  initall_dbl(dy,dim);

	  for (i = 0; i < timeslots; i++) 
	    {
	      limit = compute_gemv_dimension(sizes[i],(int)sizeof(double));
	      DBG(fprintf(stderr,"Max dimension in doubles is %d\n",limit));
	      DBG(fprintf(stderr,"Cache size is %d bytes, %d doubles\n",
			  sizes[i],sizes[i]/(int)sizeof(double)));

	      if ((type & CONSTANTITERATIONS) == 0)
		{
		  real_iterations = (i == 0 ? iterations : (int)(prefcnt/(2*limit*limit + 2*limit)));
		  if (real_iterations < LOWEST_ITERATION_COUNT)
		    real_iterations = LOWEST_ITERATION_COUNT;
		}
	      
	      refcnt = (double)real_iterations*(2*limit*limit + 2*limit);
	      opcnt = (double)real_iterations*(2*limit*limit + 2*limit);
	      DBG(printf("refcnt now %f, was %f, doing %d iterations\n",refcnt,prefcnt,real_iterations));
	      prefcnt = refcnt;
	      lda = ((type & HOLDLDA) ? dim : limit);

	      for (j = 0; j < repeat_count; j++)
		{
		  flushall(memsize);

		  TIMER_START;
		  do_dgemv(da,dx,dy,real_iterations,&limit,&lda);
		  TIMER_STOP;

		  tmicrosec = TIMER_ELAPSED_US;
		  
		  DBG(fprintf(stderr,"R: %f us. total\n",tmicrosec));	  
		  DBG(fprintf(stderr,"R: %f ns. per ref\n",(tmicrosec*1000.0) / refcnt));
		  compute_stats(i,j,refcnt,(double)sizeof(double),opcnt,
				tmicrosec,real_iterations,
				((limit*limit)+limit+limit)*(int)sizeof(double),
				((lda*lda)+lda+lda)*(int)sizeof(double),
				limit,lda);
		}
	    }
	  free(da);
	  free(dx);
	  free(dy);
	}
    }

   if (type & MATRIXMATRIX)
    {
      int lda;

      print_header(type, "GEMM Cache Test");

      if (type & SINGLEPRECISION)
	{
	  float *sa, *sb, *sc;

	  dim = compute_gemm_dimension(memsize,(int)sizeof(float));
	  DBG(fprintf(stderr,"Max dimension in bytes is %d\n",dim*(int)sizeof(float)));
	  DBG(fprintf(stderr,"Total bytes used %d out of %d\n",(3*dim*dim)*(int)sizeof(float),memsize));

	  assert(sa = (float *)malloc(dim*dim*(int)sizeof(float)));
	  assert(sb = (float *)malloc(dim*dim*(int)sizeof(float)));
	  assert(sc = (float *)malloc(dim*dim*(int)sizeof(float)));
	  memset(sa,0x00,dim*dim*(int)sizeof(float));
	  memset(sb,0x00,dim*dim*(int)sizeof(float));
	  memset(sc,0x00,dim*dim*(int)sizeof(float));

	  initall_flt(sa,dim*dim);
	  initall_flt(sa,dim*dim);
	  initall_flt(sa,dim*dim);

	  for (i = 0; i < timeslots; i++) 
	    {
	      limit = compute_gemm_dimension(sizes[i],(int)sizeof(float));
	      DBG(fprintf(stderr,"Max dimension in floats is %d\n",limit));
	      DBG(fprintf(stderr,"Cache size is %d bytes, %d floats\n",
			  sizes[i],sizes[i]/(int)sizeof(float)));
	      
	      if ((type & CONSTANTITERATIONS) == 0)
		{
		  real_iterations = (i == 0 ? iterations : (int)(prefcnt/(2*limit*limit*limit+2*limit*limit)));
		  /* For each of n*n elements, we have 2*n + 1 reads, 1 write */
		  if (real_iterations < LOWEST_ITERATION_COUNT)
		    real_iterations = LOWEST_ITERATION_COUNT;
		}

	      refcnt = (double)real_iterations*(2*limit*limit*limit+2*limit*limit);
	      opcnt = (double)real_iterations*(2*limit*limit*limit+2*limit*limit);
	      DBG(printf("refcnt now %f, was %f, doing %d iterations\n",refcnt,prefcnt,real_iterations));
	      prefcnt = refcnt;
	      lda = ((type & HOLDLDA) ? dim : limit);

	      for (j = 0; j < repeat_count; j++)
		{
		  flushall(memsize);

		  TIMER_START;
		  do_sgemm(sa,sb,sc,real_iterations,&limit,&lda);
		  TIMER_STOP;

		  tmicrosec = TIMER_ELAPSED_US;

		  DBG(fprintf(stderr,"R: %f us. total\n",tmicrosec));	  
		  DBG(fprintf(stderr,"R: %f ns. per ref\n",(tmicrosec*1000.0) / refcnt));
		  compute_stats(i,j,refcnt,(double)sizeof(float),opcnt,
				tmicrosec,real_iterations,
				(3*limit*limit)*(int)sizeof(float),
				(3*lda*lda)*(int)sizeof(float),
				limit,lda);
		}
	    }

	  free(sa);
	  free(sb);
	  free(sc);
	}
      else
	{
	  double *da, *db, *dc;

	  dim = compute_gemm_dimension(memsize,(int)sizeof(double));
	  DBG(fprintf(stderr,"Max dimension in bytes is %d\n",dim*(int)sizeof(double)));
	  DBG(fprintf(stderr,"Total bytes used %d out of %d\n",(3*dim*dim)*(int)sizeof(double),memsize));

	  assert(da = (double *)malloc(dim*dim*(int)sizeof(double)));
	  assert(db = (double *)malloc(dim*dim*(int)sizeof(double)));
	  assert(dc = (double *)malloc(dim*dim*(int)sizeof(double)));
	  memset(da,0x00,dim*dim*(int)sizeof(double));
	  memset(db,0x00,dim*dim*(int)sizeof(double));
	  memset(dc,0x00,dim*dim*(int)sizeof(double));

	  initall_dbl(da,dim*dim);
	  initall_dbl(db,dim*dim);
	  initall_dbl(dc,dim*dim);

	  for (i = 0; i < timeslots; i++) 
	    {
	      limit = compute_gemm_dimension(sizes[i],(int)sizeof(double));
	      DBG(fprintf(stderr,"Max dimension in doubles is %d\n",limit));
	      DBG(fprintf(stderr,"Cache size is %d bytes, %d doubles\n",
			  sizes[i],sizes[i]/(int)sizeof(double)))
	      
	      if ((type & CONSTANTITERATIONS) == 0)
		{
		  real_iterations = (i == 0 ? iterations : (int)(prefcnt/((limit*limit)*((2*limit)+2))));
		  /* For each of n*n elements, we have 2*n + 1 reads, 1 write */
		  if (real_iterations < LOWEST_ITERATION_COUNT)
		    real_iterations = LOWEST_ITERATION_COUNT;
		}

	      refcnt = (double)real_iterations*(limit*limit)*((2.0*limit)+2.0);
	      opcnt = (double)real_iterations*(2.0*limit*limit*limit+2.0*limit*limit);
	      DBG(printf("refcnt now %f, was %f, doing %d iterations\n",refcnt,prefcnt,real_iterations));
	      prefcnt = refcnt;
	      lda = ((type & HOLDLDA) ? dim : limit);

	      for (j = 0; j < repeat_count; j++)
		{
		  flushall(memsize);

		  TIMER_START;
		  do_dgemm(da,db,dc,real_iterations,&limit,&lda);
		  TIMER_STOP;

		  tmicrosec = TIMER_ELAPSED_US;

		  DBG(fprintf(stderr,"R: %f us. total\n",tmicrosec));	  
		  DBG(fprintf(stderr,"R: %f ns. per ref\n",(tmicrosec*1000.0) / refcnt));
		  compute_stats(i,j,refcnt,(double)sizeof(double),opcnt,
				tmicrosec,real_iterations,
				(3*limit*limit)*(int)sizeof(double),
				(3*lda*lda)*(int)sizeof(double),
				limit,lda);
		}
	    }

	  free(da);
	  free(db);
	  free(dc);
	}
    }

   flushall(0);
   exit(0);
}


