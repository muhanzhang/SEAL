/* Compare times for popcount instructions */
/* Usage:  poptest K N
   - measures the time for 1000*N popcount operations on words
     with K one bits, comparing with the macro POPCOUNTMAC.
     Compile with values for WORDSIZE and POPC with POPC values:
     0 = loop, 1 = popcount, 2 = popcountl, 3 = popcountll,
                    4 = _mm_popcnt_u32, 5 = _mm_popcnt_u32 */

#define MAXN WORDSIZE

#ifndef WORDSIZE 
#define WORDSIZE 32
#endif

#ifndef POPC
#error  Need a value for POPC
#endif

#include "gtools.h"
#ifdef __INTEL_COMPILER
#include <nmmintrin.h>
#endif

#ifndef POPCOUNTMAC
#define POPCOUNTMAC POPCOUNT
#endif

#if POPC==0
#define NEWPOPC(x,c) {c = 0; while(x){++c; x &= x-1;}}
#endif
#if POPC==1
#define NEWPOPC(x,c) {c = __builtin_popcount(x);}
#endif
#if POPC==2
#define NEWPOPC(x,c) {c = __builtin_popcountl(x);}
#endif
#if POPC==3
#define NEWPOPC(x,c) {c = __builtin_popcountll(x);}
#endif
#if POPC==4
#define NEWPOPC(x,c) {c = _mm_popcnt_u32(x);}
#endif
#if POPC==5
#define NEWPOPC(x,c) {c = _mm_popcnt_u64(x);}
#endif

static setword
ransetword(int k)     /* setword with k random bits */
{
	register setword w,rb;
	register int j;

	j = 0;
	w = 0;
	for (;;)
	{
	    if (j == k) return w;

	    rb = bit[random() % WORDSIZE];
	    if ((rb & w) == 0)
	    {
		w |= rb;
		++j;
	    }
	}
}

static double
timemac(setword *x, int n, int iters, int *sump)
{
	register int c,i,j,sum;
	register setword w;
	int it;
	double t;

	sum = 0;
	for (i = 0; i < n; ++i)
	{
	    w = x[i];
	    sum += POPCOUNTMAC(w);
	}

	t = CPUTIME;

	for (it = 0; it < iters; ++it)
	{
	    for (i = 0; i < n; ++i)
            {
                w = x[i];
                sum += POPCOUNTMAC(w);
            }
	    sum ^= it;
	}	
	    
	t = CPUTIME - t;
	*sump = sum;

	return t;
}

static double
timeold(setword *x, int n, int iters, int *sump)
{
	register int c,i,j,sum;
	register setword w;
	int it;
	double t;

	sum = 0;
	for (i = 0; i < n; ++i)
	{
	    w = x[i];
	    sum += POPCOUNT(w);
	}

	t = CPUTIME;

	for (it = 0; it < iters; ++it)
	{
	    for (i = 0; i < n; ++i)
            {
                w = x[i];
                sum += POPCOUNT(w);
            }
	    sum ^= it;
	}	
	    
	t = CPUTIME - t;
	*sump = sum;

	return t;
}

static double
timenew(setword *x, int n, int iters, int *sump)
{
	register int c,i,j,sum;
	register setword w;
	int it;
	double t;

	sum = 0;
	for (i = 0; i < n; ++i)
	{
	    w = x[i];
	    NEWPOPC(w,c);
	    sum += c;
	}

	t = CPUTIME;

	for (it = 0; it < iters; ++it)
	{
	    for (i = 0; i < n; ++i)
            {
                w = x[i];
		NEWPOPC(w,c);
                sum += c;
            }
	    sum ^= it;
	}	
	    
	t = CPUTIME - t;
	*sump = sum;

	return t;
}

static double
timenull(setword *x, int n, int iters, int *sump)
{
	register int c,i,j,sum;
	register setword w;
	int it;
	double t;

	sum = 0;
	for (i = 0; i < n; ++i)
	{
	    w = x[i];
	    NEWPOPC(w,c);
	    sum += c;
	}

	t = CPUTIME;

	for (it = 0; it < iters; ++it)
	{
	    for (i = 0; i < n; ++i)
            {
                w = x[i];
		c = w;
                sum += c;
            }
	    sum ^= it;
	}	
	    
	t = CPUTIME - t;
	*sump = sum;

	return t;
}

int
main(int argc, char *argv[])
{
	int i,k,iters;
	setword x[1000];
	double tnull,told,tnew,tmac;
	int summac,sumold,sumnew,sumnull;

	printf("WORDSIZE=%d POPC=%s  ",WORDSIZE,
          POPC==0 ? "loop" :
          POPC==1 ? "popcount" : 
          POPC==2 ? "popcountl" : 
          POPC==3 ? "popcountll" : 
          POPC==4 ? "popcnt_u32" : 
          POPC==5 ? "popcnt_u64" : "undefined");
#ifdef SETWORD_INT
        printf(" setword=unsigned int ");
#endif
#ifdef SETWORD_LONG
        printf(" setword=unsigned long ");
#endif
#ifdef SETWORD_LONGLONG
        printf(" setword=unsigned long long ");
#endif
#ifdef __SSE4_2__
	printf("__SSE4_2__ ");
#endif
#ifdef __POPCNT__
	printf("__POPCNT__ ");
#endif
#ifdef __INTEL_COMPILER
	printf("__INTEL_COMPILER ");
#endif
printf("\n");

	if (argc != 3)
        {
	    fprintf(stderr,"Usage: poptest num1bits numiters\n");
	    exit(1);
        }

	k = atoi(argv[1]);
	if (k > WORDSIZE) k = WORDSIZE;
	iters = atoi(argv[2]);

	for (i = 0; i < 1000; ++i)
	    x[i] = ransetword(k);

	tnull = timenull(x,1000,iters,&sumnull);
	tmac = timemac(x,1000,iters,&summac);
	told = timeold(x,1000,iters,&sumold);
	tnew = timenew(x,1000,iters,&sumnew);

	if (summac != sumold) printf("*** sum mismatch (mac/old)\n");
	if (sumold != sumnew) printf("*** sum mismatch (old/new)\n");

	printf("macro=%3.2f compiled=%3.2f new=%3.2f\n",
		tmac-tnull,told-tnull,tnew-tnull);

	return 0;
}
