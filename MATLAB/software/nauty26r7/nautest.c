/* Test for basic nauty functions (but not nauty itself) */

#include "naututil.h"

long seed;

int
main(int argc, char *argv[])
{
	int i,j,bad;
	setword w,ww;
        int curfile;
        FILE *f;
#ifdef CPUDEFS
        CPUDEFS
#endif

	printf("NAUTYVERSION=%s  NAUTYVERSIONID=%d  HAVE_TLS=%d\n",
		NAUTYVERSION,NAUTYVERSIONID,HAVE_TLS);
	printf("MAXN=%d  MAXM=%d  WORDSIZE=%d  NAUTY_INFINITY=%d\n",
		MAXN,MAXM,WORDSIZE,NAUTY_INFINITY);
	printf("sizes: short=%d int=%d long=%d double=%d boolean=%d setword=%d\n",
		(int)sizeof(short),(int)sizeof(int),(int)sizeof(long),
		(int)sizeof(double),(int)sizeof(boolean),(int)sizeof(setword));
        printf("CLZ=%d,%d,%d  POPCNT=%d,%d,%d;%d,%d\n",
                HAVE_CLZ,HAVE_CLZL,HAVE_CLZLL,
                HAVE_POPCNT,HAVE_POPCNTL,HAVE_POPCNTLL,HAVE_MMPOP32,HAVE_MMPOP64);
        printf("LONG_LONG_COUNTERS=%d  COUNTER_FMT=%s\n",
                LONG_LONG_COUNTERS,COUNTER_FMT);

#if SIZEOF_LONGLONG > 0
	printf("sizeof(long long)=%d\n",sizeof(long long));
#endif

	printf("defined:");
#ifdef __STDC__
	printf(" __STDC__");
#endif
#ifdef BIGNAUTY
	printf(" BIGNAUTY(obsolete!)");
#endif
#ifdef SYS_UNIX
	printf(" SYS_UNIX");
#endif
#ifdef SYS_CRAY
        printf(" SYS_CRAY");
#endif
#ifdef SETWORD_SHORT
	printf(" SETWORD_SHORT");
#endif
#ifdef SETWORD_INT
	printf(" SETWORD_INT");
#endif
#ifdef SETWORD_LONG
	printf(" SETWORD_LONG");
#endif
#ifdef SETWORD_LONGLONG
	printf(" SETWORD_LONGLONG");
#endif
	printf("\n");

	bad = 0;

	if (8*sizeof(setword) != WORDSIZE)
	{
	    printf("\n ***** NOTE:  WORDSIZE mismatch *****\n\n");
	    ++bad;
	}

	for (i = 0; i < WORDSIZE; ++i)
	{
	    w = ALLMASK(i);
	    if (POPCOUNT(w) != i)
	    {
		printf("\n ***** POPCOUNT(ALLMASK) error %d *****\n\n",i);
		++bad;
	    }
	}

	for (i = 0; i < WORDSIZE; ++i)
        {
            w = BITMASK(i);
            if (POPCOUNT(w) != WORDSIZE-i-1)
            {
                printf("\n ***** POPCOUNT(BITMASK) error %d *****\n\n",i);
                ++bad;
            }
        }

	for (i = 0; i < WORDSIZE; ++i)
	    if (POPCOUNT(ALLMASK(i)) != i)
	    {
		printf("\n ***** POPCOUNT(ALLMASK) error %d *****\n\n",i);
		++bad;
	    }

	for (i = 0; i < WORDSIZE; ++i)
            if (FIRSTBIT(BITT[i]) != i)
	    {
		printf("\n ***** FIRSTBIT(BITT) error %d *****\n\n",i);
		++bad;
	    }

	if (FIRSTBIT((setword)0) != WORDSIZE)
	{
	    printf("\n ***** FIRSTBIT(0) error *****\n\n");
	    ++bad;
	}
	
	for (i = 0; i < WORDSIZE; ++i)
            if (POPCOUNT(BITT[i]) != 1)
	    {
                printf("\n ***** POPCOUNT(BITT) error %d *****\n\n",i);
		++bad;
	    }

	for (i = 0; i < WORDSIZE; ++i)
	{
	    w = 0;
	    for (j = 1; j <= WORDSIZE; ++j)
	    {
		w |= BITT[(j*97+i)%WORDSIZE];
		if (POPCOUNT(w) != j)
		{
		    printf("\n ***** POPCOUNT(w) error %d %d *****\n\n",i,j);
		    ++bad;
		}
	    }
	}

#ifdef DOPROMPT
	curfile = 0;
	printf("DOPROMPT(stdin)=%d DOPROMPT(stdout)=%d\n",
		DOPROMPT(stdin),DOPROMPT(stdout));
#else
	printf("DOPROMPT is not defined\n");
#endif

#ifdef CPUTIME
        printf("CPUTIME = %f\n",CPUTIME);
#else
        printf("CPUTIME is not defined\n");
#endif

#ifdef INITSEED
	INITSEED;
        printf("INITSEED: seed=%ld\n",seed);
#else
        printf("INITSEED is not defined\n");
#endif

#ifdef OPENOUT
        OPENOUT(f,"nautest.txt",0);
	fprintf(f,"test\n");
#else
        printf("OPENOUT is not defined\n");
#endif

	if (!bad) printf("\nNo errors found\n");
	else      printf("\nXXXXXXX %d errors found XXXXXXX\n",bad);

	exit(0);
} 
