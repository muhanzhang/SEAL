/* sorttemplates.c version 1.1, Oct 9, 2013.
 * Author: Brendan McKay; bdm@cs.anu.edu.au
 *
 * This file contains templates for creating in-place sorting procedures
 * for different data types.  It cannot be compiled separately but
 * should be #included after defining a few preprocessor variables.
 *   SORT_OF_SORT, SORT_NAME and SORT_TYPE1 are required, and
 *   SORT_TYPE2 is needed for SORT_OF_SORT > 1.
 * 
 *   SORT_OF_SORT = 1: Creates a procedure
 *         static void SORT_NAME(SORT_TYPE1 *x, int n)
 *     which permutes x[0..n-1] so that x[0] <= ... <= x[n-1].
 *   SORT_OF_SORT = 2: Creates a procedure
 *         static void SORT_NAME(SORT_TYPE1 *x, SORT_TYPE2 *y, int n)
 *     which permutes x[0..n-1] so that x[0] <= ... <= x[n-1]
 *     and also permutes y[0..n-1] by the same permutation.
 *   SORT_OF_SORT = 3: Creates a procedure
 *         static void SORT_NAME(SORT_TYPE1 *x, SORT_TYPE2 *y, int n)
 *     which permutes x[0..n-1] so that y[x[0]] <= ... <= y[x[n-1]].
 *     The array y[] is not changed.
 * 
 *   SORT_NAME = the name of the procedure to be created
 *
 *   SORT_TYPE1 = type of the first or only array (no default)
 *   	 This can be any numeric type for SORT_OF_SORT=1,2, but
 *   	 should be an integer type for SORT_OF_SORT=3.
 *   SORT_TYPE2 = type of the second array if needed (no default)
 *       This can be any assignable type (including a structure) for
 *       SORT_OF_SORT=2, but must be a numeric type for SORT_OF_SORT=3.
 *
 *   SORT_MINPARTITION = least number of elements for using quicksort
 *           partitioning, otherwise insertion sort is used (default "11") 
 *   SORT_MINMEDIAN9 = least number of elements for using the median of 3
 *           medians of 3 for partitioning (default "320")
 *   SORT_FUNCTYPE = type of sort function (default "static void")
 *
 *   This file can be included any number of times provided the value
 *   of SORT_NAME is different each time.
 */
 
#define SORT_MEDIAN_OF_3(a,b,c) \
 ((a) <= (b) ? ((b) <= (c) ? (b) : (c) <= (a) ? (a) : (c)) \
             : ((a) <= (c) ? (a) : (c) <= (b) ? (b) : (c)))

#if !defined(SORT_OF_SORT) || !defined(SORT_NAME)
 #error Either SORT_OF_SORT or SORT_NAME is undefined
#endif

#if (SORT_OF_SORT < 1) || (SORT_OF_SORT > 3)
 #error Unknown value of SORT_OF_SORT
#endif

#ifndef SORT_TYPE1
 #error "SORT_TYPE1 must be defined before including sorttemplates.c"
#endif

#ifndef SORT_MINPARTITION
#define SORT_MINPARTITION 11
#endif

#ifndef SORT_MINMEDIAN9
#define SORT_MINMEDIAN9 320
#endif

#ifndef SORT_FUNCTYPE
#define SORT_FUNCTYPE static void
#endif

#define SORT_SWAP1(x,y) tmp1 = x; x = y; y = tmp1;
#define SORT_SWAP2(x,y) tmp2 = x; x = y; y = tmp2;

/*******************************************************************/

#if SORT_OF_SORT == 1
SORT_FUNCTYPE
SORT_NAME(SORT_TYPE1 *x, int n)
{
    int i,j;
    int a,d,ba,dc,s,nn;
    SORT_TYPE1 tmp1,v,v1,v2,v3;
    SORT_TYPE1 *x0,*xa,*xb,*xc,*xd,*xh,*xl;
    struct {SORT_TYPE1 *addr; int len;} stack[40];
    int top;

    top = 0;
    if (n > 1)
    {
        stack[top].addr = x;
        stack[top].len = n;
        ++top;
    }

    while (top > 0)
    {
	--top;
	x0 = stack[top].addr;
	nn = stack[top].len;

        if (nn < SORT_MINPARTITION)
        {
            for (i = 1; i < nn; ++i)
            {
                tmp1 = x0[i];
                for (j = i; x0[j-1] > tmp1; )
                {
                    x0[j] = x0[j-1];
                    if (--j == 0) break;
                }
                x0[j] = tmp1;
            }
	    continue;
        }

        if (nn < SORT_MINMEDIAN9)
            v = SORT_MEDIAN_OF_3(x0[0],x0[nn/2],x0[nn-1]);
        else
        {
	    v1 = SORT_MEDIAN_OF_3(x0[0],x0[1],x0[2]);
	    v2 = SORT_MEDIAN_OF_3(x0[nn/2-1],x0[nn/2],x0[nn/2+1]);
	    v3 = SORT_MEDIAN_OF_3(x0[nn-3],x0[nn-2],x0[nn-1]);
	    v = SORT_MEDIAN_OF_3(v1,v2,v3);
        }

        xa = xb = x0;  xc = xd = x0+(nn-1);
        for (;;)
        {
	    while (xb <= xc && *xb <= v)
	    {
	        if (*xb == v)
	        {
		    *xb = *xa; *xa = v; ++xa;
	        }
	        ++xb;
	    }
	    while (xc >= xb && *xc >= v)
	    {
	        if (*xc == v)
	        {
		    *xc = *xd; *xd = v; --xd;
	        }
	        --xc;
	    }
	    if (xb > xc) break;
            SORT_SWAP1(*xb,*xc);
	    ++xb;
	    --xc;
        }
    
        a = xa - x0;
        ba = xb - xa;
        if (ba > a) s = a; else s = ba;
        for (xl = x0, xh = xb-s; s > 0; --s)
        {
	    *xl = *xh; *xh = v; ++xl; ++xh;
        }
        d = xd - x0;
        dc = xd - xc;
        if (dc > nn-1-d) s = nn-1-d; else s = dc;
        for (xl = xb, xh = x0 + (nn-s); s > 0; --s)
        {
	    *xh = *xl; *xl = v; ++xl; ++xh;
        }

	if (ba > dc)
	{
	    if (ba > 1)
	    {
	        stack[top].addr = x0; stack[top].len = ba; ++top;
	    }
	    if (dc > 1)
	    {
	        stack[top].addr = x0+(nn-dc); stack[top].len = dc; ++top;
	    }
	}
	else
	{
	    if (dc > 1)
	    {
	        stack[top].addr = x0+(nn-dc); stack[top].len = dc; ++top;
	    }
	    if (ba > 1)
	    {
	        stack[top].addr = x0; stack[top].len = ba; ++top;
	    }
	}
    }
}
#endif

#if SORT_OF_SORT == 2
#ifndef SORT_TYPE2
 #error "SORT_TYPE2 must be defined before including sorttemplates.c"
#endif

SORT_FUNCTYPE
SORT_NAME(SORT_TYPE1 *x, SORT_TYPE2 *y, int n)
{
    int i,j;
    int a,d,ba,dc,s,nn;
    SORT_TYPE2 tmp2,*y0,*ya,*yb,*yc,*yd,*yl,*yh;
    SORT_TYPE1 tmp1,v,v1,v2,v3;
    SORT_TYPE1 *x0,*xa,*xb,*xc,*xd,*xh,*xl;
    struct {SORT_TYPE1 *addr; int len;} stack[40];
    int top;

    top = 0;
    if (n > 1)
    {
        stack[top].addr = x;
        stack[top].len = n;
        ++top;
    }

    while (top > 0)
    {
	--top;
	x0 = stack[top].addr;
	y0 = y + (x0-x);
	nn = stack[top].len;

        if (nn < SORT_MINPARTITION)
        {
            for (i = 1; i < nn; ++i)
            {
                tmp1 = x0[i];
                tmp2 = y0[i];
                for (j = i; x0[j-1] > tmp1; )
                {
                    x0[j] = x0[j-1];
                    y0[j] = y0[j-1];
                    if (--j == 0) break;
                }
                x0[j] = tmp1;
                y0[j] = tmp2;
            }
	    continue;
        }

        if (nn < SORT_MINMEDIAN9)
            v = SORT_MEDIAN_OF_3(x0[0],x0[nn/2],x0[nn-1]);
        else
        {
	    v1 = SORT_MEDIAN_OF_3(x0[0],x0[1],x0[2]);
	    v2 = SORT_MEDIAN_OF_3(x0[nn/2-1],x0[nn/2],x0[nn/2+1]);
	    v3 = SORT_MEDIAN_OF_3(x0[nn-3],x0[nn-2],x0[nn-1]);
	    v = SORT_MEDIAN_OF_3(v1,v2,v3);
        }

        xa = xb = x0;  xc = xd = x0+(nn-1);
        ya = yb = y0;  yc = yd = y0+(nn-1);
        for (;;)
        {
	    while (xb <= xc && *xb <= v)
	    {
	        if (*xb == v)
	        {
		    *xb = *xa; *xa = v; ++xa;
		    SORT_SWAP2(*ya,*yb); ++ya;
	        }
	        ++xb; ++yb;
	    }
	    while (xc >= xb && *xc >= v)
	    {
	        if (*xc == v)
	        {
		    *xc = *xd; *xd = v; --xd;
		    SORT_SWAP2(*yc,*yd); --yd;
	        }
	        --xc; --yc;
	    }
	    if (xb > xc) break;
            SORT_SWAP1(*xb,*xc);
            SORT_SWAP2(*yb,*yc);
	    ++xb; ++yb;
	    --xc; --yc;
        }
    
        a = xa - x0;
        ba = xb - xa;
        if (ba > a) s = a; else s = ba;
        for (xl = x0, xh = xb-s, yl = y0, yh = yb-s; s > 0; --s)
        {
	    *xl = *xh; *xh = v; ++xl; ++xh;
	    SORT_SWAP2(*yl,*yh); ++yl; ++yh;
        }
        d = xd - x0;
        dc = xd - xc;
        if (dc > nn-1-d) s = nn-1-d; else s = dc;
        for (xl = xb, xh = x0+(nn-s), yl = yb, yh = y0+(nn-s); s > 0; --s)
        {
	    *xh = *xl; *xl = v; ++xl; ++xh;
	    SORT_SWAP2(*yl,*yh); ++yl; ++yh;
        }

	if (ba > dc)
	{
	    if (ba > 1)
	    {
	        stack[top].addr = x0; stack[top].len = ba; ++top;
	    }
	    if (dc > 1)
	    {
	        stack[top].addr = x0+(nn-dc); stack[top].len = dc; ++top;
	    }
	}
	else
	{
	    if (dc > 1)
	    {
	        stack[top].addr = x0+(nn-dc); stack[top].len = dc; ++top;
	    }
	    if (ba > 1)
	    {
	        stack[top].addr = x0; stack[top].len = ba; ++top;
	    }
	}
    }
}
#endif

#if SORT_OF_SORT == 3
#ifndef SORT_TYPE2
 #error "SORT_TYPE2 must be defined before including sorttemplates.c"
#endif

SORT_FUNCTYPE
SORT_NAME(SORT_TYPE1 *x, SORT_TYPE2 *y, int n)
{
    int i,j;
    int a,d,ba,dc,s,nn;
    SORT_TYPE2 tmp2,v,v1,v2,v3;
    SORT_TYPE1 tmp1,*x0,*xa,*xb,*xc,*xd,*xh,*xl;
    struct {SORT_TYPE1 *addr; int len;} stack[40];
    int top;

    top = 0;
    if (n > 1)
    {
        stack[top].addr = x;
        stack[top].len = n;
        ++top;
    }

    while (top > 0)
    {
	--top;
	x0 = stack[top].addr;
	nn = stack[top].len;

        if (nn < SORT_MINPARTITION)
        {
            for (i = 1; i < nn; ++i)
            {
                tmp1 = x0[i];
		tmp2 = y[tmp1];
                for (j = i; y[x0[j-1]] > tmp2; )
                {
                    x0[j] = x0[j-1];
                    if (--j == 0) break;
                }
                x0[j] = tmp1;
            }
	    continue;
        }

        if (nn < SORT_MINMEDIAN9)
            v = SORT_MEDIAN_OF_3(y[x0[0]],y[x0[nn/2]],y[x0[nn-1]]);
        else
        {
	    v1 = SORT_MEDIAN_OF_3(y[x0[0]],y[x0[1]],y[x0[2]]);
	    v2 = SORT_MEDIAN_OF_3(y[x0[nn/2-1]],y[x0[nn/2]],y[x0[nn/2+1]]);
	    v3 = SORT_MEDIAN_OF_3(y[x0[nn-3]],y[x0[nn-2]],y[x0[nn-1]]);
	    v = SORT_MEDIAN_OF_3(v1,v2,v3);
        }

        xa = xb = x0;  xc = xd = x0+(nn-1);
        for (;;)
        {
	    while (xb <= xc && y[*xb] <= v)
	    {
	        if (y[*xb] == v)
	        {
		    SORT_SWAP1(*xa,*xb); ++xa;
	        }
	        ++xb;
	    }
	    while (xc >= xb && y[*xc] >= v)
	    {
	        if (y[*xc] == v)
	        {
		    SORT_SWAP1(*xc,*xd); --xd;
	        }
	        --xc;
	    }
	    if (xb > xc) break;
            SORT_SWAP1(*xb,*xc);
	    ++xb;
	    --xc;
        }
    
        a = xa - x0;
        ba = xb - xa;
        if (ba > a) s = a; else s = ba;
        for (xl = x0, xh = xb-s; s > 0; --s)
        {
	    SORT_SWAP1(*xl,*xh); ++xl; ++xh;
        }
        d = xd - x0;
        dc = xd - xc;
        if (dc > nn-1-d) s = nn-1-d; else s = dc;
        for (xl = xb, xh = x0 + (nn-s); s > 0; --s)
        {
	    SORT_SWAP1(*xl,*xh); ++xl; ++xh;
        }

	if (ba > dc)
	{
	    if (ba > 1)
	    {
	        stack[top].addr = x0; stack[top].len = ba; ++top;
	    }
	    if (dc > 1)
	    {
	        stack[top].addr = x0+(nn-dc); stack[top].len = dc; ++top;
	    }
	}
	else
	{
	    if (dc > 1)
	    {
	        stack[top].addr = x0+(nn-dc); stack[top].len = dc; ++top;
	    }
	    if (ba > 1)
	    {
	        stack[top].addr = x0; stack[top].len = ba; ++top;
	    }
	}
    }
}
#endif

#undef SORT_NAME
#undef SORT_OF_SORT
#undef SORT_TYPE1
#undef SORT_TYPE2
