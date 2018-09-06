/* gentourng.c  version 1.4; B D McKay, Jan 20, 2016 */

#define USAGE \
"gentourng [-cd#D#] [-ugsz] [-lq] n [res/mod] [file]"

#define HELPTEXT \
" Generate all tournaments of a specified class.\n\
\n\
      n    : the number of vertices\n\
   res/mod : only generate subset res out of subsets 0..mod-1\n\
\n\
     -c    : only write strongly-connected tournaments\n\
     -d#   : a lower bound for the minimum out-degree\n\
     -D#   : a upper bound for the maximum out-degree\n\
     -l    : canonically label output graphs\n\
\n\
     -u    : do not output any graphs, just generate and count them\n\
     -g    : use graph6 output (lower triangle)\n\
     -s    : use sparse6 output (lower triangle)\n\
     -z    : use digraph6 output\n\
     -h    : write a header (only with -g or -s)\n\
  Default output is upper triangle row-by-row in ascii\n\
\n\
     -q    : suppress auxiliary output\n\
\n\
  See program text for much more information.\n"


/*  Parameters:

             n    = the number of vertices (1..min(32,WORDSIZE))
             mod, res = a way to restrict the output to a subset.
                        All the graphs in G(n,mine..maxe) are divided into
                        disjoint classes C(0,mod),C(1,mod),...,C(mod-1,mod),
                        of very approximately equal size.
                        Only the class C(res,mod) is written.

	     file = a name for the output file (stdout if missing or "-")

	     All switches can be concatenated or separate.  However, the
             value of -d must be attached to the "d", and similarly for "x".

             -c    : only write connected graphs
             -d<int> : specify an upper bound for the maximum out-degree.
                     The value of the upper bound must be adjacent to
                     the "d".  Example: -d6
             -l    : canonically label output graphs

             -u    : do not output any graphs, just generate and count them
             -z    : use digraph6 output
	     -g    : use graph6 output
	     -s    : use sparse6 output
		For -g and -s, the lower triangle of the adjacency matrix
                is written as if it is an undirected graph.  Nauty tools
		like labelg do not know this format.  To read it you can
                read it as an undirected graph then complement the upper
		triangle.
	     -h    : for graph6 or sparse6 format, write a header too

	     -q    : suppress auxiliary output (except from -v)

Output formats.

  The output format is determined by the mutually exclusive switches
  -u, -z, -g and -s.  The default is ascii format.

  -u suppresses output of graphs completely.
  -z uses digraph6 output.

  -s and -g specify sparse6 and graph6 format, defined elsewhere.
  In this case a header is also written if -h is present.

OUTPROC feature.

   By defining the C preprocessor variable OUTPROC at compile time
   (for Unix the syntax is  -DOUTPROC=procname  on the cc command),
   gentourng can be made to call a procedure of your manufacture with each
   output graph instead of writing anything. Your procedure needs
   to have type void and the argument list (FILE *f, graph *g, int n).
   f is a stream open for writing, g is the graph in nauty format,
   and n is the number of vertices. Your procedure can be in a
   separate file so long as it is linked with gentourng. The global
   variables nooutput, and canonise (all type boolean) can be
   used to test for the presence of the flags -u and -l,
   respectively. If -l is present, the group size and similar
   details can be found in the global variable nauty_stats.

PRUNE feature.

   By defining the C preprocessor variable PRUNE at compile time, gentourng
   can be made to call
        int PRUNE(graph *g,int n,int maxn) 
   for each intermediate (and final) graph, and reject it if 
   the value returned is nonzero.  The arguments are:

     g      = the graph in nauty format (m=1)
     n      = the number of vertices in g
     maxn   = the number of vertices for output 
	      (the value you gave on the command line to gentourng)

   gentourng constructs the graph starting with vertex 0, then adding
   vertices 1,2,3,... in that order.  Each graph in the sequence is
   an induced subgraph of all later graphs in the sequence.

   A call is made for all orders from 1 to maxn.  In testing for
   a uniform property (such as a forbidden subgraph or forbidden
   induced subgraph) it might save time to notice that a call to
   PRUNE for n implies that the call for n-1 already passed. 

   For very fast tests, it might be worthwhile using PREPRUNE as
   well or instead. It has the same meaning but is applied earlier
   and more often.

SUMMARY

   If the C preprocessor variable SUMMARY is defined at compile time, the
   procedure SUMMARY(nauty_counter nout, double cpu) is called just before
   the program exits.  The purpose is to allow reporting of statistics
   collected by PRUNE or OUTPROC.  The values nout and cpu are the output
   count and cpu time reported on the >Z line.
   Output should be written to stderr.

INSTRUMENT feature.

   If the C preprocessor variable INSTRUMENT is defined at compile time,
   extra code is inserted to collect statistics during execution, and
   more information is written to stderr at termination.

CALLING FROM A PROGRAM

   It is possible to call gentourng from another program instead of using it
   as a stand-alone program.  The main requirement is to change the name
   of the main program to be other than "main".  This is done by defining
   the preprocessor variable GENG_MAIN.  You might also like to define
   OUTPROC to be the name of a procedure to receive the graphs. To call
   the program you need to define an argument list argv[] consistent with
   the usual one; don't forget that argv[0] is the command name and not
   the first argument.  The value of argc is the number of strings in
   argv[]; that is, one more than the number of arguments.  See the
   sample program callgeng.c.


Counts:

                all                 strong          regular
   n        tournaments           tournaments     tournaments
  
   1                     1                     1             1
   2                     1                     0             1
   3                     2                     1             1
   4                     4                     1             1
   5                    12                     6             1
   6                    56                    35             5
   7                   456                   353             3
   8                  6880                  6008            85
   9                191536                178133            15
  10               9733056               9355949         13333
  11             903753248             884464590          1223
  12          154108311168          152310149735      19434757
  13        48542114686912        48234782263293       1495297
  14     28401423719122304     28304491788158056
  15  31021002160355166848  30964247546702883729   18400989629 

**************************************************************************

    Author:   B. D. McKay, Nov 2008.
              Copyright  B. McKay (2008).  All rights reserved.
              This software is subject to the conditions and waivers
              detailed in the file nauty.h.

**************************************************************************/

#define NAUTY_PGM  3   /* 1 = geng, 2 = genbg, 3 = gentourng */

#ifndef MAXN
#define MAXN 32         /* not more than max(32,WORDSIZE) */
#endif

#if MAXN > 32
 #error "Can't have MAXN greater than 32"
#endif

#define ONE_WORD_SETS
#include "gtools.h"   /* which includes nauty.h and stdio.h */

typedef unsigned int xword;

static void (*outproc)(FILE*,graph*,int);

static FILE *outfile;           /* file for output graphs */
static int connec;              /* 1 for -c, 0 for not */
boolean graph6;                 /* presence of -g */
boolean digraph6;               /* presence of -z */
boolean sparse6;                /* presence of -s */
boolean nooutput;               /* presence of -u */
boolean canonise;               /* presence of -l */
boolean quiet;                  /* presence of -q */
boolean header;                 /* presence of -h */
statsblk nauty_stats;
static int mindeg,maxdeg,maxn,mod,res;
static boolean regular;
#define PRUNEMULT 20   /* bigger -> more even split at greater cost */
static int min_splitlevel,odometer,splitlevel,multiplicity;
static graph gcan[MAXN];

#if MAXN <= 16
static xword xbit[] = {0x0001,0x0002,0x0004,0x0008,
                       0x0010,0x0020,0x0040,0x0080,
                       0x0100,0x0200,0x0400,0x0800,
                       0x1000,0x2000,0x4000,0x8000};

#define XNEXTBIT(x) \
    ((x)&0xFF ? 7-leftbit[(x)&0xFF] : 15-leftbit[((x)>>8)&0xFF])
#define XPOPCOUNT(x) (bytecount[((x)>>8)&0xFF] + bytecount[(x)&0xFF])
#elif MAXN <= 24
static xword xbit[] = {0x000001,0x000002,0x000004,0x000008,
                       0x000010,0x000020,0x000040,0x000080,
                       0x000100,0x000200,0x000400,0x000800,
                       0x001000,0x002000,0x004000,0x008000,
		       0x010000,0x020000,0x040000,0x080000,
		       0x100000,0x200000,0x400000,0x800000};

#define XNEXTBIT(x) \
    ((x)&0xFF ? 7-leftbit[(x)&0xFF] : \
      (x)&0xFF00 ? 15-leftbit[((x)>>8)&0xFF] : 23-leftbit[((x)>>16)&0xFF])
#define XPOPCOUNT(x) (bytecount[((x)>>8)&0xFF] \
	+ bytecount[((x)>>16)&0xFF] + bytecount[(x)&0xFF])
#else
static xword xbit[] = {0x00000001,0x00000002,0x00000004,0x00000008,
                       0x00000010,0x00000020,0x00000040,0x00000080,
                       0x00000100,0x00000200,0x00000400,0x00000800,
                       0x00001000,0x00002000,0x00004000,0x00008000,
                       0x00010000,0x00020000,0x00040000,0x00080000,
                       0x00100000,0x00200000,0x00400000,0x00800000,
                       0x01000000,0x02000000,0x04000000,0x08000000,
                       0x10000000,0x20000000,0x40000000,0x80000000};

#define XNEXTBIT(x) \
    ((x)&0xFF ? 7-leftbit[(x)&0xFF] : \
      (x)&0xFF00 ? 15-leftbit[((x)>>8)&0xFF] : \
        (x)&0xFF0000 ? 23-leftbit[((x)>>16)&0xFF] : \
			  31-leftbit[((x)>>24)&0xFF])
#define XPOPCOUNT(x) (bytecount[((x)>>8)&0xFF] \
        + bytecount[((x)>>16)&0xFF] + \
	+ bytecount[((x)>>24)&0xFF] + bytecount[(x)&0xFF])
#endif

typedef struct
{
    xword lo,hi;          /* work purposes for orbit calculation */
    xword xstart[MAXN+1]; /* index into xset[] for each cardinality */
    xword *xset;          /* array of all x-sets in card order */
    xword *xcard;         /* cardinalities of all x-sets */
    xword *xinv;          /* map from x-set to index in xset */
    xword *xorb;          /* min orbit representative */
} leveldata;

static leveldata data[MAXN];      /* data[n] is data for n -> n+1 */
static nauty_counter nodes[MAXN];     /* nodes at each level */
static nauty_counter nout;

#ifdef INSTRUMENT
static unsigned long rigidnodes[MAXN],fertilenodes[MAXN];
static unsigned long a1calls,a1nauty,a1succs;
static unsigned long a2calls,a2nauty,a2uniq,a2succs;
#endif

#ifdef PLUGIN
#include PLUGIN
#endif

#ifdef OUTPROC
extern void OUTPROC(FILE*,graph*,int);
#endif
#ifdef PRUNE
extern int PRUNE(graph*,int,int);
#endif
#ifdef PREPRUNE
extern int PREPRUNE(graph*,int,int);
#endif
#ifdef SUMMARY
extern void SUMMARY(nauty_counter,double);
#endif

/************************************************************************/

void
write_ascii(FILE *f, graph *g, int n)
/* write tournament g (n vertices) to file f in ascii format */
{
	char s[MAXN*(MAXN-1)/2+2];
	int i,j;
	size_t k;

	k = 0;
	for (i = 0; i < n-1; ++i)
	    for (j = i+1; j < n; ++j)
		if ((g[i] & bit[j])) s[k++] = '1'; else s[k++] = '0';

	s[k++] = '\n';
	s[k] = '\0';

        if (fwrite(s,1,k,f) != k || ferror(f))
            gt_abort(">E write_ascii : error on writing\n");
}

/************************************************************************/

void
writeg6x(FILE *f, graph *g, int n)
/* write graph g (n vertices) to file f in graph6 format */
{
	writeg6(f,g,1,n);
}

/************************************************************************/

void
writes6x(FILE *f, graph *g, int n)
/* write graph g (n vertices) to file f in sparse6 format */
{
	writes6(f,g,1,n);
}

/************************************************************************/

void
writed6x(FILE *f, graph *g, int n)
/* write graph g (n vertices) to file f in digraph6 format */
{
	writed6(f,g,1,n);
}

/***********************************************************************/

static void
nullwrite(FILE *f, graph *g, int n)
/* don't write graph g (n vertices) to file f */
{
}

/***********************************************************************/

static boolean
isstrong(graph *g, int n)
/* test if tournament g is strongly-connected 
 * This code is strictly for tournaments only.
 */
{
        setword seen,expanded,toexpand,allbits;
        int i;

	allbits = ALLMASK(n);

        seen = bit[0] | g[0];
        expanded = bit[0];

        while (seen != allbits && (toexpand = (seen & ~expanded))) /* not == */
        {
            i = FIRSTBITNZ(toexpand);
            expanded |= bit[i];
            seen |= g[i];
        }

	if (seen != allbits) return FALSE;
 
        seen = (allbits ^ g[0]);
        expanded = bit[0];

        while (seen != allbits && (toexpand = (seen & ~expanded))) /* not == */
        {
            i = FIRSTBITNZ(toexpand);
            expanded |= bit[i];
            seen |= (g[i] ^ allbits);
        }

	return seen == allbits;
}

/**********************************************************************/

static void
gcomplement(graph *g, graph *gc, int n)
/* Take the complement of g and put it in gc */
{
	int i;
	setword all;

	all = ~(setword)BITMASK(n-1);
	for (i = 0; i < n; ++i)
	    gc[i] = g[i] ^ all ^ bit[i];
}

/**************************************************************************/

static void
makeleveldata(void)
/* make the level data for each level */
{
        long h;
        int n,dmax,dmin;
        long ncj;
        leveldata *d;
	xword *xcard,*xinv;
	xword *xset,xw,tttn,nxsets;
	xword cw;
	xword i,j;

        for (n = 1; n < maxn; ++n)
        {
            dmax = n/2; 
	    if (maxdeg < dmax) dmax = maxdeg;
            dmin = mindeg - maxn + n + 1;
	    if (dmin < 0) dmin = 0;
            ncj = 1;
	    nxsets = (dmin == 0 ? 1 : 0);
            for (j = 1; j <= dmax; ++j)
            {
                ncj = (ncj * (n-j+1)) / j;
                if (j >= dmin) nxsets += ncj;
            }
            tttn = 1L << n;

            d = &data[n];

            d->xset = xset = (xword*) calloc(nxsets,sizeof(xword));
            d->xcard = xcard = (xword*) calloc(nxsets,sizeof(xword));
            d->xinv = xinv = (xword*) calloc(tttn,sizeof(xword));
            d->xorb = (xword*) calloc(nxsets,sizeof(xword));

            if (xset==NULL || xcard==NULL || xinv==NULL || d->xorb==NULL)
            {
                fprintf(stderr,
			">E gentourng: calloc failed in makeleveldata()\n");
                exit(2);
            }

            j = 0;

            for (i = 0;; ++i)
	    {
                if ((h = XPOPCOUNT(i)) <= dmax && h >= dmin)
                {
                    xset[j] = i;
                    xcard[j] = h;
                    ++j;
                }
		if (i == (xword)((1L<<n)-1)) break;
	    }

            if (j != nxsets)
            {
                fprintf(stderr,">E gentourng: j=%u mxsets=%u\n",
                        j,(unsigned)nxsets);
                exit(2);
            }

            h = 1;
            do
                h = 3 * h + 1;
            while (h < nxsets);

            do
            {
                for (i = h; i < nxsets; ++i)
                {
                    xw = xset[i];
                    cw = xcard[i];
                    for (j = i; xcard[j-h] > cw ||
                                (xcard[j-h] == cw && xset[j-h] > xw); )
                    {
                        xset[j] = xset[j-h];
                        xcard[j] = xcard[j-h];
                        if ((j -= h) < h) break;
                    }
                    xset[j] = xw;
                    xcard[j] = cw;
                }
                h /= 3;
            }
            while (h > 0);

            for (i = 0; i < nxsets; ++i) xinv[xset[i]] = i;

            d->xstart[0] = 0;
            for (i = 1; i < nxsets; ++i)
                if (xcard[i] > xcard[i-1]) d->xstart[xcard[i]] = i;
            d->xstart[xcard[nxsets-1]+1] = nxsets;
        }
}

/**************************************************************************/

static void
userautomproc(int count, int *p, int *orbits,
              int numorbits, int stabvertex, int n)
/* form orbits on powerset of VG
   called by nauty;  operates on data[n] */
{
	xword i,j1,j2,moved,pi,pxi;
        xword lo,hi;
        xword *xorb,*xinv,*xset,w;

        xorb = data[n].xorb;
        xset = data[n].xset;
        xinv = data[n].xinv;
        lo = data[n].lo;
        hi = data[n].hi;

        if (count == 1)                         /* first automorphism */
            for (i = lo; i < hi; ++i) xorb[i] = i;

        moved = 0;
        for (i = 0; i < n; ++i)
            if (p[i] != i) moved |= xbit[i];

        for (i = lo; i < hi; ++i)
        {
            if ((w = xset[i] & moved) == 0) continue;
            pxi = xset[i] & ~moved;
            while (w)
            {
                j1 = XNEXTBIT(w);
                w ^= xbit[j1];
                pxi |= xbit[p[j1]];
            }
            pi = xinv[pxi];

            j1 = xorb[i];
            while (xorb[j1] != j1) j1 = xorb[j1];
            j2 = xorb[pi];
            while (xorb[j2] != j2) j2 = xorb[j2];

            if      (j1 < j2) xorb[j2] = xorb[i] = xorb[pi] = j1;
            else if (j1 > j2) xorb[j1] = xorb[i] = xorb[pi] = j2;
        }
}

/*****************************************************************************
*                                                                            *
*  refinex(g,lab,ptn,level,numcells,count,active,goodret,code,m,n) is a      *
*  custom version of refine() which can exit quickly if required.            *
*                                                                            *
*  Only use at level==0.                                                     *
*  goodret : whether to do an early return for code 1                        *
*  code := -1 for n-1 not max, 0 for maybe, 1 for definite                   *
*                                                                            *
*****************************************************************************/

static void
refinex(graph *g, int *lab, int *ptn, int level, int *numcells,
     int *count, set *active, boolean goodret, int *code, int m, int n)
{
        int i,c1,c2,labc1;
        setword x,lact;
        int split1,split2,cell1,cell2;
        int cnt,bmin,bmax;
        set *gptr;
        setword workset;
        int workperm[MAXN];
        int bucket[MAXN+2];

        if (n == 1)
        {
            *code = 1;
            return;
        }

        *code = 0;
	lact = *active;

        split1 = -1;
        while (*numcells < n && lact)
        {
	    TAKEBIT(split1,lact);
            
            for (split2 = split1; ptn[split2] > 0; ++split2) {}
            if (split1 == split2)       /* trivial splitting cell */
            {
                gptr = GRAPHROW(g,lab[split1],1);
                for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
                {
                    for (cell2 = cell1; ptn[cell2] > 0; ++cell2) {}
                    if (cell1 == cell2) continue;

                    c1 = cell1;
                    c2 = cell2;
                    while (c1 <= c2)
                    {
                        labc1 = lab[c1];
                        if (ISELEMENT1(gptr,labc1))
                            ++c1;
                        else
                        {
                            lab[c1] = lab[c2];
                            lab[c2] = labc1;
                            --c2;
                        }
                    }
                    if (c2 >= cell1 && c1 <= cell2)
                    {
                        ptn[c2] = 0;
                        ++*numcells;
			lact |= bit[c1];
                    }
                }
            }

            else        /* nontrivial splitting cell */
            {
                workset = 0;
                for (i = split1; i <= split2; ++i) workset |= bit[lab[i]];

                for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
                {
                    for (cell2 = cell1; ptn[cell2] > 0; ++cell2) {}
                    if (cell1 == cell2) continue;
                    i = cell1;
                    if ((x = workset & g[lab[i]]) != 0) cnt = POPCOUNT(x);
                    else                                cnt = 0;
                    count[i] = bmin = bmax = cnt;
                    bucket[cnt] = 1;
                    while (++i <= cell2)
                    {
                        if ((x = workset & g[lab[i]]) != 0)
                            cnt = POPCOUNT(x);
                        else
                            cnt = 0;

                        while (bmin > cnt) bucket[--bmin] = 0;
                        while (bmax < cnt) bucket[++bmax] = 0;
                        ++bucket[cnt];
                        count[i] = cnt;
                    }
                    if (bmin == bmax) continue;
                    c1 = cell1;
                    for (i = bmin; i <= bmax; ++i)
                        if (bucket[i])
                        {
                            c2 = c1 + bucket[i];
                            bucket[i] = c1;
                            if (c1 != cell1)
                            {
                                lact |= bit[c1];
                                ++*numcells;
                            }
                            if (c2 <= cell2) ptn[c2-1] = 0;
                            c1 = c2;
                        }
                    for (i = cell1; i <= cell2; ++i)
                        workperm[bucket[count[i]]++] = lab[i];
                    for (i = cell1; i <= cell2; ++i) lab[i] = workperm[i];
                }
            }

            if (ptn[n-2] == 0)
            {
                if (lab[n-1] == n-1)
                {
                    *code = 1;
                    if (goodret) return;
                }
                else
                {
                    *code = -1;
                    return;
                }
            }
            else
            {
                i = n - 1;
                while (TRUE)
                {
                    if (lab[i] == n-1) break;
                    --i;
                    if (ptn[i] == 0)
                    {
                        *code = -1;
                        return;
                    }
                }
            }
        }
}

/**************************************************************************/

static void
makecanon(graph *g, graph *gcan, int n)
/* gcan := canonise(g) */
{
        int lab[MAXN],ptn[MAXN],orbits[MAXN];
        static DEFAULTOPTIONS_GRAPH(options);
        setword workspace[50];

        options.getcanon = TRUE;
        options.digraph = TRUE;

        nauty(g,lab,ptn,NULL,orbits,&options,&nauty_stats,
              workspace,50,1,n,gcan);
}

/**************************************************************************/

static boolean
accept1(graph *g, int n, xword x, graph *gx, int *deg, boolean *rigid)
/* decide if n in theta(g+x) -  version for n+1 < maxn */
{
        int i;
        int lab[MAXN],ptn[MAXN],orbits[MAXN];
        int count[MAXN];
        graph h[MAXN];
        int nx,numcells,code;
        int i0,i1,degn;
        set active[MAXM];
        statsblk stats;
        static DEFAULTOPTIONS_GRAPH(options);
        setword workspace[50];

#ifdef INSTRUMENT
        ++a1calls;
#endif

        nx = n + 1;
        for (i = 0; i < n; ++i) gx[i] = g[i];
        gx[n] = 0;
        deg[n] = degn = XPOPCOUNT(x);

        for (i = 0; i < n; ++i)
	{
	    if ((xbit[i] & x))
		gx[n] |= bit[i];
	    else
	    {
		gx[i] |= bit[n];
		++deg[i];
	    }
	}

#ifdef PREPRUNE
	if (PREPRUNE(gx,n+1,maxn)) return FALSE;
#endif

        i0 = 0;
        i1 = n;
        for (i = 0; i < nx; ++i)
        {
            if (deg[i] == degn) lab[i1--] = i;
            else                lab[i0++] = i;
            ptn[i] = 1;
        }
        ptn[n] = 0;
        if (i0 == 0)
        {
            numcells = 1;
            active[0] = bit[0];
        }
        else
        {
            numcells = 2;
            active[0] = bit[0] | bit[i1+1];
            ptn[i1] = 0;
        }
        refinex(gx,lab,ptn,0,&numcells,count,active,FALSE,&code,1,nx);

        if (code < 0) return FALSE;

        if (numcells == nx)
        {
            *rigid = TRUE;
#ifdef INSTRUMENT
            ++a1succs;
#endif
            return TRUE;
        }

        options.getcanon = TRUE;
        options.digraph = TRUE;
        options.defaultptn = FALSE;
	options.userautomproc = userautomproc;
     /*
        if (!regular || nx != maxn-1) options.userautomproc = userautomproc;
	else                          options.userautomproc = NULL;
     */

        active[0] = 0;
#ifdef INSTRUMENT
        ++a1nauty;
#endif
        nauty(gx,lab,ptn,active,orbits,&options,&stats,workspace,50,1,nx,h);

        if (orbits[lab[n]] == orbits[n])
        {
            *rigid = stats.numorbits == nx;
#ifdef INSTRUMENT
            ++a1succs;
#endif
            return TRUE;
        }
        else
            return FALSE;
}

/**************************************************************************/

static boolean
hitinvar(graph *g, int *invar, int n)
/* make hitting invariant
 *    return FALSE if n-1 not maximal else return TRUE */
{
        setword x,y,z;
        int inv,i,v,d;

        for (v = n-1; v >= 0; --v)
        {
            inv = 0;
            x = y = g[v];
	    while (y)
	    {
		TAKEBIT(i,y);
		z = x & g[i];
		d = POPCOUNT(z);
		if (d > inv) inv = d;
	    }

            invar[v] = inv;
            if (v < n-1 && inv > invar[n-1]) return FALSE;
        }
        return TRUE;
}

/**************************************************************************/

static boolean
accept2(graph *g, int n, xword x, graph *gx, int *deg, boolean nuniq)
/* decide if n in theta(g+x)  --  version for n+1 == maxn */
{
        int i;
        int lab[MAXN],ptn[MAXN],orbits[MAXN];
        int degx[MAXN],invar[MAXN];
        setword vmax,gv,gxn;
        int qn,qv;
        int count[MAXN];
        int nx,numcells,code;
        int degn,i0,i1,j,j0,j1;
        set active[MAXM];
        statsblk stats;
        static DEFAULTOPTIONS_GRAPH(options);
        setword workspace[50];
	boolean cheapacc;

#ifdef INSTRUMENT
        ++a2calls;
        if (nuniq) ++a2uniq;
#endif
        nx = n + 1;
        gxn = 0;

        for (i = 0; i < n; ++i)
	{
	    if ((xbit[i] & x))
	    {
		gxn |= bit[i];
		gx[i] = g[i];
		degx[i] = deg[i];
	    }
	    else
	    {
		gx[i] = g[i] | bit[n];
		degx[i] = deg[i] + 1;
	    }
	}
	gx[n] = gxn;
        degx[n] = degn = XPOPCOUNT(x);
	    
#ifdef PREPRUNE
        if (PREPRUNE(gx,n+1,maxn)) return FALSE;
#endif

        if (nuniq)
        {
#ifdef INSTRUMENT
            ++a2succs;
#endif
            if (canonise) makecanon(gx,gcan,nx);
            return TRUE;
        }

        i0 = 0;
        i1 = n;
        for (i = 0; i < nx; ++i)
        {
            if (degx[i] == degn) lab[i1--] = i;
            else                 lab[i0++] = i;
            ptn[i] = 1;
        }
        ptn[n] = 0;
        if (i0 == 0)
        {
            numcells = 1;
            active[0] = bit[0];
            if (!hitinvar(gx,invar,nx)) return FALSE;
            qn = invar[n];
            j0 = 0;
            j1 = n;
            while (j0 <= j1)
            {
                j = lab[j0];
                qv = invar[j];
                if (qv < qn)
                    ++j0;
                else
                {
                    lab[j0] = lab[j1];
                    lab[j1] = j;
                    --j1;
                }
            }
            if (j0 > 0)
            {
                if (j0 == n)
                {
#ifdef INSTRUMENT
                    ++a2succs;
#endif
                    if (canonise) makecanon(gx,gcan,nx);
                    return TRUE;
                }
                ptn[j1] = 0;
                ++numcells;
                active[0] |= bit[j0];
            }
        }
        else
        {
            numcells = 2;
            ptn[i1] = 0;
            active[0] = bit[0] | bit[i1+1];

            vmax = 0;
            for (i = i1+1; i < nx; ++i) vmax |= bit[lab[i]];

            gv = gx[n] & vmax;
            qn = POPCOUNT(gv);

            j0 = i1+1;
            j1 = n;
            while (j0 <= j1)
            {
                j = lab[j0];
                gv = gx[j] & vmax;
                qv = POPCOUNT(gv);
                if (qv > qn)
                    return FALSE;
                else if (qv < qn)
                    ++j0;
                else
                {
                    lab[j0] = lab[j1];
                    lab[j1] = j;
                    --j1;
                }
            }
            if (j0 > i1+1)
            {
                if (j0 == n)
                {
#ifdef INSTRUMENT
                    ++a2succs;
#endif
                    if (canonise) makecanon(gx,gcan,nx);
                    return TRUE;
                }
                ptn[j1] = 0;
                ++numcells;
                active[0] |= bit[j0];
            }
        }

        refinex(gx,lab,ptn,0,&numcells,count,active,TRUE,&code,1,nx);

        if (code < 0) return FALSE;

	cheapacc = FALSE;
	if (code > 0) cheapacc = TRUE;
    
        if (cheapacc)
        {
#ifdef INSTRUMENT
            ++a2succs;
#endif
            if (canonise) makecanon(gx,gcan,nx);
            return TRUE;
        }

        options.getcanon = TRUE;
        options.digraph = TRUE;
        options.defaultptn = FALSE;

        active[0] = 0;
#ifdef INSTRUMENT
        ++a2nauty;
#endif
        nauty(gx,lab,ptn,active,orbits,&options,&stats,workspace,50,1,nx,gcan);

        if (orbits[lab[n]] == orbits[n])
        {
#ifdef INSTRUMENT
            ++a2succs;
#endif
            if (canonise) makecanon(gx,gcan,nx);
            return TRUE;
        }
        else
            return FALSE;
}

/**************************************************************************/

static void
genextend(graph *g, int n, int *deg, boolean rigid)
/* extend from n to n+1 -- version for general graphs */
{
        xword x,dlow,dhigh,dcrit;
        xword *xset,*xcard,*xorb;
	xword i,imin,imax;
        int nx,xc,j,dmax;
        int xlb,xub,xlbx,xubx;
        graph gx[MAXN];
        int degx[MAXN];
        boolean rigidx;
	boolean subconnec;

#ifdef INSTRUMENT
        boolean haschild;

        haschild = FALSE;
        if (rigid) ++rigidnodes[n];
#endif
        nx = n + 1;
        ++nodes[n];

	if (regular && nx == maxn)
	{
	    x = 0;
	    for (i = 0; i < n; ++i)
		if (deg[i] == maxdeg) x |= xbit[i];

            if (accept2(g,n,x,gx,deg,FALSE))
            {
#ifdef PRUNE
                if (!PRUNE(gx,nx,maxn))
#endif
		{
#ifdef INSTRUMENT
                    haschild = TRUE;
#endif
		    ++nout;
                    (*outproc)(outfile,canonise ? gcan : gx,nx);
		}
            }
#ifdef INSTRUMENT
            if (haschild) ++fertilenodes[n];
#endif
	    return;
	}

        dmax = deg[n-1];

	xlb = mindeg + n + 1 - maxn;
	if (0 > xlb) xlb = 0;
	xub = dmax+1;
	if (n/2 < xub) xub = n/2;
	if (maxdeg < xub) xub = maxdeg;

	if (xlb > xub) return;

        dlow = dcrit = dhigh = 0;
        for (i = 0; i < n; ++i)
	{
            if (deg[i] == dmax) dlow |= xbit[i];
	    if (deg[i] == maxdeg) dhigh |= xbit[i];
	    if (deg[i] == mindeg + n - maxn) dcrit |= xbit[i];
	}

        if (XPOPCOUNT(dhigh) > xlb) xlb = XPOPCOUNT(dhigh);
	if (n-XPOPCOUNT(dcrit) < xub) xub = n - XPOPCOUNT(dcrit);
	if (xub == dmax+1 && XPOPCOUNT(dlow)+dmax >= n) --xub;
        if (xlb > xub) return;

#ifdef PRUNE 
        if (PRUNE(g,n,maxn)) return; 
#endif 

        imin = data[n].xstart[xlb];
        imax = data[n].xstart[xub+1];
        xset = data[n].xset;
        xcard = data[n].xcard;
        xorb = data[n].xorb;

        if (nx == maxn)
	{
	    subconnec = connec && isstrong(g,n);
            for (i = imin; i < imax; ++i)
            {
                if (!rigid && xorb[i] != i) continue;
                x = xset[i];
                xc = xcard[i];
                if (xc == dmax+1 && (x & dlow) != 0) continue;
 		if ((dhigh & ~x) != 0) continue;
		if ((dcrit & x) != 0) continue;

                if (accept2(g,n,x,gx,deg,
                            xc < dmax || (xc == dmax && (x & dlow) == 0)))
                    if (!connec || (subconnec && x != 0) || isstrong(gx,nx))
                    {
#ifdef PRUNE
                        if (!PRUNE(gx,nx,maxn))
#endif
			{
#ifdef INSTRUMENT
                            haschild = TRUE;
#endif
			    ++nout;
                            (*outproc)(outfile,canonise ? gcan : gx,nx);
			}
                    }
            }
	}
        else
	{
            for (i = imin; i < imax; ++i)
            {
                if (!rigid && xorb[i] != i) continue;
                x = xset[i];
                xc = xcard[i];
                if (xc == dmax+1 && (x & dlow) != 0) continue;
 		if ((dhigh & ~x) != 0) continue;
		if ((dcrit & x) != 0) continue;
                if (nx == splitlevel)
                {
                    if (odometer-- != 0) continue;
                    odometer = mod - 1;
                }

                for (j = 0; j < n; ++j) degx[j] = deg[j];
                xlbx = mindeg+nx+1-maxn;
		if (xlbx < 0) xlbx = 0;
		xubx = xc + 1;
		if (maxdeg < xubx) xubx = maxdeg;
		if (nx/2 < xubx) xubx = nx/2;
                if (xlbx > xubx) continue;

                data[nx].lo = data[nx].xstart[xlbx];
                data[nx].hi = data[nx].xstart[xubx+1];
                if (accept1(g,n,x,gx,degx,&rigidx))
                {
#ifdef INSTRUMENT
                    haschild = TRUE;
#endif
                    genextend(gx,nx,degx,rigidx);
                }
            }
	}

        if (n == splitlevel-1 && n >= min_splitlevel
                && nodes[n] >= multiplicity)
            --splitlevel;
#ifdef INSTRUMENT
        if (haschild) ++fertilenodes[n];
#endif
}

/**************************************************************************/
/**************************************************************************/

int
#ifdef GENG_MAIN
GENG_MAIN(int argc, char *argv[])
#else
main(int argc, char *argv[])
#endif
{
        char *arg;
        boolean badargs,gotd,gotD,gotf,gotmr;
	boolean secret,connec1;
	char *outfilename,sw;
        int i,j,argnum;
        graph g[1];
        int deg[1];
	int splitlevinc;
        double t1,t2;
	char msg[201];

	HELP; PUTVERSION;
	nauty_check(WORDSIZE,1,MAXN,NAUTYVERSIONID);

	if (MAXN > 32 || MAXN > WORDSIZE || MAXN > 8*sizeof(xword))
	{
	    fprintf(stderr,
		    "gentourng: incompatible MAXN, WORDSIZE, or xword\n");
	    fprintf(stderr,"--See notes in program source\n");
	    exit(1);
	}

        badargs = FALSE;
	graph6 = FALSE;
	digraph6 = FALSE;
	sparse6 = FALSE;
        nooutput = FALSE;
        canonise = FALSE;
	header = FALSE;
	outfilename = NULL;
	secret = FALSE;
	connec1 = FALSE;

        maxdeg = MAXN; 
	mindeg = 0;
	splitlevinc = 0;
	
	gotd = gotD = gotf = gotmr = FALSE;

        argnum = 0;
        for (j = 1; !badargs && j < argc; ++j)
        {
            arg = argv[j];
            if (arg[0] == '-' && arg[1] != '\0')
            {
		++arg;
		while (*arg != '\0')
		{
		    sw = *arg++;
		         SWBOOLEAN('u',nooutput)
		    else SWBOOLEAN('g',graph6)
		    else SWBOOLEAN('z',digraph6)
		    else SWBOOLEAN('s',sparse6)
		    else SWBOOLEAN('l',canonise)
		    else SWBOOLEAN('h',header)
		    else SWBOOLEAN('q',quiet)
		    else SWBOOLEAN('c',connec1)
		    else SWBOOLEAN('$',secret)
		    else SWINT('d',gotd,mindeg,"gentourng -d")
		    else SWINT('D',gotD,maxdeg,"gentourng -D")
#ifdef PLUGIN_SWITCHES
PLUGIN_SWITCHES
#endif
                    else badargs = TRUE;
		}
            }
            else if (arg[0] == '-' && arg[1] == '\0')
		gotf = TRUE;
	    else
            {
                if (argnum == 0)
                {
                    if (sscanf(arg,"%d",&maxn) != 1) badargs = TRUE;
		    ++argnum;
		}
		else if (gotf)
		    badargs = TRUE;
		else
		{
		    if (!gotmr)
		    {
			if (sscanf(arg,"%d/%d",&res,&mod) == 2)
                        { 
                            gotmr = TRUE; 
                            continue; 
                        }
		    }
		    if (!gotf)
		    {
			outfilename = arg;
			gotf = TRUE;
			continue;
		    }
                }
            }
        }

        if (argnum == 0)
            badargs = TRUE;
        else if (maxn < 1 || maxn > MAXN)
        {
            fprintf(stderr,">E gentourng: n must be in the range 1..%d\n",MAXN);
            badargs = TRUE;
        }

        if (!gotmr)
        {
            mod = 1;
            res = 0;
        }

	if (connec1) connec = 1;
	else         connec = 0;

        if (maxdeg >= maxn) maxdeg = maxn - 1;
	if (mindeg < 0) mindeg = 0;

	if (!badargs &&
	     (maxdeg < mindeg || 2*maxdeg < maxn-1 || 2*mindeg > maxn-1))
	{
            fprintf(stderr,">E gentourng: impossible degree bounds\n");
            badargs = TRUE;
        }

	if (connec && mindeg < 1 && maxn > 1) mindeg = 1;
        if (connec && maxdeg == maxn-1 && maxn > 1) maxdeg = maxn - 2;

        if (!badargs && (res < 0 || res >= mod))
        {
            fprintf(stderr,">E gentourng: must have 0 <= res < mod\n");
            badargs = TRUE;
        }

        if (badargs)
        {
            fprintf(stderr,">E Usage: %s\n",USAGE);
	    GETHELP;
            exit(1);
        }

	if ((graph6!=0) + (sparse6!=0) + (digraph6!=0) + (nooutput!=0) > 1)
	    gt_abort(">E gentourng: -ungzs are incompatible\n");

#ifdef OUTPROC
        outproc = OUTPROC;
#else
        if      (nooutput) outproc = nullwrite;
	else if (sparse6)  outproc = writes6x;
        else if (graph6)   outproc = writeg6x;
        else if (digraph6) outproc = writed6x;
	else               outproc = write_ascii;
#endif

#ifdef PLUGIN_INIT
PLUGIN_INIT
#endif

	for (i = 0; i < maxn; ++i)  nodes[i] = 0;

	if (nooutput)
	    outfile = stdout;
	else if (!gotf || outfilename == NULL)
	{
            outfilename = "stdout";
            outfile = stdout;
        }
        else if ((outfile = fopen(outfilename,"w")) == NULL)
        {
            fprintf(stderr,
                  ">E gentourng: can't open %s for writing\n",outfilename);
            gt_abort(NULL);
        }

        multiplicity = PRUNEMULT * mod;

	if (!quiet)
	{
	    msg[0] = '\0';
	    if (strlen(argv[0]) > 75)
		fprintf(stderr,">A %s",argv[0]);
	    else
		CATMSG1(">A %s",argv[0]);
	   
	    CATMSG2(" -%s%s",connec1  ? "c" : "",canonise ? "l" : "");
	 /*
	    if (mod > 1)
		CATMSG2("X%dx%d",splitlevinc,multiplicity);
	 */
	    CATMSG3("d%dD%d n=%d",mindeg,maxdeg,maxn);
	    if (mod > 1) CATMSG2(" class=%d/%d",res,mod);
	    CATMSG0("\n");
	    fputs(msg,stderr);
	    fflush(stderr);
	}

        g[0] = 0;
        deg[0] = 0;

        t1 = CPUTIME;

	if (header)
	{
	    if (sparse6)
	    {
                writeline(outfile,SPARSE6_HEADER);
		fflush(outfile);
	    }
	    else if (!nooutput)
	    {
		writeline(outfile,GRAPH6_HEADER);
		fflush(outfile);
	    }
	}

        nout = 0;

        if (maxn == 1)
        {
            if (res == 0)
            {
                ++nout;
                (*outproc)(outfile,g,1);
            }
        }
        else
        {
            makeleveldata();

	    if (maxn >= 14 && mod > 1)     splitlevel = maxn - 4;
	    else if (maxn >= 6 && mod > 1) splitlevel = maxn - 3;
	    else                           splitlevel = -1;

	    splitlevel += splitlevinc;
	    if (splitlevel > maxn - 1) splitlevel = maxn - 1;
	    if (splitlevel < 3) splitlevel = -1;

	    min_splitlevel = 6;
	    odometer = secret ? -1 : res;
	    regular = mindeg == maxdeg;
	    if (regular) connec = FALSE; /* All reg tourns are strong! */

            genextend(g,1,deg,TRUE);
	}
        t2 = CPUTIME;

#ifdef INSTRUMENT
        fprintf(stderr,"\n>N node counts\n");
        for (i = 1; i < maxn; ++i)
	{
            fprintf(stderr," level %2d: \n",i);
	    fprintf(stderr,COUNTER_FMT " (%lu rigid, %lu fertile)\n",
                           nodes[i],rigidnodes[i],fertilenodes[i]);
	}
        fprintf(stderr,">A1 %lu calls to accept1, %lu nauty, %lu succeeded\n",
                        a1calls,a1nauty,a1succs);
        fprintf(stderr,
             ">A2 %lu calls to accept2, %lu nuniq, %lu nauty, %lu succeeded\n",
                        a2calls,a2uniq,a2nauty,a2succs);
        fprintf(stderr,"\n");
#endif

#ifdef SUMMARY
	SUMMARY(nout,t2-t1);
#endif

	if (!quiet)
	{
            fprintf(stderr,">Z " COUNTER_FMT " graphs generated in %3.2f sec\n",
                    nout,t2-t1);
	}

#ifdef GENG_MAIN
	for (i = 1; i < maxn; ++i)
	{
	    free(data[i].xorb);
	    free(data[i].xset);
	    free(data[i].xinv);
	    free(data[i].xcard);
	}
	return 0;
#else
        exit(0);
#endif
}
