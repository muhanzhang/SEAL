/* genbg.c : version 2.4; B D McKay, 20 Jan 2016. */

/* TODO: consider colour swaps */

#define USAGE \
"genbg [-c -ugs -vq -lzF] [-Z#] [-D#] [-A] [-d#|-d#:#] [-D#|-D#:#] n1 n2 \n\
                [mine[:maxe]] [res/mod] [file]"

#define HELPTEXT \
" Find all bicoloured graphs of a specified class.\n\
\n\
  n1   : the number of vertices in the first class\n\
  n2   : the number of vertices in the second class\n\
 mine:maxe : a range for the number of edges\n\
              #:0 means '# or more' except in the case 0:0\n\
  res/mod : only generate subset res out of subsets 0..mod-1\n\
   file : the name of the output file (default stdout)\n\
  -c    : only write connected graphs\n\
  -z    : all the vertices in the second class must have\n\
          different neighbourhoods\n\
  -F    : the vertices in the second class must have at least two\n\
          neighbours of degree at least 2\n\
  -L    : there is no vertex in the first class whose removal leaves\n\
          the vertices in the second class unreachable from each other\n\
  -Z#   : two vertices in the second class may have at most # common nbrs\n\
  -A    : no vertex in the second class has a neighbourhood whith is a\n\
          subset of another vertex in the second class\n\
  -D#   : specify an upper bound for the maximum degree\n\
          Example: -D6.  You can also give separate maxima for the\n\
          two parts, for example: -D5:6\n\
  -d#   : specify a lower bound for the minimum degree.\n\
          Again, you can specify it separately for the two parts: -d1:2\n\
  -g    : use graph6 format for output (default)\n\
  -s    : use sparse6 format for output\n\
  -a    : use Greechie diagram format for output\n\
  -u    : do not output any graphs, just generate and count them\n\
  -v    : display counts by number of edges to stderr\n\
  -l    : canonically label output graphs (using the 2-part colouring)\n\
\n\
  -q    : suppress auxiliary output\n\
\n\
  See program text for much more information.\n"

/*
Output formats.

  If -n is absent, any output graphs are written in graph6 format.

  If -n is present, any output graphs are written in nauty format.

   For a graph of n vertices, the output consists of n+1 long ints
   (even if a setword is shorter than a long int). The first contains
   n, and the others contain the adjacency matrix.  Long int i of
   the adjacency matrix (0 <= i <= n-1) is the set of neighbours of
   vertex i, cast from setword to long int.

PRUNE feature.

   By defining the C preprocessor variables PRUNE1 and/or PRUNE2 at
   compile time, you can filter the output of the program efficiently.
   The value of the variable is a function name with parameter list
   (graph *g, int *deg, int n1, int n2, int maxn2)

   The function will be called for each intermediate graph generated
   by the program, including output graphs.  The parameters are:
      g   = the graph in nauty format
      deg = an array giving the degrees of the vertices
      n1  = the number of vertices in the first colour class
            (same as the n1 parameter on the command line)
      n2  = the number of vertices in the second colour class
	    (this will always be at least 1)
      maxn2 = the value of n2 on the command line
   If n2=maxn2, the graph has the output size.

   If the function returns a non-zero value, neither this graph nor
   any of its descendants will be written to the output.

   PRUNE1 and PRUNE2 are functionally equivalent, but placed at different
   points in the program.  Essentially, use PRUNE1 for fast tests that
   eliminate many cases and PRUNE2 for slow tests that eliminate few.
   If in doubt, try it both ways and choose the fastest.
   You can use both PRUNE1 and PRUNE2 if you wish, in which case note
   that the PRUNE1 test has already been passed before PRUNE2 is applied.

   Vertices 0..n1-1 are always present.  The program works by successively
   adding one more vertex to the second colour class.  Vertex n1+n2-1 is
   always the last one that has been added, and (except for n2=1) the
   subgraph formed by deleting n1+n2-1 has previously been passed by the
   pruning function.

   Note that neither PRUNE1 nor PRUNE2 are called with n2=0, even if that
   is the output level.

   If -c is specified, the connectivity test has NOT been performed yet
   at the time the pruning function is called.  However the simplicity
   test indicated by -z HAS been performed if -z is specified.

OUTPROC feature.

   By defining the C preprocessor variable OUTPROC at compile time
   (for Unix the syntax is -DOUTPROC=procname on the cc command),
   genbg can be made to call a procedure of your manufacture with
   each output graph instead of writing anything. Your procedure
   needs to have type void and the argument list (FILE *f, graph *g,
   int n1, int n2). f is a stream open for writing (in fact, in the
   current version it is always stdout), g is the graph in nauty
   format, and n1,n2 are the numbers of vertices on each side. Your
   procedure can be in a separate file so long as it is linked with
   genbg. The global variables nooutput, nautyformat and canonise
   (all type boolean) can be used to test for the presence of the
   flags -u, -n and -l, respectively.

   For backward compatibility, it is possible to instead define OUTPROC1 to
   be the name of a procedure with argument list  (FILE *f, graph *g, int n).

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

**************************************************************************

    Author:   B. D. McKay, Oct 1994.     bdm@cs.anu.edu.au
              Copyright  B. McKay (1994-2008).  All rights reserved.
              This software is subject to the conditions and waivers
              detailed in the file nauty.h.
    1 May 2003 : fixed PRUNE feature
   13 Sep 2003 : added Greechie output, all outprocs have n1,n2
    9 Oct 2003 : changed -l to respect partition
   11 Apr 2007 : make >A line more atomic
   29 Aug 2008 : include PLUGIN insertion
   29 Nov 2008 : slight improvement of connectivity testing
   27 Jul 2011 : fixed error in PRUNE1 found by Stephen Hartke
    8 Jan 2012 : add antichains -A suggested by Andrew Juell
   23 Jan 2013 : fix splitlevinc initialization
   16 Feb 2014 : add a missing call to PRUNE2
   20 Jan 2016 : changed bigint to nauty_counter

**************************************************************************/

#define NAUTY_PGM  2   /* 1 = geng, 2 = genbg, 3 = gentourng */
#undef MAXN
#define MAXN WORDSIZE

#ifndef MAXN1
#define MAXN1 24        /* not more than 30 */
#endif
#define ONE_WORD_SETS
#include "gtools.h"   /* which includes nauty.h and stdio.h */

static void (*outproc)(FILE*,graph*,int,int);
#ifdef OUTPROC
extern void OUTPROC(FILE*,graph*,int,int);
#endif

#ifdef SUMMARY
extern void SUMMARY(nauty_counter,double);
#endif

static FILE *outfile;           /* file for output graphs */
static boolean connec;          /* presence of -c */
static boolean verbose;         /* presence of -v */
static boolean simple;          /* presence of -z */
boolean nautyformat;            /* presence of -n */
boolean nooutput;               /* presence of -u */
boolean canonise;               /* presence of -l */
boolean graph6;                 /* presence of -g */
boolean sparse6;                /* presence of -s */
boolean greout;                 /* presence of -a */
boolean quiet;                  /* presence of -q */
boolean footfree;               /* presence of -F */
boolean cutfree;                /* presence of -L */
boolean antichain;              /* presence of -A */
int class1size;                 /* same as n1 */
int maxcommon;                  /* -1 or value of -Z */
static int maxdeg1,maxdeg2,n1,maxn2,mine,maxe,nprune,mod,res,curres;
static int mindeg1,mindeg2;
static graph gcan[MAXN];
static int xval[MAXN];   /* x-bit version of second class, xval[0..] */

#if MAXN1 <= 16
static int xbit[] = {0x0001,0x0002,0x0004,0x0008,
                   0x0010,0x0020,0x0040,0x0080,
                   0x0100,0x0200,0x0400,0x0800,
                   0x1000,0x2000,0x4000,0x8000};

#define XNEXTBIT(x) \
    ((x)&0xFF ? 7-leftbit[(x)&0xFF] : 15-leftbit[((x)>>8)&0xFF])
#define XPOPCOUNT(x) (bytecount[((x)>>8)&0xFF] + bytecount[(x)&0xFF])
#elif MAXN1 <= 24
static int xbit[] = {0x000001,0x000002,0x000004,0x000008,
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
static int xbit[] = {0x00000001,0x00000002,0x00000004,0x00000008,
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
    int ne,dmax;         /* values used for xlb,xub calculation */
    int xlb,xub;         /* saved bounds on extension degree */
    int lo,hi;           /* work purposes for orbit calculation */
    int *xorb;           /* min orbit representative */
} leveldata;

#define IFLE1BITS(ww)  if (!((ww)&((ww)-1)))

static leveldata data[MAXN];      /* data[n] is data for n -> n+1 */
static nauty_counter ecount[1+MAXN*MAXN/4];  /* counts by number of edges */
static int xstart[MAXN+1];  /* index into xset[] for each cardinality */
static int *xset;           /* array of all x-sets in card order */
static int *xcard;          /* cardinalities of all x-sets */
static int *xinv;           /* map from x-set to index in xset */

#ifdef INSTRUMENT
static long nodes[MAXN],rigidnodes[MAXN],fertilenodes[MAXN];
static long a1calls,a1nauty,a1succs;
static long a2calls,a2nauty,a2uniq,a2succs;
#endif

#ifdef SPLITTEST
static unsigned long splitcases = 0;
#endif

#ifdef PLUGIN
#include PLUGIN
#endif

/************************************************************************/

void
writeny(FILE *f, graph *g, int n1, int n2) 
/* write graph g (n1+n2 vertices) to file f in y format */
{
    static char ybit[] = {32,16,8,4,2,1};
    char s[(MAXN*(MAXN-1)/2 + 5)/6 + 4];
    int i,j,k;
    char y,*sp;
    int n;

    n = n1 + n2;
    sp = s;
    *(sp++) = 0x40 | n;
    y = 0x40;

    k = -1;
    for (j = 1; j < n; ++j)
    for (i = 0; i < j; ++i)
    {
        if (++k == 6)
        {
            *(sp++) = y;
            y = 0x40;
            k = 0;
        }
        if (g[i] & bit[j]) y |= ybit[k];
    }
    if (n >= 2) *(sp++) = y;
    *(sp++) = '\n';
    *sp = '\0';

    if (fputs(s,f) == EOF || ferror(f))
    {
        fprintf(stderr,">E writeny : error on writing file\n");
        exit(2);
    }
}

/************************************************************************/
 
void
writeg6x(FILE *f, graph *g, int n1, int n2)
/* write graph g (n1+n2 vertices) to file f in graph6 format */
{
    writeg6(f,g,1,n1+n2);
}

/************************************************************************/

void
writes6x(FILE *f, graph *g, int n1, int n2)
/* write graph g (n1+n2 vertices) to file f in graph6 format */
{
    writes6(f,g,1,n1+n2);
}

/************************************************************************/

void
writegre(FILE *f, graph *g, int n1, int n2)
/* write graph g (n1+n2 vertices) to file f in Greechie diagram format */
{
    static char atomname[] = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ\
abcdefghijklmnopqrstuvwxyz!\"#$%&'()*-/:;<=>?@[\\]^_`{|}~";
    char grestr[MAXN*MAXN+MAXN+5];
    int i,j,k;
    setword gi;

    k = 0;
    for (i = n1; i < n1+n2; ++i)
    {
        if (i > n1) grestr[k++] = ',';
        gi = g[i];
        while (gi)
        {
            TAKEBIT(j,gi);
            grestr[k++] = atomname[j];
        }
    }
    grestr[k++] = '.';
    grestr[k++] = '\n';
    grestr[k] = '\0';
    if (fputs(grestr,f) == EOF || ferror(f))
    {
        fprintf(stderr,">E genbg : error on writing file\n");
        gt_abort(NULL);
    }
}

/***********************************************************************/

static void
nullwrite(FILE *f, graph *g, int n1, int n2)
/* don't write graph g (n1+n2 vertices) to file f */
{
}

/***********************************************************************/

#ifdef OUTPROC1

static void
write12(FILE *f, graph *g, int n1, int n2)
/* pass to OUTPROC1 */
{
    OUTPROC1(f,g,1,n1+n2);
}

#endif

/***********************************************************************/

void
writenauty(FILE *f, graph *g, int n1, int n2)
/* write graph g (n1+n2 vertices) to file f in nauty format */
{
    int nn;

    nn = n1+n2;

    if (fwrite((char *)&nn,sizeof(int),(size_t)1,f) != 1 ||
          fwrite((char*)g,sizeof(setword),(size_t)nn,f) != nn)
    {
        fprintf(stderr,">E writenauty : error on writing file\n");
        exit(2);
    }
}

/*********************************************************************/

static void
fragments(int *x, int nx, int *frag, int *nfrag)
/* For each v in union(x[0..nx-1]), find the components of the
   hypergraph x - v and add them to frag if there are more than one. */
/* This implementation is shocking.  Improve it! */
{
    int allx,i,j,v;
    int vbit,nw,w[MAXN];
    boolean done;

    allx = 0;
    for (i = 0; i < nx; ++i) allx |= x[i];

    *nfrag = 0;
    
    while (allx)
    {
        v = XNEXTBIT(allx);
        vbit = xbit[v];
        allx &= ~vbit;

        for (i = 0; i < nx; ++i) w[i] = x[i] & ~vbit;
        nw = nx;

        done = FALSE;
        while (!done && nw > 1)
        {
            done = TRUE;
            for (i = 0; i < nw-1; ++i)
            for (j = i+1; j < nw; )
                if ((w[i] & w[j]) != 0)
                {
                    w[i] |= w[j];
                    w[j] = w[nw-1];
                    --nw;
                    done = FALSE;
                }
                else
                    ++j;
        }
    
        if (nw > 1)
            for (i = 0; i < nw; ++i)
                frag[(*nfrag)++] = w[i];
    }
}
 
/*********************************************************************/

static boolean
isconnected(graph *g, int n)
/* test if g is connected */
{
    setword seen,expanded,toexpand,allbits;
    int i;

    allbits = ALLMASK(n);

    expanded = bit[n-1];
    seen = expanded | g[n-1];

    while (seen != allbits && (toexpand = (seen & ~expanded))) /* not == */
    {   
        i = FIRSTBITNZ(toexpand);
        expanded |= bit[i];
        seen |= g[i];
    }   

    return  seen == allbits;
}

/**************************************************************************/

static boolean
distinvar(graph *g, int *invar, int n1, int n2)
 /* make distance invariant/
    exit immediately FALSE if n-1 not maximal else exit TRUE
    Note: only invar[n1..n1+n2-1] set */
{
    int w,n;
    setword workset,frontier;
    setword sofar;
    int inv,d,v;

    n = n1 + n2;
    for (v = n-1; v >= n1; --v)
    {
        inv = 0;
        sofar = frontier = bit[v];
        for (d = 1; frontier != 0; ++d)
        {
            workset = 0;
            inv += POPCOUNT(frontier) ^ (0x57 + d);
            while (frontier)
            {
                w = FIRSTBITNZ(frontier);
                frontier &= ~bit[w];
                workset |= g[w];
            }
            frontier = workset & ~sofar;
            sofar |= frontier;
        }
        invar[v] = inv;
        if (v < n-1 && inv > invar[n-1]) return FALSE;
    }
    return TRUE;
}

/**************************************************************************/

static void
makeleveldata(void)
/* make the level data for each level */
{
    int i,j,h;
    int nn,nxsets,tttn;
    long ncj;
    leveldata *d;
    int xw,cw;

    nn = maxdeg2 <= n1 ? maxdeg2 : n1;
    ncj = nxsets = 1;
    for (j = 1; j <= nn; ++j)
    {
        ncj = (ncj * (n1 - j + 1)) / j;
        nxsets += ncj;
    }

    tttn = 1 << n1;
    xset = (int*) ALLOCS(nxsets,sizeof(int));
    xcard = (int*) ALLOCS(nxsets,sizeof(int));
    xinv = (int*) ALLOCS(tttn,sizeof(int));
    if (xset==NULL || xcard==NULL || xinv==NULL)
    {
        fprintf(stderr,">E genbg: malloc failed in makeleveldata()\n");
        exit(2);
    }

    j = 0;
    for (i = 0; i < tttn; ++i)
        if ((h = XPOPCOUNT(i)) <= maxdeg2)
        {
            xset[j] = i;
            xcard[j] = h;
            ++j;
        }

    if (j != nxsets)
    {
        fprintf(stderr,">E genbg: j=%d mxsets=%d\n",j,nxsets);
        exit(2);
    }

    /* The following is not SORT_OF_SORT 1, 2 or 3 */

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

    xstart[0] = 0;
    for (i = 1; i < nxsets; ++i)
        if (xcard[i] > xcard[i-1]) xstart[xcard[i]] = i;
    xstart[xcard[nxsets-1]+1] = nxsets;

    for (i = 0; i < maxn2; ++i)
    {

        d = &data[i];

        d->xorb = (int*) ALLOCS(nxsets,sizeof(int));

        if (d->xorb==NULL)
        {
            fprintf(stderr,">E genbg: malloc failed in makeleveldata()\n");
            exit(2);
        }

        d->ne = d->dmax = d->xlb = d->xub = -1;
    }
}

/**************************************************************************/

static UPROC
userautomproc(int count, int *p, int *orbits, int numorbits,
          int stabvertex, int n)
/* Automorphism procedure called by nauty
   Form orbits on powerset of VG
   Operates on data[n-n1] */
{
    int i,j1,j2;
    int moved,pxi,pi;
    int w,lo,hi;
    int *xorb;

    xorb = data[n-n1].xorb;
    lo = data[n-n1].lo;
    hi = data[n-n1].hi;

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
            w &= ~xbit[j1];
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
    int *count, set *active, boolean goodret,
    int *code, int m, int n)
{
    int i,c1,c2,labc1;
    setword x;
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
    split1 = -1;
    while (*numcells < n && ((split1 = nextelement(active,1,split1)) >= 0
                         || (split1 = nextelement(active,1,-1)) >= 0))
    {
        DELELEMENT1(active,split1);
        for (split2 = split1; ptn[split2] > 0; ++split2)
        {}
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
                    ADDELEMENT1(active,c1);
                }
            }
        }

        else        /* nontrivial splitting cell */
        {
            workset = 0;
            for (i = split1; i <= split2; ++i)
                workset |= bit[lab[i]];

            for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
            {
                for (cell2 = cell1; ptn[cell2] > 0; ++cell2) {}
                if (cell1 == cell2) continue;
                i = cell1;
                if ((x = workset & g[lab[i]]))     /* not == */
                    cnt = POPCOUNT(x);
                else
                    cnt = 0;
                count[i] = bmin = bmax = cnt;
                bucket[cnt] = 1;
                while (++i <= cell2)
                {
                    if ((x = workset & g[lab[i]])) /* not == */
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
                            ADDELEMENT1(active,c1);
                            ++*numcells;
                        }
                        if (c2 <= cell2) ptn[c2-1] = 0;
                        c1 = c2;
                    }
                for (i = cell1; i <= cell2; ++i)
                    workperm[bucket[count[i]]++] = lab[i];
                for (i = cell1; i <= cell2; ++i)
                    lab[i] = workperm[i];
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
            while (1)
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
makecanon(graph *g, graph *gcan, int n1, int n2)
/* gcan := canonise(g) */
{
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    setword active[1];
    int i;
    statsblk stats;
    static DEFAULTOPTIONS_GRAPH(options);
    setword workspace[50];

    options.writemarkers = FALSE;
    options.writeautoms = FALSE;
    options.getcanon = TRUE;
    options.defaultptn = FALSE;

    for (i = 0; i < n1+n2; ++i)
    {
        lab[i] = i;
        ptn[i] = 1;
    }
    ptn[n1-1] = ptn[n1+n2-1] = 0;
    EMPTYSET(active,1);
    ADDELEMENT(active,0);
    ADDELEMENT(active,n1);

    nauty(g,lab,ptn,active,orbits,&options,&stats,
                                         workspace,50,1,n1+n2,gcan);
}

/**************************************************************************/

static boolean
accept1(graph *g, int n2, int x, graph *gx, int *deg, boolean *rigid)
 /* decide if n2 in theta(g+x) -- version for n2+1 < maxn2 */
{
    int i,n;
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    int count[MAXN];
    graph h[MAXN];
    int xw;
    int nx,numcells,code;
    int i0,i1,degn;
    set active[MAXM];
    statsblk stats;
    static DEFAULTOPTIONS_GRAPH(options);
    setword workspace[50];

#ifdef INSTRUMENT
    ++a1calls;
#endif

    n = n1 + n2;
    nx = n + 1;
    for (i = 0; i < n; ++i) gx[i] = g[i];
    gx[n] = 0;
    deg[n] = degn = XPOPCOUNT(x);

    xw = x;
    while (xw)
    {
        i = XNEXTBIT(xw);
        xw &= ~xbit[i];
        gx[i] |= bit[n];
        gx[n] |= bit[i];
        ++deg[i];
    }
#ifdef PRUNE1
    if (PRUNE1(gx,deg,n1,n2+1,maxn2)) return FALSE;
#endif

    for (i = 0; i < n1; ++i)
    {
        lab[i] = i;
        ptn[i] = 1;
    }
    ptn[n1-1] = 0;

    i0 = n1;
    i1 = n;
    for (i = n1; i < nx; ++i)
    {
        if (deg[i] == degn) lab[i1--] = i;
        else                lab[i0++] = i;
        ptn[i] = 1;
    }

    ptn[n] = 0;

    if (i0 == n1)
    {
        numcells = 2;
        active[0] = bit[0] | bit[n1];
    }
    else
    {
        numcells = 3;
        active[0] = bit[0] | bit[n1] | bit[i1+1];
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
#ifdef PRUNE2
        if (PRUNE2(gx,deg,n1,n2+1,maxn2)) return FALSE;
#endif
        return TRUE;
    }

    options.writemarkers = FALSE;
    options.writeautoms = FALSE;
    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    options.userautomproc = userautomproc;

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
#ifdef PRUNE2
        if (PRUNE2(gx,deg,n1,n2+1,maxn2)) return FALSE;
#endif
        return TRUE;
    }
    else
        return FALSE;
}

/**************************************************************************/

static boolean
accept2(graph *g, int n2, int x, graph *gx, int *deg, boolean nuniq)
/* decide if n in theta(g+x) -- version for n+1 == maxn */
{
    int i,n;
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    int degx[MAXN],invar[MAXN];
    setword vmax,gv;
    int qn,qv;
    int count[MAXN];
    int xw;
    int nx,numcells,code;
    int degn,i0,i1,j,j0,j1;
    set active[MAXM];
    statsblk stats;
    static DEFAULTOPTIONS_GRAPH(options);
    setword workspace[50];

#ifdef INSTRUMENT
    ++a2calls;
    if (nuniq) ++a2uniq;
#endif
    n = n1 + n2;
    nx = n + 1;
    for (i = 0; i < n; ++i)
    {
        gx[i] = g[i];
        degx[i] = deg[i];
    }
    gx[n] = 0;
    degx[n] = degn = XPOPCOUNT(x);

    xw = x;
    while (xw)
    {
        i = XNEXTBIT(xw);
        xw &= ~xbit[i];
        gx[i] |= bit[n];
        gx[n] |= bit[i];
        ++degx[i];
    }
#ifdef PRUNE1
    if (PRUNE1(gx,degx,n1,n2+1,maxn2)) return FALSE;
#endif

    if (nuniq)
    {
#ifdef INSTRUMENT
        ++a2succs;
#endif
#ifdef PRUNE2
        if (PRUNE2(gx,degx,n1,n2+1,maxn2)) return FALSE;
#endif
        if (canonise) makecanon(gx,gcan,n1,n2+1);
        return TRUE;
    }

    for (i = 0; i < n1; ++i)
    {
        lab[i] = i;
        ptn[i] = 1;
    }
    ptn[n1-1] = 0;

    i0 = n1;
    i1 = n;
    for (i = n1; i < nx; ++i)
    {
        if (degx[i] == degn) lab[i1--] = i;
        else                 lab[i0++] = i;
        ptn[i] = 1;
    }

    ptn[n] = 0;

    if (i0 == n1)
    {
        numcells = 2;
        active[0] = bit[0] | bit[n1];

        if (!distinvar(gx,invar,n1,n2+1)) return FALSE;
        qn = invar[n];
        j0 = n1;
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
        if (j0 > n1)
        {
            if (j0 == n)
            {
#ifdef INSTRUMENT
                ++a2succs;
#endif
#ifdef PRUNE2
                if (PRUNE2(gx,degx,n1,n2+1,maxn2)) return FALSE;
#endif
                if (canonise) makecanon(gx,gcan,n1,n2+1);
                return TRUE;
            }
            ptn[j1] = 0;
            ++numcells;
            active[0] |= bit[j0];
        }
    }
    else
    {
        numcells = 3;
        ptn[i1] = 0;
        active[0] = bit[0] | bit[n1] | bit[i1+1];

        vmax = 0;
        j = MAXN;
        for (i = 0; i < n1; ++i)
            if (degx[i] < j && degx[i] > 0)
            {
                j = degx[i];
                vmax = bit[i];
            }
            else if (degx[i] == j)
                vmax |= bit[i];

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
#ifdef PRUNE2
                if (PRUNE2(gx,degx,n1,n2+1,maxn2)) return FALSE;
#endif
                if (canonise) makecanon(gx,gcan,n1,n2+1);
                return TRUE;
            }
            ptn[j1] = 0;
            ++numcells;
            active[0] |= bit[j0];
        }
    }

    refinex(gx,lab,ptn,0,&numcells,count,active,TRUE,&code,1,nx);

    if (code < 0)
        return FALSE;
    else if (code > 0 || numcells >= nx-4)
    {
#ifdef INSTRUMENT
        ++a2succs;
#endif
#ifdef PRUNE2
        if (PRUNE2(gx,degx,n1,n2+1,maxn2)) return FALSE;
#endif
        if (canonise) makecanon(gx,gcan,n1,n2+1);
        return TRUE;
    }

    options.writemarkers = FALSE;
    options.writeautoms = FALSE;
    options.getcanon = TRUE;
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
#ifdef PRUNE2
    if (PRUNE2(gx,degx,n1,n2+1,maxn2)) return FALSE;
#endif
        if (canonise) makecanon(gx,gcan,n1,n2+1);
        return TRUE;
    }
    else
        return FALSE;
}

/**************************************************************************/

static void
xbnds(int n2, int ne, int dmax)
/* find bounds on degree for vertex n2
   Store answer in data[*].*  */
{
    int xlb,xub,m;

    xlb = n2 == 0 ? (connec ? 1 : 0) : dmax;
    if (xlb < mindeg2) xlb = mindeg2;
    m = mine - ne - (maxn2 - n2 -1)*maxdeg2;
    if (m > xlb) xlb = m;

    xub = maxdeg2;
    m = (maxe - ne) / (maxn2 - n2);
    if (m < xub) xub = m;

    data[n2].ne = ne;
    data[n2].dmax = dmax;
    data[n2].xlb = xlb;
    data[n2].xub = xub;
}

/**************************************************************************/

static void
genextend(graph *g, int n2, int *deg, int ne, boolean rigid, int xlb, int xub)
/* extend from n2 to n2+1 */
{
    int x,y,d;
    int *xorb,xc;
    int nx,i,j,imin,imax,dmax;
    int xlbx,xubx,n;
    graph gx[MAXN];
    int degx[MAXN];
    boolean rigidx;
    int dneed,need,nfeet,hideg,deg1,ft[MAXN],nfrag,frag[MAXN];

#ifdef INSTRUMENT
    boolean haschild;

    haschild = FALSE;
    ++nodes[n2];
    if (rigid) ++rigidnodes[n2];
#endif

    n = n1 + n2;
    nx = n2 + 1;
    dmax = deg[n-1];

    d = 0;
    dneed = mindeg1 - maxn2 + n2;
    need = 0;
    hideg = 0;
    deg1 = 0;
    for (i = 0; i < n1; ++i)
    {
        if (deg[i] == maxdeg1) d |= xbit[i];
        if (deg[i] <= dneed) need |= xbit[i];
        if (deg[i] >= 2) hideg |= xbit[i];
        if (deg[i] == 1) deg1 |= xbit[i];
    }

    if (xlb < XPOPCOUNT(need)) xlb = XPOPCOUNT(need);
    if (xlb > xub) return;

    imin = xstart[xlb];
    imax = xstart[xub+1];
    xorb = data[n2].xorb;

    if (nx == maxn2)
    {
        if (footfree)
        {
            nfeet = 0;
            for (j = 0; j < n2; ++j)
            {
                x = xval[j] & hideg;
                IFLE1BITS(x) ft[nfeet++] = xval[j] & deg1;
            }
        }
        if (cutfree) fragments(xval,n2,frag,&nfrag);

        for (i = imin; i < imax; ++i)
        {
            if (!rigid && xorb[i] != i) continue;
            x = xset[i];
            xc = xcard[i];
            if ((x & d) != 0) continue;
            if ((need & ~x) != 0) continue;

            if (simple)
            {
                for (j = n2; --j >= 0;)
                    if (x == xval[j]) break;
                if (j >= 0) continue;
            }
            if (maxcommon >= 0)
            {
                for (j = n2; --j >= 0;)
                {
                    y = x & xval[j];
                    if (XPOPCOUNT(y) > maxcommon) break;
                }
                if (j >= 0) continue;
            }
	    if (antichain)
	    {
		for (j = 0; j < n2; ++j)
		    if ((xval[j] & ~x) == 0) break;
		if (j < n2) continue;
	    }
            if (footfree)
            {
                y = x & (hideg | deg1);
                IFLE1BITS(y) continue;
                for (j = 0; j < nfeet; ++j)
                    if ((x & ft[j]) == 0) break;
                if (j < nfeet) continue;
            }
            if (cutfree)
            {
                y = x & (hideg | deg1);
                IFLE1BITS(y) continue;
                for (j = 0; j < nfrag; ++j)
                    if ((x & frag[j]) == 0) break;
                if (j < nfrag) continue;
            }

            xval[n2] = x;

            if (nx == nprune)
            {
                if (curres == 0) curres = mod;
#ifdef SPLITTEST
                --curres;
                ++splitcases;
                continue;
#else
                if (--curres != 0) continue;
#endif
            }
            if (accept2(g,n2,x,gx,deg,xc > dmax))
                if (!connec || isconnected(gx,n+1))
                {
                    ++ecount[ne+xc];
#ifdef INSTRUMENT
                    haschild = TRUE;
#endif
                    (*outproc)(outfile,canonise ? gcan : gx,n1,n2+1);
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
            if ((x & d) != 0) continue;
            if ((need & ~x) != 0) continue;

            if (nx == nprune)
            {
                if (curres == 0) curres = mod;
#ifdef SPLITTEST
                --curres;
                ++splitcases;
                continue;
#else
                if (--curres != 0) continue;
#endif
            }

            if (simple)
            {
                for (j = n2; --j >= 0;)
                    if (x == xval[j]) break;
                if (j >= 0) continue;
            }
            if (maxcommon >= 0)
            {
                for (j = n2; --j >= 0;)
                {
                    y = x & xval[j];
                    if (XPOPCOUNT(y) > maxcommon) break;
                }
                if (j >= 0) continue;
            }
	    if (antichain)
	    {
		for (j = 0; j < n2; ++j)
		    if ((xval[j] & ~x) == 0) break;
		if (j < n2) continue;
	    }
            xval[n2] = x;

            for (j = 0; j < n; ++j) degx[j] = deg[j];
            if (data[nx].ne != ne+xc || data[nx].dmax != xc)
                xbnds(nx,ne+xc,xc);
            xlbx = data[nx].xlb;
            xubx = data[nx].xub;
            if (xlbx > xubx) continue;

            data[nx].lo = xstart[xlbx];
            data[nx].hi = xstart[xubx+1];
            if (accept1(g,n2,x,gx,degx,&rigidx))
            {
#ifdef INSTRUMENT
                haschild = TRUE;
#endif
                genextend(gx,nx,degx,ne+xc,rigidx,xlbx,xubx);
            }
        }
    }
#ifdef INSTRUMENT
    if (haschild) ++fertilenodes[n2];
#endif
}

/**************************************************************************/
/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *arg;
    boolean badargs,gotD,gote,gotf,gotmr,gotZ,gotd,gotX;
    long Dval1,Dval2;
    long dval1,dval2;
    int i,j,imin,imax,argnum,sw;
    int splitlevinc;
    graph g[MAXN1];
    int deg[MAXN1];
    nauty_counter nout;
    double t1,t2;
    char *outfilename;
    char msg[201];

    HELP; PUTVERSION;
    nauty_check(WORDSIZE,1,MAXN,NAUTYVERSIONID);

    if (MAXN > WORDSIZE || MAXN1 > 8*sizeof(int)-2)
    {
        fprintf(stderr,"genbg: incompatible MAXN, MAXN1 or WORDSIZE\n");
        fprintf(stderr,"--See notes in program source\n");
        exit(1);
    }

    badargs = FALSE;
    connec = FALSE;
    verbose = FALSE;
    nautyformat = FALSE;
    nooutput = FALSE;
    canonise = FALSE;
    greout = FALSE;
    simple = FALSE;
    graph6 = FALSE;
    sparse6 = FALSE;
    quiet = FALSE;
    footfree = FALSE;
    cutfree = FALSE;

    gote = FALSE;
    gotf = FALSE;
    gotmr = FALSE;
    gotD = FALSE;
    gotd = FALSE;
    gotZ = FALSE;
    gotX = FALSE;
    outfilename = NULL;

    maxdeg1 = maxdeg2 = MAXN;
    mindeg1 = mindeg2 = 0;

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
                     SWBOOLEAN('n',nautyformat)
                else SWBOOLEAN('u',nooutput)
                else SWBOOLEAN('q',quiet)
                else SWBOOLEAN('v',verbose)
                else SWBOOLEAN('z',simple)
                else SWBOOLEAN('F',footfree)
                else SWBOOLEAN('L',cutfree)
                else SWBOOLEAN('A',antichain)
                else SWBOOLEAN('l',canonise)
                else SWBOOLEAN('c',connec)
                else SWBOOLEAN('a',greout)
                else SWBOOLEAN('g',graph6)
                else SWBOOLEAN('s',sparse6)
                else SWINT('Z',gotZ,maxcommon,"genbg -Z")
                else SWINT('X',gotX,splitlevinc,"geng -X")
                else SWRANGE('D',":-",gotD,Dval1,Dval2,"genbg -D")
                else SWRANGE('d',":-",gotd,dval1,dval2,"genbg -d")
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
                if (sscanf(arg,"%d",&n1) != 1) badargs = TRUE;
                ++argnum;
            }
            else if (argnum == 1)
            {
                if (sscanf(arg,"%d",&maxn2) != 1) badargs = TRUE;
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
                if (!gote)
                {
                    if (sscanf(arg,"%d:%d",&mine,&maxe) == 2
                     || sscanf(arg,"%d-%d",&mine,&maxe) == 2)
                    {
                        gote = TRUE;
                        if (maxe == 0 && mine > 0) maxe = MAXN*MAXN/4;
                        continue;
                    }
                    else if (sscanf(arg,"%d",&mine) == 1)
                    {
                        gote = TRUE;
                        maxe = mine;
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

    if (argnum < 2)
        badargs = TRUE;
    else if (n1 < 1 || maxn2 < 0 || n1 > MAXN1 || n1+maxn2 > MAXN)
    {
        fprintf(stderr,
           ">E genbg: must have n1=1..%d, n1+n2=1..%d\n",MAXN1,MAXN);
        badargs = TRUE;
    }

    if (!gote)
    {
        mine = 0;
        maxe = n1 * maxn2;
    }

    if (!gotmr)
    {
        mod = 1;
        res = 0;
    }
    else if (argnum == 5 || argnum > 6)
        badargs = TRUE;

    if (gotd)
    {
        mindeg1 = dval1;
        mindeg2 = dval2;
    }
    if (gotD)
    {
        maxdeg1 = Dval1;
        maxdeg2 = Dval2;
    }
    if (maxdeg1 > maxn2) maxdeg1 = maxn2;
    if (maxdeg2 > n1) maxdeg2 = n1;
    if (connec && mine < n1+maxn2-1) mine = n1 + maxn2 - 1;
    if (connec && mindeg1 == 0) mindeg1 = 1;
    if (connec && mindeg2 == 0) mindeg2 = 1;
    if (maxe > n1*maxdeg1) maxe =  n1*maxdeg1;
    if (maxe > maxn2*maxdeg2) maxe =  maxn2*maxdeg2;
    if (mine < n1*mindeg1) mine = n1*mindeg1;
    if (mine < maxn2*mindeg2) mine = maxn2*mindeg2;

    if (!badargs && (mine > maxe || maxe < 0 || maxdeg1 < 0 || maxdeg2 < 0))
    {
        fprintf(stderr,">E genbg: impossible mine,maxe,maxdeg values\n");
        badargs = TRUE;
    }

    if (!gotZ) maxcommon = -1;

    if (!badargs && (mine > maxe || maxe < 0 || maxdeg1 < 0 || maxdeg2 < 0))
    {
        fprintf(stderr,
                ">E genbg: impossible mine,maxe,mindeg,maxdeg values\n");
        badargs = TRUE;
    }

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if ((nautyformat!=0) + (graph6!=0) + (greout!=0)
                         + (sparse6!=0) + (nooutput!=0) > 1)
        gt_abort(">E genbg: -ungsa are incompatible\n");

#ifdef OUTPROC
    outproc = OUTPROC;
#else
#ifdef OUTPROC1
    outproc = write12;
#endif
    if (nautyformat)   outproc = writenauty;
    else if (nooutput) outproc = nullwrite;
    else if (sparse6)  outproc = writes6x;
    else if (greout)   outproc = writegre;
    else               outproc = writeg6x;
#endif

#ifdef PLUGIN_INIT
PLUGIN_INIT
#endif

    for (i = 0; i <= maxe; ++i) ecount[i] = 0;

    if (nooutput)
        outfile = stdout;
    else if (!gotf || outfilename == NULL)
    {
        outfilename = "stdout";
        outfile = stdout;
    }
    else if ((outfile = fopen(outfilename,
                    nautyformat ? "wb" : "w")) == NULL)
    {
        fprintf(stderr,
              ">E genbg: can't open %s for writing\n",outfilename);
        gt_abort(NULL);
    }

/*
    if (!quiet)
    {
        fprintf(stderr,">A %s n=%d+%d e=%d:%d d=%d:%d D=%d:%d ",
                       argv[0],n1,maxn2,mine,maxe,
                       mindeg1,mindeg2,maxdeg1,maxdeg2);
        if (simple) fprintf(stderr,"z");
        if (footfree) fprintf(stderr,"F");
        if (connec) fprintf(stderr,"c");
        if (maxcommon >= 0) fprintf(stderr,"Z%d",maxcommon);
        if (mod > 1) fprintf(stderr," class=%d/%d",res,mod);
        fprintf(stderr,"\n");
    }
*/

    if (!quiet)
    {
        msg[0] = '\0';
        if (strlen(argv[0]) > 75)
            fprintf(stderr,">A %s",argv[0]);
        else
            CATMSG1(">A %s",argv[0]);
        CATMSG4(" n=%d+%d e=%d:%d",n1,maxn2,mine,maxe);
        CATMSG4(" d=%d:%d D=%d:%d ",mindeg1,mindeg2,maxdeg1,maxdeg2);
        if (simple) CATMSG0("z");
        if (footfree) CATMSG0("F");
        if (antichain) CATMSG0("A");
        if (connec) CATMSG0("c");
        if (maxcommon >= 0) CATMSG1("Z%d",maxcommon);
        if (cutfree) CATMSG0("L");
        if (mod > 1) CATMSG2(" class=%d/%d",res,mod);
        CATMSG0("\n");
        fputs(msg,stderr);
        fflush(stderr);
    }

    class1size = n1;

    for (i = 0; i < n1; ++i)
    {
        g[i] = 0;
        deg[i] = 0;
    }

    t1 = CPUTIME;

    if (maxn2 == 0)
    {
        if (res == 0)
        {
            ++ecount[0];
            (*outproc)(outfile,g,n1,0);
        }
    }
    else
    {
        makeleveldata();
        curres = res;
        if (mod <= 1)        nprune = 0;
        else if (maxn2 >= 6) nprune = maxn2 - 2;
        else if (maxn2 >= 3) nprune = maxn2 - 1;
        else                 nprune = maxn2;

        if (gotX)
        {
            nprune += splitlevinc;
            if (nprune > maxn2) nprune = maxn2;
            if (nprune < 0) nprune = 0;
        }

        xbnds(0,0,0);
        imin = xstart[data[0].xlb];
        imax = xstart[data[0].xub+1];

        for (i = imin; i < imax; ++i)
            data[0].xorb[i] = -1;

        for (i = data[0].xlb; i <= data[0].xub; ++i)
            data[0].xorb[xstart[i]] = xstart[i];

        genextend(g,0,deg,0,FALSE,data[0].xlb,data[0].xub);
    }
    t2 = CPUTIME;

    nout = 0;
    for (i = 0; i <= maxe; ++i) nout += ecount[i];

    if (verbose)
        for (i = 0; i <= maxe; ++i)
            if (ecount[i] != 0)
            {
                fprintf(stderr,">C " COUNTER_FMT " graphs with %d edges\n",
                     ecount[i],i);
            }

#ifdef INSTRUMENT
    fprintf(stderr,"\n>N node counts\n");
    for (i = 0; i < maxn2; ++i)
        fprintf(stderr," level %2d: %7ld (%ld rigid, %ld fertile)\n",
                        i,nodes[i],rigidnodes[i],fertilenodes[i]);
    fprintf(stderr,">A1 %ld calls to accept1, %ld nauty, %ld succeeded\n",
                    a1calls,a1nauty,a1succs);
    fprintf(stderr,
         ">A2 %ld calls to accept2, %ld nuniq, %ld nauty, %ld succeeded\n",
                    a2calls,a2uniq,a2nauty,a2succs);
    fprintf(stderr,"\n");
#endif

#ifdef SPLITTEST
    fprintf(stderr,">Z %lu splitting cases at level %d; cpu=%3.2f sec\n",
            splitcases,nprune,t2-t1);
#else
#ifdef SUMMARY
    SUMMARY(&nout,t2-t1);
#endif

    if (!quiet)
    {
        fprintf(stderr,">Z " COUNTER_FMT " graphs generated in %3.2f sec\n",
                nout,t2-t1);
    }
#endif

    exit(0);
}
