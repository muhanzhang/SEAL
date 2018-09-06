/* TODO:  insert new timings
 *        add chordal graphs 
 *        add perfect graphs
 *        add complements for ordinary graphs */

/* geng.c  version 2.9; B D McKay, Jan 2016. */

#define USAGE \
"geng [-cCmtfbd#D#] [-uygsnh] [-lvq] \n\
              [-x#X#] n [mine[:maxe]] [res/mod] [file]"

#define HELPTEXT \
" Generate all graphs of a specified class.\n\
\n\
      n    : the number of vertices\n\
 mine:maxe : a range for the number of edges\n\
              #:0 means '# or more' except in the case 0:0\n\
   res/mod : only generate subset res out of subsets 0..mod-1\n\
\n\
     -c    : only write connected graphs\n\
     -C    : only write biconnected graphs\n\
     -t    : only generate triangle-free graphs\n\
     -f    : only generate 4-cycle-free graphs\n\
     -b    : only generate bipartite graphs\n\
                (-t, -f and -b can be used in any combination)\n\
     -m    : save memory at the expense of time (only makes a\n\
                difference in the absence of -b, -t, -f and n <= 28).\n\
     -d#   : a lower bound for the minimum degree\n\
     -D#   : a upper bound for the maximum degree\n\
     -v    : display counts by number of edges\n\
     -l    : canonically label output graphs\n\
\n\
     -u    : do not output any graphs, just generate and count them\n\
     -g    : use graph6 output (default)\n\
     -s    : use sparse6 output\n\
     -y    : use the obsolete y-format instead of graph6 format\n\
     -h    : for graph6 or sparse6 format, write a header too\n\
\n\
     -q    : suppress auxiliary output (except from -v)\n\
\n\
  See program text for much more information.\n"


/*  Parameters:

             n    = the number of vertices (1..min(32,WORDSIZE))
             mine = the minimum number of edges (no bounds if missing)
             maxe = the maximum number of edges (same as mine if missing)
                    0 means "infinity" except in the case "0-0"
             mod, res = a way to restrict the output to a subset.
                        All the graphs in G(n,mine..maxe) are divided into
                        disjoint classes C(0,mod),C(1,mod),...,C(mod-1,mod),
                        of very approximately equal size.
                        Only the class C(res,mod) is written.

                        If the -x or -X switch is used, they must have the 
                        same value for different values of res; otherwise 
                        the partitioning may not be valid.  In this case
                        (-x,-X with constant value), the usual relationships 
                        between modulo classes are obeyed; for example 
                        C(3,4) = C(3,8) union C(7,8).  This is not true
                        if 3/8 and 7/8 are done with -x or -X values
                        different from those used for 3/4.

             file = a name for the output file (stdout if missing or "-")

             All switches can be concatenated or separate.  However, the
             value of -d must be attached to the "d", and similarly for "x".

             -c    : only write connected graphs
             -C    : only write biconnected graphs
             -t    : only generate triangle-free graphs
             -f    : only generate 4-cycle-free graphs
             -b    : only generate bipartite graphs
                        (-t, -f and -b can be used in any combination)
             -m    : save memory at expense of time (only makes a
                        difference in the absence of -b, -t, -f and n <= 30).
             -d<int> : specify an upper bound for the maximum degree.
                     The value of the upper bound must be adjacent to
                     the "d".  Example: -d6
             -v    : display counts by number of edges
             -l    : canonically label output graphs

             -u    : do not output any graphs, just generate and count them
             -g    : use graph6 output (default)
             -s    : use sparse6 output
             -n    : use nauty format instead of graph6 format for output
             -y    : use the obsolete y-format instead of graph6 format
             -h    : for graph6 or sparse6 format, write a header too

             -q    : suppress auxiliary output (except from -v)

             -x<int> : specify a parameter that determines how evenly
                     the res/mod facility splits the graphs into subsets.
                     High values mean more even splitting at slight cost
                     to the total time.  The default is 20*mod, and the
                     the legal minimum is 3*mod.  More information is given 
                     under "res/mod" above.
             -X<lev> : move the initial splitting level higher by <lev>,
                     in order to force more even splitting at the cost
                     of speed.  Default is -X0.  More information is given
                     under "res/mod" above.

Output formats.

  The output format is determined by the mutually exclusive switches
  -u, -n, -y, -g and -s.  The default is -g.

  -u suppresses output of graphs completely.

  -s and -g specify sparse6 and graph6 format, defined elsewhere.
  In this case a header is also written if -h is present.

  If -y is present, graphs will be written in y-format.
  y-format is obsolete and only provided for backwards compatibility.

    Each graph occupies one line with a terminating newline.
    Except for the newline, each byte has the format  01xxxxxx, where
    each "x" represents one bit of data.
    First byte:  xxxxxx is the number of vertices n
    Other ceiling(n(n-1)/12) bytes:  These contain the upper triangle of
    the adjacency matrix in column major order.  That is, the entries
    appear in the order (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),(0,4),... .
    The bits are used in left to right order within each byte.
    Any unused bits on the end are set to zero.

  If -n is present, any output graphs are written in nauty format.

    For a graph of n vertices, the output consists of one int giving
    the number of vertices, and n setwords containing the adjacency
    matrix.  Note that this is system dependent (i.e. don't use it).
    It will not work properly if the output is to stdout and your
    system distinguishes binary and text files.

OUTPROC feature.

   By defining the C preprocessor variable OUTPROC at compile time
   (for Unix the syntax is -DOUTPROC=procname on the cc command),
   geng can be made to call a procedure of your manufacture with
   each output graph instead of writing anything. Your procedure
   needs to have type void and the argument list (FILE *f, graph
   *g, int n). f is a stream open for writing, g is the graph in
   nauty format, and n is the number of vertices. Your procedure
   can be in a separate file so long as it is linked with geng. The
   global variables sparse6, graph6, quiet, nooutput, nautyformat,
   yformat and canonise (all type boolean) can be used to test
   for the presence of the flags -s, -g, -q, -u, -n, -y and -l,
   respectively. If -l is present, the group size and similar
   details can be found in the global variable nauty_stats.

PRUNE feature.

   By defining the C preprocessor variable PRUNE at compile time, geng
   can be made to call
        int PRUNE(graph *g,int n,int maxn) 
   for each intermediate (and final) graph, and reject it if 
   the value returned is nonzero.  The arguments are:

     g      = the graph in nauty format (m=1)
     n      = the number of vertices in g
     maxn   = the number of vertices for output 
              (the value you gave on the command line to geng)

   geng constructs the graph starting with vertex 0, then adding
   vertices 1,2,3,... in that order.  Each graph in the sequence is
   an induced subgraph of all later graphs in the sequence.

   A call is made for all orders from 1 to maxn.  In testing for
   a uniform property (such as a forbidden subgraph or forbidden
   induced subgraph) it might save time to notice that a call to
   PRUNE for n implies that the call for n-1 already passed. 

   For very fast tests, it might be worthwhile using PREPRUNE as
   well.  It has the same meaning but is applied earlier and more
   often.
  
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

   It is possible to call geng from another program instead of using it
   as a stand-alone program.  The main requirement is to change the name
   of the main program to be other than "main".  This is done by defining
   the preprocessor variable GENG_MAIN.  You might also like to define
   OUTPROC to be the name of a procedure to receive the graphs. To call
   the program you need to define an argument list argv[] consistent with
   the usual one; don't forget that argv[0] is the command name and not
   the first argument.  The value of argc is the number of strings in
   argv[]; that is, one more than the number of arguments.  See the
   sample program callgeng.c.

**************************************************************************

Sample performance statistics.

    Here we give some graph counts and execution times on a Linux
    Pentium III running at 550 MHz.  Times are with the -u option
    (generate but don't write); add 3-5 microseconds per graph for
    output to a file.  Add another 0.2-0.3 microseconds per graph 
    if you specify connectivity (-c), or 0.6-0.7 microseconds per
    graph if you specific biconnectivity (-C).


      General Graphs                     C3-free Graphs (-t)

      1              1                    1              1
      2              2                    2              2
      3              4                    3              3
      4             11                    4              7
      5             34                    5             14
      6            156                    6             38
      7           1044                    7            107
      8          12346   0.11 sec         8            410
      9         274668   1.77 sec         9           1897
     10       12005168   1.22 min        10          12172   0.21 sec
     11     1018997864   1.72 hr         11         105071   1.49 sec
     12   165091172592    285 hr         12        1262180   15.9 sec
     13 50502031367952  ~10 years        13       20797002   4.08 min
    These can be done in about half      14      467871369   1.50 hr
    the time by setting the edge limit   15    14232552452   45.6 hr
    half way then adding complements.    16   581460254001   79 days
                                         17 31720840164950


     C4-free Graphs  (-f)              (C3,C4)-free Graphs (-tf)

      1             1                      1            1
      2             2                      2            2 
      3             4                      3            3 
      4             8                      4            6
      5            18                      5           11
      6            44                      6           23
      7           117                      7           48 
      8           351                      8          114 
      9          1230                      9          293
     10          5069   0.11 sec          10          869
     11         25181   0.48 sec          11         2963   0.10 sec
     12        152045   2.67 sec          12        12066   0.36 sec
     13       1116403   18.0 sec          13        58933   1.50 sec
     14       9899865   2.50 min          14       347498   7.76 sec
     15     104980369   25.7 min          15      2455693   50.9 sec
     16    1318017549   5.33 hr           16     20592932   6.79 min
     17   19427531763   82.6 hr           17    202724920   1.11 hr
     18  333964672216   62 days           18   2322206466   12.7 hr
     19 6660282066936                     19  30743624324   168 hr
                                          20 468026657815   110 days

                   Old value was wrong:  18   2142368552  (The program was
                   ok, but somehow we tabulated the answer incorrectly.)


     Bipartite Graphs (-b)            C4-free Bipartite Graphs (-bf)

      1              1                      1            1
      2              2                      2            2
      3              3                      3            3
      4              7                      4            6
      5             13                      5           10
      6             35                      6           21
      7             88                      7           39
      8            303                      8           86
      9           1119                      9          182
     10           5479   0.11 sec          10          440
     11          32303   0.59 sec          11         1074
     12         251135   3.99 sec          12         2941   0.15 sec
     13        2527712   35.1 sec          13         8424   0.43 sec
     14       33985853   7.22 min          14        26720   1.37 sec
     15      611846940   2.05 hr           15        90883   4.30 sec
     16    14864650924   48.9 hr           16       340253   14.9 sec
     17   488222721992   70 days           17      1384567   57.1 sec
     18 21712049275198                     18      6186907   4.01 min
                                           19     30219769   18.4 min
                                           20    161763233   1.57 hr
                                           21    946742190   8.85 hr
                                           22   6054606722   56.2 hr
                                           23  42229136988   16.6 days
                                           24 320741332093   121 days

If you know any more of these counts, please tell me.

**************************************************************************

Hints:

To make all the graphs of order n, without restriction on type,
it is fastest to make them up to binomial(n,2)/2 edges and append
the complement of those with strictly less than binomial(n,2)/2 edges.

If it is necessary to split the computation into pieces, it is more
efficient to use the res/mod feature than to split by numbers of edges.

The memory requirements are exponential in n if no maxdeg is given.
For maxdeg=D, the requirements are still exponential (but smaller)
for general graphs, but O(n^D) for other graphs.  Use -m to reduce
the memory requirements to O(n^D) for general graphs too, at the
cost of a small increase in cpu time.

**************************************************************************

    Author:   B. D. McKay, Sep 1991 and many later dates.
              Copyright  B. McKay (1991-2013).  All rights reserved.
              This software is subject to the conditions and waivers
              detailed in the file nauty.h.

    Changes:  Nov 18, 1991 : added -d switch
                             fixed operation for n=16
              Nov 26, 1991 : added OUTPROC feature
              Nov 29, 1991 : -c implies mine >= n-1
              Jan  8, 1992 : make writeny() not static
              Jan 10, 1992 : added -n switch
              Feb  9, 1992 : fixed case of n=1
              Feb 16, 1992 : changed mine,maxe,maxdeg testing
              Feb 19, 1992 : added -b, -t and -u options
                             documented OUTPROC and added external
                                 declaration for it.
              Feb 20, 1992 : added -v option
              Feb 22, 1992 : added INSTRUMENT compile-time option
              Feb 23, 1992 : added xbnds() for more effective pruning
              Feb 24, 1992 : added -l option
              Feb 25, 1992 : changed writenauty() to use fwrite()
              Mar 11, 1992 : completely revised many parts, incl
                             new refinement procedure for fast rejection,
                             distance invariant for regular graphs
              May 19, 1992 : modified userautomproc slightly.  xorb[]
                             is no longer idempotent but it doesn't matter.
                             Speed-up of 2-5% achieved.
              June 5, 1993 : removed ";" after "CPUDEFS" to avoid illegal
                             empty declaration.
              Nov 24, 1994 : tested for 0 <= res < mod

              Apr 13, 1996 : Major overhaul.  Renamed "geng".
                             Changed argument syntax.
                             Removed 16-vertex limit.
                             Added -s, -m, -x.  Allowed combinations.
                             Replaced code for non-general graphs.
                             Very many small changes.
              Jul 12, 1996 : Changed semantics of -x and res/mod.
                             Changed >A line and added fflush()/
                             All switches can be concatenated or not.
              Aug 16, 1996 : Added -X switch and PRUNE() feature.
                             Fixed case of argument 0-0.
              Sep 22, 1996 : Improved 1-2% by tweaking refinex().
              Jan 21, 1997 : Renamed to geng.  
                             Changed -s to -f, and added -sghq.
              Sep  7, 1997 : Fixed WORDSIZE=16 problems.
              Sep 22, 1997 : Use "wb" open for nautyformat.
              Jan 26, 1998 : Added SUMMARY feature.
              Mar  4, 1998 : Added -C.
              Mar 12, 1998 : Moved stats to nauty_stats.
              Jan  1, 2000 : Changed -d to -D and added -d.
              Feb 24, 2000 : Raised limit to 32 vertices.
              Mar  3, 2000 : Made some counts into unsigned long.
                             (Includes first arg to SUMMARY.)
              Mar 12, 2000 : Used bigint for counts that may exceed 2^32.
                             Now all counts from very long runs are ok.
              Oct 12, 2000 : Changed maxef[32] to 92 after confirmation
                             from Yang Yuansheng.  The old value of 93 was
                             valid but 92 is slightly more efficient.
              Nov 16, 2000 : Used fuction prototypes.
              Jul 31, 2001 : Added PREPRUNE
               May 7, 2004 : Complete all function prototypes
              Nov 24, 2004 : Force -m for very large sizes
                             Add -bf automatically if generating trees
              Apr 1, 2007  : Write >A in one fputs() to try to reduce
                             mixing of outputs in multi-process pipes.
              Sep 19, 2007 : Force -m for n > 28 regardless of word size.
              Nov 29, 2008 : Slightly improved connectivity testing.
              Mar 3,  2015 : Improve maxdeg tweaking.
              Jan 18, 2016 : Replace bigint by nauty_counter.

**************************************************************************/

#define NAUTY_PGM  1   /* 1 = geng, 2 = genbg, 3 = gentourng */

#ifndef MAXN
#define MAXN 32         /* not more than max(32,WORDSIZE) */
#endif

#if MAXN > 32
 #error "Can't have MAXN greater than 32"
#endif

#define ONE_WORD_SETS
#include "gtools.h"   /* which includes nauty.h and stdio.h */

#if MAXN < 32
typedef int xword;   /* Must be as large as MAXN bits, and
                    must be unsigned if equal to MAXN bits */
#else
typedef unsigned int xword;
#endif

static void (*outproc)(FILE*,graph*,int);

static FILE *outfile;           /* file for output graphs */
static int connec;              /* 1 for -c, 2 for -C, 0 for neither */
static boolean bipartite;       /* presence of -b */
static boolean trianglefree;    /* presence of -t */
static boolean squarefree;      /* presence of -f */
static boolean savemem;         /* presence of -m */
static boolean verbose;         /* presence of -v */
boolean nautyformat;            /* presence of -n */
boolean yformat;                /* presence of -y */
boolean graph6;                 /* presence of -g */
boolean sparse6;                /* presence of -s */
boolean nooutput;               /* presence of -u */
boolean canonise;               /* presence of -l */
boolean quiet;                  /* presence of -q */
boolean header;                 /* presence of -h */
statsblk nauty_stats;
static int mindeg,maxdeg,maxn,mine,maxe,mod,res;
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
    int ne,dmax;          /* values used for xlb,xub calculation */
    int xlb,xub;          /* saved bounds on extension degree */
    xword lo,hi;          /* work purposes for orbit calculation */
    xword xstart[MAXN+1]; /* index into xset[] for each cardinality */
    xword *xset;          /* array of all x-sets in card order */
    xword *xcard;         /* cardinalities of all x-sets */
    xword *xinv;          /* map from x-set to index in xset */
    xword *xorb;          /* min orbit representative */
    xword *xx;            /* (-b, -t, -s, -m) candidate x-sets */
                      /*   note: can be the same as xcard */
    xword xlim;           /* number of x-sets in xx[] */
} leveldata;

static leveldata data[MAXN];      /* data[n] is data for n -> n+1 */
static nauty_counter ecount[1+MAXN*(MAXN-1)/2];  /* counts by number of edges */
static nauty_counter nodes[MAXN];     /* nodes at each level */

#ifdef INSTRUMENT
static unsigned long rigidnodes[MAXN],fertilenodes[MAXN];
static unsigned long a1calls,a1nauty,a1succs;
static unsigned long a2calls,a2nauty,a2uniq,a2succs;
#endif

/* The numbers below are actual maximum edge counts.  The apparently
   anomolous value of 92 for maxef[32] has been confirmed independently
   by Yang Yuansheng (as well as all the smaller maxef[] values).
   geng works correctly with any upper bounds.
   To extend an upper bound upwards: (n-1, E) -> (n, E + floor(2*E/(n-2))).
*/

static int maxeb[] =     /* max edges for -b */
 {0,0,1,2,4, 6,9,12,16,20, 25,30,36,42,49, 56,64,72,81,90,
  100,110,121,132,144, 156,169,182,196,210, 225,240,256};
static int maxet[] =     /* max edges for -t */
 {0,0,1,2,4, 6,9,12,16,20, 25,30,36,42,49, 56,64,72,81,90,
  100,110,121,132,144, 156,169,182,196,210, 225,240,256};
static int maxef[] =     /* max edges for -f */
 {0,0,1,3,4, 6,7,9,11,13,  16,18,21,24,27, 30,33,36,39,42,
  46,50,52,56,59, 63,67,71,76,80, 85,90,92};
static int maxeft[] =    /* max edges for -ft */
 {0,0,1,2,3, 5,6,8,10,12,  15,16,18,21,23, 26,28,31,34,38,
  41,44,47,50,54, 57,61,65,68,72, 76,80,85};
static int maxebf[] =    /* max edges for -bf */
 {0,0,1,2,3, 4,6,7,9,10,   12,14,16,18,21, 22,24,26,29,31,
  34,36,39,42,45, 48,52,53,56,58, 61,64,67}; 

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
writeny(FILE *f, graph *g, int n)
/* write graph g (n vertices) to file f in y format */
{
    static char ybit[] = {32,16,8,4,2,1};
    char s[(MAXN*(MAXN-1)/2 + 5)/6 + 4];
    int i,j,k;
    char y,*sp;

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

/***********************************************************************/

static void
nullwrite(FILE *f, graph *g, int n)
/* don't write graph g (n vertices) to file f */
{
}

/***********************************************************************/

void
writenauty(FILE *f, graph *g, int n)
/* write graph g (n vertices) to file f in nauty format */
{
    int nn;

    nn = n;

    if (fwrite((char *)&nn,sizeof(int),(size_t)1,f) != 1 ||
          fwrite((char*)g,sizeof(setword),(size_t)n,f) != n)
    {
        fprintf(stderr,">E writenauty : error on writing file\n");
        exit(2);
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

/**********************************************************************/
 
static boolean
isbiconnected(graph *g, int n)
/* test if g is biconnected */
{
    int sp,v,w;
    setword sw;
    setword visited;
    int numvis,num[MAXN],lp[MAXN],stack[MAXN];
 
    if (n <= 2) return FALSE;
 
    visited = bit[0];
    stack[0] = 0;
    num[0] = 0;
    lp[0] = 0;
    numvis = 1;
    sp = 0;
    v = 0;
 
    for (;;)
    {
        if ((sw = g[v] & ~visited))           /* not "==" */
        {
            w = v;
            v = FIRSTBITNZ(sw);       /* visit next child */
            stack[++sp] = v;
            visited |= bit[v];
            lp[v] = num[v] = numvis++;
            sw = g[v] & visited & ~bit[w];
            while (sw)
            {
                w = FIRSTBITNZ(sw);
                sw &= ~bit[w];
                if (num[w] < lp[v])  lp[v] = num[w];
            }
        }
        else
        {
            w = v;                  /* back up to parent */
            if (sp <= 1)          return numvis == n;
            v = stack[--sp];
            if (lp[w] >= num[v])  return FALSE;
            if (lp[w] < lp[v])    lp[v] = lp[w];
        }
    }
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

/**********************************************************************/

static boolean
distinvar(graph *g, int *invar, int n)
/* make distance invariant
   return FALSE if n-1 not maximal else return TRUE */
{
    int w;
    setword workset,frontier;
    setword sofar;
    int inv,d,v;

    for (v = n-1; v >= 0; --v)
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
                frontier ^= bit[w];
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
makexgraph(graph *g, xword *h, int n)
/* make x-format graph from nauty format graph */
{
    setword gi;
    int i,j;
    xword hi;

    for (i = 0; i < n; ++i)
    {
        hi = 0;
        gi = g[i];
        while (gi)
        {
            j = FIRSTBITNZ(gi);
            gi ^= bit[j];
            hi |= xbit[j];
        }
        h[i] = hi;
    }
}

/**************************************************************************/

static void
make0graph(graph *g, xword *h, int n)
/* make x-format graph without edges */
{
    int i;

    for (i = 0; i < n; ++i) h[i] = 0;
}

/**************************************************************************/

static void
makebgraph(graph *g, xword *h, int n)
/* make x-format graph of different colour graph */
{
    setword seen1,seen2,expanded,w;
    setword restv;
    xword xseen1,xseen2;
    int i;

    restv = 0;
    for (i = 0; i < n; ++i) restv |= bit[i];

    seen1 = seen2 = 0;
    expanded = 0;

    while (TRUE)
    {
        if ((w = ((seen1 | seen2) & ~expanded)) == 0)
        {
            xseen1 = 0;
            w = seen1;
            while (w)
            {
                i = FIRSTBITNZ(w);
                w ^= bit[i];
                xseen1 |= xbit[i];
            }
            xseen2 = 0;
            w = seen2;
            while (w)
            {
                i = FIRSTBITNZ(w);
                w ^= bit[i];
                xseen2 |= xbit[i];
            }

            w = seen1;
            while (w)
            {
                i = FIRSTBITNZ(w);
                w ^= bit[i];
                h[i] = xseen2;
            }
            w = seen2;
            while (w)
            {
                i = FIRSTBITNZ(w);
                w ^= bit[i];
                h[i] = xseen1;
            }

            restv &= ~(seen1 | seen2);
            if (restv == 0) return;
            i = FIRSTBITNZ(restv);
            seen1 = bit[i];
            seen2 = 0;
        }
        else
            i = FIRSTBITNZ(w);

        expanded |= bit[i];
        if (bit[i] & seen1) seen2 |= g[i];
        else                seen1 |= g[i];
    }
}

/**************************************************************************/
 
static void
makeb6graph(graph *g, xword *h, int n)
/* make x-format bipartite girth 6 graph */
{
    setword w,x;
    xword hi;
    int i,j;

    makebgraph(g,h,n);

    for (i = 0; i < n; ++i)
    {
        w = g[i];
        x = 0;
        while (w)
        {
            j = FIRSTBITNZ(w);
            w ^= bit[j];
            x |= g[j];
        }
        x &= ~bit[i];
        hi = h[i];
        while (x)
        {
            j = FIRSTBITNZ(x);
            x ^= bit[j];
            hi |= xbit[j];
        }
        h[i] = hi;
    }
}

/**************************************************************************/
 
static void
makesgraph(graph *g, xword *h, int n)
/* make x-format square graph */
{
    setword w,x;
    xword hi;
    int i,j;

    for (i = 0; i < n; ++i)
    {
        w = g[i];
        x = 0;
        while (w)
        {
            j = FIRSTBITNZ(w);
            w ^= bit[j];
            x |= g[j];
        }
        x &= ~bit[i];
        hi = 0;
        while (x)
        {
            j = FIRSTBITNZ(x);
            x ^= bit[j];
            hi |= xbit[j];
        }
        h[i] = hi;
    }
}

/**************************************************************************/ 
 
static void 
makeg5graph(graph *g, xword *h, int n)
/* make x-format girth-5 graph */
{
    setword w,x; 
    xword hi;
    int i,j;
 
    for (i = 0; i < n; ++i)
    { 
        w = g[i]; 
        x = g[i];
        while (w) 
        {
            j = FIRSTBITNZ(w);
            w ^= bit[j];
            x |= g[j];
        } 
        x &= ~bit[i]; 
        hi = 0; 
        while (x) 
        { 
            j = FIRSTBITNZ(x); 
            x ^= bit[j]; 
            hi |= xbit[j]; 
        } 
        h[i] = hi; 
    } 
} 

/**************************************************************************/  

static void
makeleveldata(boolean restricted)
/* make the level data for each level */
{
    long h;
    int n,nn;
    long ncj;
    leveldata *d;
    xword *xcard,*xinv;
    xword *xset,xw,tttn,nxsets;
    xword cw;
    xword i,j;

    for (n = 1; n < maxn; ++n)
    {
        nn = maxdeg <= n ? maxdeg : n;
        ncj = nxsets = 1;
        for (j = 1; j <= nn; ++j)
        {
            ncj = (ncj * (n - j + 1)) / j;
            nxsets += ncj;
        }
        tttn = 1L << n;

        d = &data[n];

        d->ne = d->dmax = d->xlb = d->xub = -1;

        if (restricted)
        {
            d->xorb = (xword*) calloc(nxsets,sizeof(xword));
            d->xx = (xword*) calloc(nxsets,sizeof(xword));
            if (d->xorb == NULL || d->xx == NULL)
            {
                fprintf(stderr,
                   ">E geng: calloc failed in makeleveldata()\n");
                exit(2);
            }
            continue;   /* <--- NOTE THIS! */
        }

        d->xset = xset = (xword*) calloc(nxsets,sizeof(xword));
        d->xcard = xcard = (xword*) calloc(nxsets,sizeof(xword));
        d->xinv = xinv = (xword*) calloc(tttn,sizeof(xword));
        d->xorb = (xword*) calloc(nxsets,sizeof(xword));
        d->xx = d->xcard;

        if (xset==NULL || xcard==NULL || xinv==NULL || d->xorb==NULL)
        {
            fprintf(stderr,">E geng: calloc failed in makeleveldata()\n");
            exit(2);
        }

        j = 0;

        for (i = 0;; ++i)
        {
            if ((h = XPOPCOUNT(i)) <= maxdeg)
            {
                xset[j] = i;
                xcard[j] = h;
                ++j;
            }
            if (i == (xword)((1L<<n)-1)) break;
        }

        if (j != nxsets)
        {
            fprintf(stderr,">E geng: j=%u mxsets=%u\n",
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

/**************************************************************************/

static void
userautomprocb(int count, int *p, int *orbits,
          int numorbits, int stabvertex, int n)
/* form orbits on powerset of VG
   called by nauty;  operates on data[n] */
{
    xword j1,j2,moved,pi,pxi,lo,hi,x;
    xword i,*xorb,*xx,w,xlim,xlb;

    xorb = data[n].xorb;
    xx = data[n].xx;
    xlim = data[n].xlim;

    if (count == 1)                         /* first automorphism */
    {
        j1 = 0;
        xlb = data[n].xlb;

        for (i = 0; i < xlim; ++i)
        {
            x = xx[i];
            if (XPOPCOUNT(x) >= xlb)
            {
                xx[j1] = x;
                xorb[j1] = j1;
                ++j1;
            }
        }
        data[n].xlim = xlim = j1;
    }

    moved = 0;
    for (i = 0; i < n; ++i)
        if (p[i] != i) moved |= xbit[i];

    for (i = 0; i < xlim; ++i)
    {
        if ((w = xx[i] & moved) == 0) continue;
        pxi = xx[i] & ~moved;
        while (w)
        {
            j1 = XNEXTBIT(w);
            w ^= xbit[j1];
            pxi |= xbit[p[j1]];
        }
        /* pi = position of pxi */

        lo = 0;
        hi = xlim - 1;

        for (;;)
        {
            pi = (lo + hi) >> 1;
            if (xx[pi] == pxi) break;
            else if (xx[pi] < pxi) lo = pi + 1;
            else                   hi = pi - 1;
        }

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
    xword xw;
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

    xw = x;
    while (xw)
    {
        i = XNEXTBIT(xw);
        xw ^= xbit[i];
        gx[i] |= bit[n];
        gx[n] |= bit[i];
        ++deg[i];
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
        return TRUE;
    }
    else
        return FALSE;
}

/**************************************************************************/

static boolean
accept1b(graph *g, int n, xword x, graph *gx, int *deg, boolean *rigid,
     void (*makeh)(graph*,xword*,int))
/* decide if n in theta(g+x)  --  version for n+1 < maxn */
{
    int i,v;
    xword z,hv,bitv,ixx;
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    int count[MAXN];
    graph gc[MAXN];
    xword h[MAXN],xw,jxx,kxx,*xx;
    int nx,numcells,code;
    int i0,i1,degn,xubx;
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

    xw = x;
    while (xw)
    {
        i = XNEXTBIT(xw);
        xw ^= xbit[i];
        gx[i] |= bit[n];
        gx[n] |= bit[i];
        ++deg[i];
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

    (*makeh)(gx,h,nx);
    xx = data[nx].xx;
    xubx = data[nx].xub;

    xx[0] = 0;
    kxx = 1;
    for (v = 0; v < nx; ++v)
    {
        bitv = xbit[v];
        hv = h[v];
        jxx = kxx;
        for (ixx = 0; ixx < jxx; ++ixx)
        if ((hv & xx[ixx]) == 0)
        {
            z = xx[ixx] | bitv;
            if (XPOPCOUNT(z) <= xubx) xx[kxx++] = z;
        }
    }
    data[nx].xlim = kxx;

    if (numcells == nx)
    {
        *rigid = TRUE;
#ifdef INSTRUMENT
        ++a1succs;
#endif
        return TRUE;
    }

    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    options.userautomproc = userautomprocb;

    active[0] = 0;
#ifdef INSTRUMENT
    ++a1nauty;
#endif
    nauty(gx,lab,ptn,active,orbits,&options,&stats,workspace,50,1,nx,gc);

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
accept2(graph *g, int n, xword x, graph *gx, int *deg, boolean nuniq)
/* decide if n in theta(g+x)  --  version for n+1 == maxn */
{
    int i;
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    int degx[MAXN],invar[MAXN];
    setword vmax,gv;
    int qn,qv;
    int count[MAXN];
    xword xw;
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
        xw ^= xbit[i];
        gx[i] |= bit[n];
        gx[n] |= bit[i];
        ++degx[i];
    }

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

        if (!distinvar(gx,invar,nx)) return FALSE;
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
    if (code > 0 || numcells >= nx-4)
        cheapacc = TRUE;
    else if (numcells == nx-5)
    {
        for (j1 = nx-2; j1 >= 0 && ptn[j1] > 0; --j1) {}
        if (nx - j1 != 5) cheapacc = TRUE;
    }
    else
    {
        j1 = nx;
        j0 = 0;
        for (i1 = 0; i1 < nx; ++i1)
        {
            --j1;
            if (ptn[i1] > 0)
            {
                ++j0;
                while (ptn[++i1] > 0) {}
            }
        }
        if (j1 <= j0 + 1) cheapacc = TRUE;
    }
    
    if (cheapacc)
    {
#ifdef INSTRUMENT
        ++a2succs;
#endif
        if (canonise) makecanon(gx,gcan,nx);
        return TRUE;
    }

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
        if (canonise) makecanon(gx,gcan,nx);
        return TRUE;
    }
    else
        return FALSE;
}

/**************************************************************************/

static void
xbnds(int n, int ne, int dmax)
/* find bounds on extension degree;  store answer in data[*].*  */
{
    int xlb,xub,d,nn,m,xc;

    xlb = n == 1 ? 0 : (dmax > (2*ne + n - 2)/(n - 1) ?
                        dmax : (2*ne + n - 2)/(n - 1));
    xub = n < maxdeg ? n : maxdeg;

    for (xc = xub; xc >= xlb; --xc)
    {
        d = xc;
        m = ne + d;
        for (nn = n+1; nn < maxn; ++nn)
        {
            if (d < (2*m + nn - 2)/(nn - 1)) d = (2*m + nn - 2)/(nn - 1);
            m += d;
        }
        if (d > maxdeg || m > maxe) xub = xc - 1;
        else                        break;
    }

    if (ne + xlb < mine)
        for (xc = xlb; xc <= xub; ++xc)
        {
            m = ne + xc;
            for (nn = n + 1; nn < maxn; ++nn)
                m += maxdeg < nn ? maxdeg : nn;
            if (m < mine) xlb = xc + 1;
            else          break;
        }

    data[n].ne = ne;
    data[n].dmax = dmax;
    data[n].xlb = xlb;
    data[n].xub = xub;
}

/**************************************************************************/

static void
spaextend(graph *g, int n, int *deg, int ne, boolean rigid,
      int xlb, int xub, void (*makeh)(graph*,xword*,int))
/* extend from n to n+1 -- version for restricted graphs */
{
    xword x,d,dlow;
    xword xlim,*xorb;
    int xc,nx,i,j,dmax,dcrit,xlbx,xubx;
    graph gx[MAXN];
    xword *xx,ixx;
    int degx[MAXN];
    boolean rigidx;

#ifdef INSTRUMENT
    boolean haschild;

    haschild = FALSE;
    if (rigid) ++rigidnodes[n];
#endif
    ++nodes[n];

    nx = n + 1;
    dmax = deg[n-1];
    dcrit = mindeg - maxn + n;
    d = dlow = 0;
    for (i = 0; i < n; ++i)
    {
        if (deg[i] == dmax) d |= xbit[i];
        if (deg[i] == dcrit) dlow |= xbit[i];
    }

    if (xlb == dmax && XPOPCOUNT(d) + dmax > n) ++xlb;
    if (nx == maxn && xlb < mindeg) xlb = mindeg;
    if (xlb > xub) return;

#ifdef PRUNE
    if (PRUNE(g,n,maxn)) return;
#endif

    xorb = data[n].xorb;
    xx = data[n].xx;
    xlim = data[n].xlim;

    if (nx == maxn)
    {
        for (ixx = 0; ixx < xlim; ++ixx)
        {
            x = xx[ixx];
            xc = XPOPCOUNT(x);
            if (xc < xlb || xc > xub) continue;
            if ((rigid || xorb[ixx] == ixx) 
                && (xc > dmax || (xc == dmax && (x & d) == 0))
                && (dlow & ~x) == 0)
            {
                if (accept2(g,n,x,gx,deg,
                    xc > dmax+1 || (xc == dmax+1 && (x & d) == 0))
                    && (!connec ||
                          (connec==1 && isconnected(gx,nx)) ||
                          (connec>1 && isbiconnected(gx,nx))))
                {
#ifdef PRUNE
                    if (!PRUNE(gx,nx,maxn))
#endif
                    {
#ifdef INSTRUMENT
                        haschild = TRUE;
#endif
                        ++ecount[ne+xc];
                        (*outproc)(outfile,canonise ? gcan : gx,nx);
                    }
                }
            }
        }
    }
    else
    {
        for (ixx = 0; ixx < xlim; ++ixx)
        {
            if (nx == splitlevel)
            {
                if (odometer-- != 0) continue;
                odometer = mod - 1;
            }
            x = xx[ixx];
            xc = XPOPCOUNT(x);
            if (xc < xlb || xc > xub) continue;
            if ((rigid || xorb[ixx] == ixx)
                && (xc > dmax || (xc == dmax && (x & d) == 0))
                && (dlow & ~x) == 0)
            {
                for (j = 0; j < n; ++j) degx[j] = deg[j];
                if (data[nx].ne != ne+xc || data[nx].dmax != xc)
                    xbnds(nx,ne+xc,xc);

                xlbx = data[nx].xlb;
                xubx = data[nx].xub;
                if (xlbx <= xubx
                    && accept1b(g,n,x,gx,degx,&rigidx,makeh))
                {
#ifdef INSTRUMENT
                    haschild = TRUE;
#endif
                    spaextend(gx,nx,degx,ne+xc,rigidx,xlbx,xubx,makeh);
                }
            }
        }
        if (n == splitlevel - 1 && n >= min_splitlevel
                                && nodes[n] >= multiplicity)
            --splitlevel;
    }
#ifdef INSTRUMENT
    if (haschild) ++fertilenodes[n];
#endif
}

/**************************************************************************/

static void
genextend(graph *g, int n, int *deg, int ne, boolean rigid, int xlb, int xub)
/* extend from n to n+1 -- version for general graphs */
{
    xword x,d,dlow;
    xword *xset,*xcard,*xorb;
    xword i,imin,imax;
    int nx,xc,j,dmax,dcrit;
    int xlbx,xubx;
    graph gx[MAXN];
    int degx[MAXN];
    boolean rigidx;

#ifdef INSTRUMENT
    boolean haschild;

    haschild = FALSE;
    if (rigid) ++rigidnodes[n];
#endif
    ++nodes[n];

    nx = n + 1;
    dmax = deg[n-1];
    dcrit = mindeg - maxn + n;
    d = dlow = 0;
    for (i = 0; i < n; ++i)
    {
        if (deg[i] == dmax) d |= xbit[i];
        if (deg[i] == dcrit) dlow |= xbit[i];
    }

    if (xlb == dmax && XPOPCOUNT(d) + dmax > n) ++xlb;
    if (nx == maxn && xlb < mindeg) xlb = mindeg;
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
        for (i = imin; i < imax; ++i)
        {
            if (!rigid && xorb[i] != i) continue;
            x = xset[i];
            xc = xcard[i];
            if (xc == dmax && (x & d) != 0) continue;
            if ((dlow & ~x) != 0) continue;

            if (accept2(g,n,x,gx,deg,
                        xc > dmax+1 || (xc == dmax+1 && (x & d) == 0)))
                if (!connec || (connec==1 && isconnected(gx,nx))
                            || (connec>1 && isbiconnected(gx,nx)))
                {
#ifdef PRUNE
                    if (!PRUNE(gx,nx,maxn))
#endif
                    {
#ifdef INSTRUMENT
                        haschild = TRUE;
#endif
                        ++ecount[ne+xc];
                        (*outproc)(outfile,canonise ? gcan : gx,nx);
                    }
                }
        }
    else
        for (i = imin; i < imax; ++i)
        {
            if (!rigid && xorb[i] != i) continue;
            x = xset[i];
            xc = xcard[i];
            if (xc == dmax && (x & d) != 0) continue;
            if ((dlow & ~x) != 0) continue;
            if (nx == splitlevel)
            {
                if (odometer-- != 0) continue;
                odometer = mod - 1;
            }

            for (j = 0; j < n; ++j) degx[j] = deg[j];
            if (data[nx].ne != ne+xc || data[nx].dmax != xc)
                xbnds(nx,ne+xc,xc);
            xlbx = data[nx].xlb;
            xubx = data[nx].xub;
            if (xlbx > xubx) continue;

            data[nx].lo = data[nx].xstart[xlbx];
            data[nx].hi = data[nx].xstart[xubx+1];
            if (accept1(g,n,x,gx,degx,&rigidx))
            {
#ifdef INSTRUMENT
                haschild = TRUE;
#endif
                genextend(gx,nx,degx,ne+xc,rigidx,xlbx,xubx);
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

#ifdef GENG_MAIN
int
GENG_MAIN(int argc, char *argv[])
#else
int
main(int argc, char *argv[])
#endif
{
    char *arg;
    boolean badargs,gote,gotmr,gotf,gotd,gotD,gotx,gotX;
    boolean secret,connec1,connec2,safe,sparse;
    char *outfilename,sw;
    int i,j,argnum;
    graph g[1];
    int tmaxe,deg[1];
    nauty_counter nout;
    int splitlevinc;
    xword testxword;
    double t1,t2;
    char msg[201];

    HELP; PUTVERSION;
    nauty_check(WORDSIZE,1,MAXN,NAUTYVERSIONID);

    testxword = (xword)(-1);
    if (MAXN > 32 || MAXN > WORDSIZE || MAXN > 8*sizeof(xword)
        || (MAXN == 8*sizeof(xword) && testxword < 1))
    {
        fprintf(stderr,"geng: incompatible MAXN, WORDSIZE, or xword\n");
        fprintf(stderr,"--See notes in program source\n");
        exit(1);
    }

    badargs = FALSE;
    trianglefree = FALSE;
    bipartite = FALSE;
    squarefree = FALSE;
    verbose = FALSE;
    nautyformat = FALSE;
    yformat = FALSE;
    graph6 = FALSE;
    sparse6 = FALSE;
    savemem = FALSE;
    nooutput = FALSE;
    canonise = FALSE;
    header = FALSE;
    outfilename = NULL;
    secret = FALSE;
    safe = FALSE;
    connec1 = connec2 = FALSE;

    maxdeg = MAXN; 
    mindeg = 0;
    
    gotX = gotx = gotd = gotD = gote = gotmr = gotf = FALSE;

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
                else SWBOOLEAN('g',graph6)
                else SWBOOLEAN('s',sparse6)
                else SWBOOLEAN('t',trianglefree)
                else SWBOOLEAN('f',squarefree)
                else SWBOOLEAN('b',bipartite)
                else SWBOOLEAN('v',verbose)
                else SWBOOLEAN('l',canonise)
                else SWBOOLEAN('y',yformat)
                else SWBOOLEAN('h',header)
                else SWBOOLEAN('m',savemem)
                else SWBOOLEAN('c',connec1)
                else SWBOOLEAN('C',connec2)
                else SWBOOLEAN('q',quiet)
                else SWBOOLEAN('$',secret)
                else SWBOOLEAN('S',safe)
                else SWINT('d',gotd,mindeg,"geng -d")
                else SWINT('D',gotD,maxdeg,"geng -D")
                else SWINT('x',gotx,multiplicity,"geng -x")
                else SWINT('X',gotX,splitlevinc,"geng -X")
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
                if (!gote)
                {
                    if (sscanf(arg,"%d:%d",&mine,&maxe) == 2
                     || sscanf(arg,"%d-%d",&mine,&maxe) == 2)
                    {
                        gote = TRUE;
                        if (maxe == 0 && mine > 0) maxe = MAXN*(MAXN-1)/2;
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

    if (argnum == 0)
        badargs = TRUE;
    else if (maxn < 1 || maxn > MAXN)
    {
        fprintf(stderr,">E geng: n must be in the range 1..%d\n",MAXN);
        badargs = TRUE;
    }

    if (!gotmr)
    {
        mod = 1;
        res = 0;
    }

    if (!gote)
    {
        mine = 0;
        maxe = (maxn*maxn - maxn) / 2;
    }

    if (connec1 && mindeg < 1 && maxn > 1) mindeg = 1;
    if (connec2 && mindeg < 2 && maxn > 2) mindeg = 2;
    if (maxdeg >= maxn) maxdeg = maxn - 1;
    if (maxe > maxn*maxdeg / 2) maxe = maxn*maxdeg / 2;
    if (maxdeg > maxe) maxdeg = maxe;
    if (mindeg < 0) mindeg = 0;
    if (mine < (maxn*mindeg+1) / 2) mine = (maxn*mindeg+1) / 2;
    if (maxdeg > 2*maxe - mindeg*(maxn-1)) maxdeg = 2*maxe - mindeg*(maxn-1);

    if (!badargs && (mine > maxe || maxe < 0 || maxdeg < 0))
    {
        fprintf(stderr,
                ">E geng: impossible mine,maxe,mindeg,maxdeg values\n");
        badargs = TRUE;
    }

    if (!badargs && (res < 0 || res >= mod))
    {
        fprintf(stderr,">E geng: must have 0 <= res < mod\n");
        badargs = TRUE;
    }

    if      (connec2) connec = 2;
    else if (connec1) connec = 1;
    else              connec = 0;

    if (connec && mine < maxn-1) mine = maxn - 2 + connec;

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if ((nautyformat!=0) + (yformat!=0) + (graph6!=0)
                         + (sparse6!=0) + (nooutput!=0) > 1)
        gt_abort(">E geng: -uyngs are incompatible\n");

#ifdef OUTPROC
    outproc = OUTPROC;
#else
    if (nautyformat)   outproc = writenauty;
    else if (yformat)  outproc = writeny;
    else if (nooutput) outproc = nullwrite;
    else if (sparse6)  outproc = writes6x;
    else               outproc = writeg6x;
#endif

#ifdef PLUGIN_INIT
PLUGIN_INIT
#endif

    for (i = 0; i <= maxe; ++i) ecount[i] = 0;
    for (i = 0; i < maxn; ++i)  nodes[i] = 0;

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
              ">E geng: can't open %s for writing\n",outfilename);
        gt_abort(NULL);
    }

    if (bipartite)
        if (squarefree)  tmaxe = maxebf[maxn];
        else             tmaxe = maxeb[maxn];
    else if (trianglefree)
        if (squarefree)  tmaxe = maxeft[maxn];
        else             tmaxe = maxet[maxn];
    else if (squarefree) tmaxe = maxef[maxn];
    else                 tmaxe = (maxn*maxn - maxn) / 2;

    if (safe) ++tmaxe;

    if (maxe > tmaxe) maxe = tmaxe;

    if (gotx)
    {
        if (multiplicity < 3 * mod || multiplicity > 999999999)
            gt_abort(">E geng: -x value must be in [3*mod,10^9-1]\n");
    }
    else
    {
        multiplicity = PRUNEMULT * mod;
	if (multiplicity / PRUNEMULT != mod)
	    gt_abort(">E geng: mod value is too large\n");
    }

    if (!gotX) splitlevinc = 0;

    if (!quiet)
    {
        msg[0] = '\0';
        if (strlen(argv[0]) > 75)
            fprintf(stderr,">A %s",argv[0]);
        else
            CATMSG1(">A %s",argv[0]);
       
        CATMSG6(" -%s%s%s%s%s%s",
            connec2      ? "C" : connec1 ? "c" : "",
            trianglefree ? "t" : "",
            squarefree   ? "f" : "",
            bipartite    ? "b" : "",
            canonise     ? "l" : "",
            savemem      ? "m" : "");
        if (mod > 1)
            CATMSG2("X%dx%d",splitlevinc,multiplicity);
        CATMSG4("d%dD%d n=%d e=%d",mindeg,maxdeg,maxn,mine);
        if (maxe > mine) CATMSG1("-%d",maxe);
        if (mod > 1) CATMSG2(" class=%d/%d",res,mod);
        CATMSG0("\n");
        fputs(msg,stderr);
        fflush(stderr);
    }

    g[0] = 0;
    deg[0] = 0;

    sparse = bipartite || squarefree || trianglefree || savemem;

    t1 = CPUTIME;

    if (header)
    {
        if (sparse6)
        {
            writeline(outfile,SPARSE6_HEADER);
            fflush(outfile);
        }
        else if (!yformat && !nautyformat && !nooutput)
        {
            writeline(outfile,GRAPH6_HEADER);
            fflush(outfile);
        }
    }

    if (maxn == 1)
    {
        if (res == 0)
        {
            ++ecount[0];
            (*outproc)(outfile,g,1);
        }
    }
    else
    {
        if (maxn > 28 || maxn+4 > 8*sizeof(xword))
            savemem = sparse = TRUE;
        if (maxn == maxe+1 && connec)
            bipartite = squarefree = sparse = TRUE;  /* trees */

        makeleveldata(sparse);

        if (maxn >= 14 && mod > 1)     splitlevel = maxn - 4;
        else if (maxn >= 6 && mod > 1) splitlevel = maxn - 3;
        else                           splitlevel = -1;

        if (splitlevel > 0) splitlevel += splitlevinc;
        if (splitlevel > maxn - 1) splitlevel = maxn - 1;
        if (splitlevel < 3) splitlevel = -1;

        min_splitlevel = 6;
        odometer = secret ? -1 : res;

        if (maxe >= mine && 
                (mod <= 1 || (mod > 1 && (splitlevel > 2 || res == 0))))
        {
            xbnds(1,0,0);   
            if (sparse)
            {
                data[1].xx[0] = 0;
                if (maxdeg > 0) data[1].xx[1] = xbit[0];
                data[1].xlim = data[1].xub + 1;
            }

            if (bipartite)
                if (squarefree)
                    spaextend(g,1,deg,0,TRUE,
                                    data[1].xlb,data[1].xub,makeb6graph);
                else
                    spaextend(g,1,deg,0,TRUE,
                                    data[1].xlb,data[1].xub,makebgraph);
            else if (trianglefree)
                if (squarefree)
                    spaextend(g,1,deg,0,TRUE,
                                    data[1].xlb,data[1].xub,makeg5graph);
                else
                    spaextend(g,1,deg,0,TRUE,
                                    data[1].xlb,data[1].xub,makexgraph);
            else if (squarefree)
                spaextend(g,1,deg,0,TRUE,
                                    data[1].xlb,data[1].xub,makesgraph);
            else if (savemem)
                spaextend(g,1,deg,0,TRUE,
                                    data[1].xlb,data[1].xub,make0graph);
            else
                genextend(g,1,deg,0,TRUE,data[1].xlb,data[1].xub);
        }
    }
    t2 = CPUTIME;

    nout = 0;
    for (i = 0; i <= maxe; ++i) nout += ecount[i];

    if (verbose)
    {
        for (i = 0; i <= maxe; ++i)
            if (ecount[i] != 0)
            {
		fprintf(stderr,">C " COUNTER_FMT " graphs with %d edges\n",
                     ecount[i],i);
            }
    }

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
        if (sparse)
        {
            free(data[i].xorb);
            free(data[i].xx);
        }
        else
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
