/*****************************************************************************
*                                                                            *
*  Auxiliary source file for version 2.6 of nauty.                           *
*                                                                            *
*   Copyright (1984-2013) Brendan McKay.  All rights reserved.               *
*   Subject to waivers and disclaimers in nauty.h.                           *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       10-Nov-87 : final changes for version 1.2                            *
*        5-Dec-87 : renamed to version 1.3 (no changes to this file)         *
*       28-Sep-88 : renamed to version 1.4 (no changes to this file)         *
*       23-Mar-89 : changes for version 1.5 :                                *
*                   - added procedure refine1()                              *
*                   - changed type of ptn from int* to nvector* in fmptn()   *
*                   - declared level in breakout()                           *
*                   - changed char[] to char* in a few places                *
*                   - minor rearrangement in bestcell()                      *
*       31-Mar-89 : - added procedure doref()                                *
*        5-Apr-89 : - changed MAKEEMPTY uses to EMPTYSET                     *
*       12-Apr-89 : - changed writeperm() and fmperm() to not use MARKing    *
*        5-May-89 : - redefined MASH to gain about 8% efficiency             *
*       18-Oct-90 : changes for version 1.6 :                                *
*                   - improved line breaking in writeperm()                  *
*       10-Nov-90 : - added dummy routine nautil_null()                      *
*       27-Aug-92 : changes for version 1.7 :                                *
*                   - made linelength <= 0 mean no line breaks               *
*        5-Jun-93 : renamed to version 1.7+ (no changes to this file)        *
*       18-Aug-93 : renamed to version 1.8 (no changes to this file)         *
*       17-Sep-93 : renamed to version 1.9 (no changes to this file)         *
*       29-Jun-95 : changes for version 1.10 :                               *
*                   - replaced loop in nextelement() to save reference past  *
*                     end of array (thanks to Kevin Maylsiak)                *
*       11-Jul-96 : changes for version 2.0 :                                *
*                   - added alloc_error()                                    *
*                   - added dynamic allocation                               *
*       21-Oct-98 : use 077777 in place of INFINITY for CLEANUP()            *
*        9-Jan-00 : added nautil_check()                                     *
*       12-Feb-00 : did a little formating of the code                       *
*       28-May-00 : added nautil_freedyn()                                   *
*       16-Aug-00 : added OLDNAUTY behaviour                                 *
*       16-Nov-00 : moved graph-specific things to naugraph.c                *
*                   use function prototypes, remove UPROC, nvector           *
*       22-Apr-01 : added code for compilation into Magma                    *
*                   removed nautil_null()                                    *
*                   removed EXTDEFS and included labelorg                    *
*       21-Nov-01 : use NAUTYREQUIRED in nautil_check()                      *
*       26-Jun-02 : revised permset() to avoid fetch past the end of         *
*                     the array (thanks to Jan Kieffer)                      *
*       17-Nov-03 : changed INFINITY to NAUTY_INFINITY                       *
*       14-Sep-04 : extended prototypes to recursive functions               *
*       23-Nov-06 : replave targetcell() by maketargetcell()                 *
*       10-Dec-06 : remove BIGNAUTY                                          *
*       10-Dec-10 : remove shortish and permutation types                    *
*       11-May-10 : use sorttemplates.c                                      *
*       27-Mar-11 : add writegroupsize()                                     *
*       15-Jan-12 : add TLS_ATTR attributes                                  *
*       16-Sep-12 : small change to objoin(), more efficient for sparse case *
*       22-Sep-12 : change documentation of orbjoin()                        *
*       18-Jan-12 : changes for version 2.6 :                                *
*                 - declare nauty_kill_request                               *
*                                                                            *
*****************************************************************************/

#define ONE_WORD_SETS
#include "nauty.h"
#ifdef NAUTY_IN_MAGMA
#include "io.e"
#endif

    /* macros for hash-codes: */
    /* Don't use NAUTY_INFINITY here as that would make the canonical
     * labelling depend on whether BIGNAUTY is in operation */
#define MASH(l,i) ((((l) ^ 065435) + (i)) & 077777)
    /* : expression whose long value depends only on long l and int/long i.
         Anything goes, preferably non-commutative. */

#define CLEANUP(l) ((int)((l) % 077777))
    /* : expression whose value depends on long l and is less than 077777
         when converted to int then short.  Anything goes. */

#if  MAXM==1
#define M 1
#else
#define M m
#endif

#if !MAXN
DYNALLSTAT(int,workperm,workperm_sz);
#else
static TLS_ATTR int workperm[MAXN];
#endif

int labelorg = 0;   /* no TLS_ATTR on purpose */
volatile int nauty_kill_request = 0;   /* no TLS_ATTR on purpose */

/* aproto: header new_nauty_protos.h */

/*****************************************************************************
*                                                                            *
*  nextelement(set1,m,pos) = the position of the first element in set set1   *
*  which occupies a position greater than pos.  If no such element exists,   *
*  the value is -1.  pos can have any value less than n, including negative  *
*  values.                                                                   *
*                                                                            *
*  GLOBALS ACCESSED: none                                                    *
*                                                                            *
*****************************************************************************/

int
nextelement(set *set1, int m, int pos)
{
    setword setwd;

#if  MAXM==1
    if (pos < 0) setwd = set1[0];
    else         setwd = set1[0] & BITMASK(pos);

    if (setwd == 0) return -1;
    else            return FIRSTBITNZ(setwd);
#else
    int w;

    if (pos < 0)
    {
        w = 0;
        setwd = set1[0];
    }
    else
    {
        w = SETWD(pos);
        setwd = set1[w] & BITMASK(SETBT(pos));
    }

    for (;;)
    {
        if (setwd != 0) return  TIMESWORDSIZE(w) + FIRSTBITNZ(setwd);
        if (++w == m) return -1;
        setwd = set1[w];
    }

#endif
}

/*****************************************************************************
*                                                                            *
*  permset(set1,set2,m,perm)  defines set2 to be the set                     *
*  {perm[i] | i in set1}.                                                    *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>,leftbit<r>                                       *
*                                                                            *
*****************************************************************************/

void
permset(set *set1, set *set2, int m, int *perm)
{
    setword setw;
    int pos,w,b;

    EMPTYSET(set2,m);

#if  MAXM==1
    setw = set1[0];
    while (setw  != 0)
    {
        TAKEBIT(b,setw);
        pos = perm[b];
        ADDELEMENT(set2,pos);
    }
#else
    for (w = 0; w < m; ++w)
    {
        setw = set1[w];
        while (setw != 0)
        {
            TAKEBIT(b,setw);
            pos = perm[TIMESWORDSIZE(w)+b];
            ADDELEMENT(set2,pos);
        }
    }
#endif
}

/*****************************************************************************
*                                                                            *
*  putstring(f,s) writes the nul-terminated string s to file f.              *
*                                                                            *
*****************************************************************************/

void
putstring(FILE *f, char *s)
{
    while (*s != '\0')
    {
        PUTC(*s,f);
        ++s;
    }
}

/*****************************************************************************
*                                                                            *
*  itos(i,s) converts the int i to a nul-terminated decimal character        *
*  string s.  The value returned is the number of characters excluding       *
*  the nul.                                                                  *
*                                                                            *
*  GLOBALS ACCESSED: NONE                                                    *
*                                                                            *
*****************************************************************************/

int
itos(int i, char *s)
{
    int digit,j,k;
    char c;
    int ans;

    if (i < 0)
    {
        k = 0;
        i = -i;
        j = 1;
        s[0] = '-';
    }
    else
    {
        k = -1;
        j = 0;
    }

    do
    {
        digit = i % 10;
        i = i / 10;
        s[++k] = (char)(digit + '0');
    }
    while (i);

    s[k+1] = '\0';
    ans = k + 1;

    for (; j < k; ++j, --k)
    {
        c = s[j];
        s[j] = s[k];
        s[k] = c;
    }

    return ans;
}

/*****************************************************************************
*                                                                            *
*  orbits represents a partition of {0,1,...,n-1}, by orbits[i] = the        *
*  smallest element in the same cell as i.  map[] is any array with values   *
*  in {0,1,...,n-1}.  orbjoin(orbits,map,n) joins the cells of orbits[]      *
*  together to the minimum extent such that for each i, i and map[i] are in  *
*  the same cell.  The function value returned is the new number of cells.   *
*                                                                            *
*  GLOBALS ACCESSED: NONE                                                    *
*                                                                            *
*****************************************************************************/

int
orbjoin(int *orbits, int *map, int n)
{
    int i,j1,j2;

    for (i = 0; i < n; ++i)
    if (map[i] != i)
    {
        j1 = orbits[i];
        while (orbits[j1] != j1) j1 = orbits[j1];
        j2 = orbits[map[i]];
        while (orbits[j2] != j2) j2 = orbits[j2];

        if (j1 < j2)      orbits[j2] = j1;
        else if (j1 > j2) orbits[j1] = j2;
    }

    j1 = 0;
    for (i = 0; i < n; ++i)
        if ((orbits[i] = orbits[orbits[i]]) == i) ++j1;

    return j1;
}

/*****************************************************************************
*                                                                            *
*  writeperm(f,perm,cartesian,linelength,n) writes the permutation perm to   *
*  the file f.  The cartesian representation (i.e. perm itself) is used if   *
*  cartesian != FALSE; otherwise the cyclic representation is used.  No      *
*  more than linelength characters (not counting '\n') are written on each   *
*  line, unless linelength is ridiculously small.  linelength<=0 causes no   *
*  line breaks at all to be made.  The global int labelorg is added to each  *
*  vertex number.                                                            *
*                                                                            *
*  GLOBALS ACCESSED: itos(),putstring()                                      *
*                                                                            *
*****************************************************************************/

void
writeperm(FILE *f, int *perm, boolean cartesian, int linelength, int n)
{
    int i,k,l,curlen,intlen;
    char s[30];

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n,"writeperm");
#endif

    /* CONDNL(x) writes end-of-line and 3 spaces if x characters
       won't fit on the current line. */
#define CONDNL(x) if (linelength>0 && curlen+(x)>linelength)\
              {putstring(f,"\n   ");curlen=3;}

    curlen = 0;
    if (cartesian)
    {
        for (i = 0; i < n; ++i)
        {
            intlen = itos(perm[i]+labelorg,s);
            CONDNL(intlen+1);
            PUTC(' ',f);
            putstring(f,s);
            curlen += intlen + 1;
        }
        PUTC('\n',f);
    }
    else
    {
        for (i = n; --i >= 0;) workperm[i] = 0;

        for (i = 0; i < n; ++i)
        {
            if (workperm[i] == 0 && perm[i] != i)
            {
                l = i;
                intlen = itos(l+labelorg,s);
                if (curlen > 3) CONDNL(2*intlen+4);
                PUTC('(',f);
                do
                {
                    putstring(f,s);
                    curlen += intlen + 1;
                    k = l;
                    l = perm[l];
                    workperm[k] = 1;
                    if (l != i)
                    {
                        intlen = itos(l+labelorg,s);
                        CONDNL(intlen+2);
                        PUTC(' ',f);
                    }
                }
                while (l != i);
                PUTC(')',f);
                ++curlen;
            }
        }

        if (curlen == 0) putstring(f,"(1)\n");
        else             PUTC('\n',f);
    }
}

/*****************************************************************************
*                                                                            *
*  fmperm(perm,fix,mcr,m,n) uses perm to construct fix and mcr.  fix         *
*  contains those points are fixed by perm, while mcr contains the set of    *
*  those points which are least in their orbits.                             *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>                                                  *
*                                                                            *
*****************************************************************************/

void
fmperm(int *perm, set *fix, set *mcr, int m, int n)
{
    int i,k,l;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n,"writeperm");
#endif

    EMPTYSET(fix,m);
    EMPTYSET(mcr,m);

    for (i = n; --i >= 0;) workperm[i] = 0;

    for (i = 0; i < n; ++i)
        if (perm[i] == i)
        {
            ADDELEMENT(fix,i);
            ADDELEMENT(mcr,i);
        }
        else if (workperm[i] == 0)
        {
            l = i;
            do
            {
                k = l;
                l = perm[l];
                workperm[k] = 1;
            }
            while (l != i);

            ADDELEMENT(mcr,i);
        }
}

/*****************************************************************************
*                                                                            *
*  fmptn(lab,ptn,level,fix,mcr,m,n) uses the partition at the specified      *
*  level in the partition nest (lab,ptn) to make sets fix and mcr.  fix      *
*  represents the points in trivial cells of the partition, while mcr        *
*  represents those points which are least in their cells.                   *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>                                                  *
*                                                                            *
*****************************************************************************/

void
fmptn(int *lab, int *ptn, int level, set *fix, set *mcr, int m, int n)
{
    int i,lmin;

    EMPTYSET(fix,m);
    EMPTYSET(mcr,m);

    for (i = 0; i < n; ++i)
        if (ptn[i] <= level)
        {
            ADDELEMENT(fix,lab[i]);
            ADDELEMENT(mcr,lab[i]);
        }
        else
        {
            lmin = lab[i];
            do
                if (lab[++i] < lmin) lmin = lab[i];
            while (ptn[i] > level);
            ADDELEMENT(mcr,lmin);
        }
}

#define SORT_TYPE1 int
#define SORT_TYPE2 int
#define SORT_OF_SORT 2
#define SORT_NAME sortparallel
#include "sorttemplates.c"

/*****************************************************************************
*                                                                            *
*  doref(g,lab,ptn,level,numcells,qinvar,invar,active,code,refproc,          *
*        invarproc,mininvarlev,maxinvarlev,invararg,digraph,m,n)             *
*  is used to perform a refinement on the partition at the given level in    *
*  (lab,ptn).  The number of cells is *numcells both for input and output.   *
*  The input active is the active set for input to the refinement procedure  *
*  (*refproc)(), which must have the argument list of refine().              *
*  active may be arbitrarily changed.  invar is used for working storage.    *
*  First, (*refproc)() is called.  Then, if invarproc!=NULL and              *
*  |mininvarlev| <= level <= |maxinvarlev|, the routine (*invarproc)() is    *
*  used to compute a vertex-invariant which may refine the partition         *
*  further.  If it does, (*refproc)() is called again, using an active set   *
*  containing all but the first fragment of each old cell.  Unless g is a    *
*  digraph, this guarantees that the final partition is equitable.  The      *
*  arguments invararg and digraph are passed to (*invarproc)()               *
*  uninterpretted.  The output argument code is a composite of the codes     *
*  from all the calls to (*refproc)().  The output argument qinvar is set    *
*  to 0 if (*invarproc)() is not applied, 1 if it is applied but fails to    *
*  refine the partition, and 2 if it succeeds.                               *
*  See the file nautinv.c for a further discussion of vertex-invariants.     *
*  Note that the dreadnaut I command generates a call to  this procedure     *
*  with level = mininvarlevel = maxinvarlevel = 0.                           *
*                                                                            *
*****************************************************************************/

void
doref(graph *g, int *lab, int *ptn, int level, int *numcells,
     int *qinvar, int *invar, set *active, int *code,
     void (*refproc)(graph*,int*,int*,int,int*,int*,set*,int*,int,int),
     void (*invarproc)(graph*,int*,int*,int,int,int,int*,
                                                 int,boolean,int,int),
     int mininvarlev, int maxinvarlev, int invararg,
     boolean digraph, int m, int n)
{
    int pw;
    int i,cell1,cell2,nc,tvpos,minlev,maxlev;
    long longcode;
    boolean same;

#if !MAXN 
    DYNALLOC1(int,workperm,workperm_sz,n,"doref"); 
#endif

    if ((tvpos = nextelement(active,M,-1)) < 0) tvpos = 0;

    (*refproc)(g,lab,ptn,level,numcells,invar,active,code,M,n);

    minlev = (mininvarlev < 0 ? -mininvarlev : mininvarlev);
    maxlev = (maxinvarlev < 0 ? -maxinvarlev : maxinvarlev);
    if (invarproc != NULL && *numcells < n
                        && level >= minlev && level <= maxlev)
    {
        (*invarproc)(g,lab,ptn,level,*numcells,tvpos,invar,invararg,
                                                             digraph,M,n);
        EMPTYSET(active,m);
        for (i = n; --i >= 0;) workperm[i] = invar[lab[i]];
        nc = *numcells;
        for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
        {
            pw = workperm[cell1];
            same = TRUE;
            for (cell2 = cell1; ptn[cell2] > level; ++cell2)
                if (workperm[cell2+1] != pw) same = FALSE;

            if (same) continue;

            sortparallel(workperm+cell1, lab+cell1, cell2-cell1+1);

            for (i = cell1 + 1; i <= cell2; ++i)
                if (workperm[i] != workperm[i-1])
                {
                    ptn[i-1] = level;
                    ++*numcells;
                    ADDELEMENT(active,i);
                }
        }

        if (*numcells > nc)
        {
            *qinvar = 2;
            longcode = *code;
            (*refproc)(g,lab,ptn,level,numcells,invar,active,code,M,n);
            longcode = MASH(longcode,*code);
            *code = CLEANUP(longcode);
        }
        else
            *qinvar = 1;
    }
    else
        *qinvar = 0;
}

/*****************************************************************************
*                                                                            *
*  maketargetcell(g,lab,ptn,level,tcell,tcellsize,&cellpos,                  *
*                 tc_level,digraph,hint,targetcell,m,n)                      *
*  calls targetcell() to determine the target cell at the specified level    *
*  in the partition nest (lab,ptn).  It must be a nontrivial cell (if not,   *
*  the first cell.  The intention of hint is that, if hint >= 0 and there    *
*  is a suitable non-trivial cell starting at position hint in lab,          *
*  that cell is chosen.                                                      *
*  tc_level and digraph are input options.                                   *
*  When a cell is chosen, tcell is set to its contents, *tcellsize to its    *
*  size, and cellpos to its starting position in lab.                        *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>                                                  *
*                                                                            *
*****************************************************************************/

void
maketargetcell(graph *g, int *lab, int *ptn, int level, set *tcell,
       int *tcellsize, int *cellpos, int tc_level, boolean digraph,
       int hint,
       int (*targetcell)(graph*,int*,int*,int,int,boolean,int,int,int),
       int m, int n)
{
    int i,j,k;

    i = (*targetcell)(g,lab,ptn,level,tc_level,digraph,hint,m,n);
    for (j = i + 1; ptn[j] > level; ++j) {}

    *tcellsize = j - i + 1;

    EMPTYSET(tcell,m);
    for (k = i; k <= j; ++k) ADDELEMENT(tcell,lab[k]);

    *cellpos = i;
}

/*****************************************************************************
*                                                                            *
*  shortprune(set1,set2,m) ANDs the contents of set set2 into set set1.      *
*                                                                            *
*  GLOBALS ACCESSED: NONE                                                    *
*                                                                            *
*****************************************************************************/

void
shortprune(set *set1, set *set2, int m)
{
    int i;

    for (i = 0; i < M; ++i) INTERSECT(set1[i],set2[i]);
}

/*****************************************************************************
*                                                                            *
*  breakout(lab,ptn,level,tc,tv,active,m) operates on the partition at       *
*  the specified level in the partition nest (lab,ptn).  It finds the        *
*  element tv, which is in the cell C starting at index tc in lab (it had    *
*  better be) and splits C in the two cells {tv} and C\{tv}, in that order.  *
*  It also sets the set active to contain just the element tc.               *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>                                                  *
*                                                                            *
*****************************************************************************/

void
breakout(int *lab, int *ptn, int level, int tc, int tv,
     set *active, int m)
{
    int i,prev,next;

    EMPTYSET(active,m);
    ADDELEMENT(active,tc);

    i = tc;
    prev = tv;

    do
    {
        next = lab[i];
        lab[i++] = prev;
        prev = next;
    }
    while (prev != tv);

    ptn[tc] = level;
}

/*****************************************************************************
*                                                                            *
*  longprune(tcell,fix,bottom,top,m) removes zero or elements of the set     *
*  tcell.  It is assumed that addresses bottom through top-1 contain         *
*  contiguous pairs of sets (f1,m1),(f2,m2), ... .  tcell is intersected     *
*  with each mi such that fi is a subset of fix.                             *
*                                                                            *
*  GLOBALS ACCESSED: NONE                                                    *
*                                                                            *
*****************************************************************************/

void
longprune(set *tcell, set *fix, set *bottom, set *top, int m)
{
    int i;

    while (bottom < top)
    {
        for (i = 0; i < M; ++i)
            if (NOTSUBSET(fix[i],bottom[i])) break;
        bottom += M;

        if (i == M)
            for (i = 0; i < M; ++i) INTERSECT(tcell[i],bottom[i]);
        bottom += M;
    }
}

/*****************************************************************************
*                                                                            *
*  writegroupsize(f,gpsize1,gpsize2) writes a real number gpsize1*10^gpsize2 *
*  It is assumed that gpsize1 >= 1 and that gpsize1 equals an integer in the *
*  case that gpsize2==0.  These assumptions are true for group sizes         *
*  computed by nauty.                                                        *
*                                                                            *
*****************************************************************************/

void
writegroupsize(FILE *f, double gpsize1, int gpsize2)
{
    if (gpsize2 == 0)
        fprintf(f,"%.0f",gpsize1+0.1);
    else
    {   
        while (gpsize1 >= 10.0)
        {   
            gpsize1 /= 10.0;
            ++gpsize2;
        }
        fprintf(f,"%14.12fe%d",gpsize1,gpsize2);
    }
}

/*****************************************************************************
*                                                                            *
*  nautil_check() checks that this file is compiled compatibly with the      *
*  given parameters.   If not, call exit(1).                                 *
*                                                                            *
*****************************************************************************/

void
nautil_check(int wordsize, int m, int n, int version)
{
    if (wordsize != WORDSIZE)
    {
        fprintf(ERRFILE,"Error: WORDSIZE mismatch in nautil.c\n");
        exit(1);
    }

#if MAXN
    if (m > MAXM)
    {
        fprintf(ERRFILE,"Error: MAXM inadequate in nautil.c\n");
        exit(1);
    }

    if (n > MAXN)
    {
        fprintf(ERRFILE,"Error: MAXN inadequate in nautil.c\n");
        exit(1);
    }
#endif

    if (version < NAUTYREQUIRED)
    {
        fprintf(ERRFILE,"Error: nautil.c version mismatch\n");
        exit(1);
    }
}

/*****************************************************************************
*                                                                            *
*  alloc_error() writes a message and exits.  Used by DYNALLOC? macros.      *
*                                                                            *
*****************************************************************************/

void
alloc_error(const char *s)
{
    fprintf(ERRFILE,"Dynamic allocation failed: %s\n",s);
    exit(2);
}

/*****************************************************************************
*                                                                            *
*  nautil_freedyn() - free the dynamic memory in this module                 *
*                                                                            *
*****************************************************************************/

void
nautil_freedyn(void)
{
#if !MAXN
    DYNFREE(workperm,workperm_sz);
#endif
}
