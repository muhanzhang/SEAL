/*****************************************************************************
*                                                                            *
* miscellaneous utilities for use with nauty 2.6.                            *
* None of these procedures are needed by nauty, but all are by dreadnaut.    *
*                                                                            *
*   Copyright (1984-2016) Brendan McKay.  All rights reserved.               *
*   Subject to waivers and disclaimers in nauty.h.                           *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       10-Nov-87 : final changes for version 1.2                            *
*        5-Dec-87 : changes made for version 1.3 :                           *
*                   - added procedures readinteger() and readstring()        *
*                   - replaced all uses of fscanf() by appropriate uses      *
*                     of readinteger() or readstring()                       *
*                   - "N:" is now illegal in readgraph() if N is too large   *
*                     or too small                                           *
*       28-Sep-88 : renamed to version 1.4 (no changes to this file)         *
*       23-Mar-89 : changes for version 1.5 :                                *
*                   - declared key in hash()                                 *
*                   - changed file name to naututil.c                        *
*       29-Mar-89 : - declared workperm[] and workset[], and modified        *
*                     many routines to use them.                             *
*                   - added putmapping()                                     *
*                   - reworked some code in mathon() and rangraph()          *
*        3-Apr-89 : - changed putorbits() to use O(n) algorithm              *
*        5-Apr-89 : - modifed readgraph() to not require fresh line          *
*                   - changed MAKEEMPTY uses to EMPTYSET uses                *
*       26-Apr-89 : - moved ptncode() and equitable() to nautaux.c           *
*                   - added putquotient()                                    *
*       18-Aug-89 : - modified putset() to use "i:j" syntax optionally       *
*                   - modified putorbits() to use putset()                   *
*                   - modified calling sequence for mathon()                 *
*       19-Aug-90 : - changed delimeter arg of copycomment to int            *
*       14-Oct-90 : renamed to version 1.6 (no changes to this file)         *
*       23-Jan-91 : changes for version 1.7 :                                *
*                   - fixed bug in complement()                              *
*       27-Aug-92 : - made linelength <= 0 mean no line breaks               *
*        5-Jun-93 : renamed to version 1.7+ (no changes to this file)        *
*       18-Aug-93 : renamed to version 1.8 (no changes to this file)         *
*       17-Sep-93 : renamed to version 1.9 (no changes to this file)         *
*       13-Jul-96 : changes for version 2.0 :                                *
*                   - added dynamic allocation                               *
*                   - added limit parameter to readstring                    *
*                   - added readvperm() and sublabel()                       *
*       31-Aug-96 : - converted from RAN to KRAN                             *
*        6-Feb-97 : - corrected DYNALLOC1 call in putmapping                 *
*       10-Dec-97 : - KRAN now initialises automatically                     *
*        9-Jan-00 : - added naututil_check()                                 *
*       12-Feb-00 : - some minor code formatting                             *
*       16-Nov-00 : - changes as listed in nauty.h                           *
*       23-Apr-01 : changes for version 2.1 :                                *
*                   - removed EXTDEFS                                        *
*        2-Jun-01 : - added converse()                                       *
*       21-Nov-01 : use NAUTYREQUIRED in naututil_check()                    *
*       11-Apr-03 : changes for version 2.2 :                                *
*                   - added rangraph2()                                      *
*       17-Nov-03 : changed INFINITY to NAUTY_INFINITY                       *
*       10-Dec-06 : removed BIGNAUTY                                         *
*        4-Nov-09 : added readgraph_sg, putgraph_sg, putcanon_sg             *
*       10-Nov-09 : removed shortish and permutation types                   *
*       14-Nov-09 : added relabel_sg(), putdegs_sg(), sublabel_sg()          *
*       19-Nov-09 : added individualise()                                    *
*       20-Nov-09 : added hashgraph_sg(), listhash(), hashgraph()            *
*       19-Dec-09 : added ranreg_sg(), rangraph2_sg()                        *
*       21-May-10 : conform to type changes in sparsegraph                   *
*        5-Jun-10 : add mathon_sg()                                          *
*       10-Jun-10 : add putquotient_sg() and complement_sg()                 *
*       26-Jan-11 : fix nde error in sublabel_sg()                           *
*       15-Jan-12 : add TLS_ATTR attributes                                  *
*        3-Mar-12 : add putorbitsplus() and putset_firstbold()               *
*                 : write orbit sizes if not trivial                         *
*       17-Mar-12 : move seed definition to naututil.h                       *
*       20-Sep-12 : allow quoted strings in readstring()                     *
*       20-Sep-12 : the first argument of ungetc is int, not char            *
*        4-Mar-13 : remove a side-effect issue in setinter()                 *
*       17-Dec-15 : add readgraph_swg() and update putgraph_sg()             *
*       22-Jan-16 : add readintger_sl() and getint_sl()                      *
*       29-Feb-16 : add subpartition()                                       *
*        6-Apr-16 : add countcells(), make subpartition return a count       *
*                                                                            *
*****************************************************************************/

#define ONE_WORD_SETS
#include "naututil.h"    /* which includes nauty.h, nautinv.h and stdio.h */
#include "nausparse.h"

#if  MAXM==1
#define M 1
#else
#define M m
#endif

#if !MAXN
DYNALLSTAT(int,workperm,workperm_sz);
DYNALLSTAT(set,workset,workset_sz);
#else
static TLS_ATTR int workperm[MAXN+2];   /* used for scratch purposes */
static TLS_ATTR set workset[MAXM];      /* used for scratch purposes */
#endif

#define ECHUNKSIZE 1000  /* should be a multiple of 2 */
typedef struct echunk {struct echunk *next; int edge[ECHUNKSIZE];} echunk;
static TLS_ATTR echunk first_echunk = {NULL,{0}};
typedef struct echunkw {struct echunkw *next; \
 struct {int v1,v2; sg_weight wt;} edge[ECHUNKSIZE];} echunkw;
static TLS_ATTR echunkw first_echunkw = {NULL,{0,0,0}};

#ifdef  NLMAP
#define GETNW(c,f) do c = getc(f); while (c==' '||c=='\t')
#define GETNWC(c,f) do c = getc(f); while (c==' '||c==','||c=='\t')
#define GETNWL(c,f) do c = getc(f); while (c==' '||c=='\n'||c=='\t')
#else
#define GETNW(c,f) do c = getc(f); while (c==' '||c=='\t'||c=='\r')
#define GETNWC(c,f) do c = getc(f); while (c==' '||c==','||c=='\t'||c=='\r')
#define GETNWL(c,f) do c = getc(f); while (c==' '||c=='\n'||c=='\t'||c=='\r')
#endif

#define ISDIGIT(c) ((c) >= '0' && (c) <= '9')

static const long fuzz1[] = {1984625421L, 971524688L,1175081625L, 377165387L};
static const long fuzz2[] = {2001381726L,1615243355L, 191176436L,1212176501L};
#define FUZZ1(x) ((x) ^ fuzz1[(x)&3L])
#define FUZZ2(x) ((x) ^ fuzz2[(x)&3L])

#define SORT_OF_SORT 1
#define SORT_NAME sort1int
#define SORT_TYPE1 int
#include "sorttemplates.c"   /* define sort1int(a,n) */

/*****************************************************************************
*                                                                            *
*  setinter(set1,set2,m) = the number of elements in the intersection of     *
*  the sets set1 and set2.                                                   *
*                                                                            *
*****************************************************************************/

int
setinter(set *set1, set *set2, int m)
{
    setword x;

#if  MAXM==1
    if ((x = *set1 & *set2) != 0) return POPCOUNT(x);
    else                          return 0;
#else
    int count,i;

    count = 0;
    for (i = m; --i >= 0;)
    {
        if ((x = (*set1 & *set2)) != 0) count += POPCOUNT(x);
        ++set1;
	++set2;
    }

    return count;
#endif
}

/*****************************************************************************
*                                                                            *
*  setsize(set1,m) = the number of elements in the set set1.                 *
*                                                                            *
*****************************************************************************/

int
setsize(set *set1, int m)
{

#if  MAXM==1
    if (set1 != 0) return POPCOUNT(*set1);
    else           return 0;
#else
    int count,i;
    setword x;

    count = 0;
    for (i = m; --i >= 0;)
        if ((x = *set1++) != 0) count += POPCOUNT(x);

    return count;
#endif
}

/*****************************************************************************
*                                                                            *
*  flushline(f) reads from f until '\n' or EOF.                              *
*  If non-trivial characters are skipped, the user is informed.              *
*                                                                            *
*****************************************************************************/

void
flushline(FILE *f)
{
    boolean msg;
    int c;

    msg = FALSE;

    while ((c = getc(f)) != EOF && c != '\n')
        if (msg)
            PUTC((char)c,ERRFILE);
        else if (c != ' ' && c != '\t' && c != '\f' &&
                                  c != '\r' && c != ',')
        {
            msg = TRUE;
            fprintf(ERRFILE,"input skipped : '%c",(char)c);
        }
    if (msg) fprintf(ERRFILE,"'\n\n");
}

/*****************************************************************************
*                                                                            *
*  readinteger(f,&i) reads an optionally-signed integer from f, preceded by  *
*  any amount of white space.  The function value returned is TRUE if an     *
*  integer was found, FALSE otherwise.                                       *
*                                                                            *
*****************************************************************************/

boolean
readinteger(FILE *f, int *p)
{
    int c,ans,minus;

    GETNWL(c,f);
    if (!ISDIGIT(c) && c != '-' && c != '+')
    {
        if (c != EOF) ungetc(c,f);
        return FALSE;
    }

    minus = c == '-';
    ans = (c == '-' || c == '+' ? 0 : c - '0');

    c = getc(f);
    while (ISDIGIT(c))
    {
        ans *= 10;
        ans += c - '0';
        c = getc(f);
    }

    if (c != EOF) ungetc(c,f);

    *p = (minus ? -ans : ans);
    return TRUE;
}

/*****************************************************************************
*                                                                            *
*  readinteger_sl(f,&i) reads an optionally-signed integer from f, preceded  *
*  by any amount of white space except newlines.  The function value         *
*  returned is TRUE if an integer was found, FALSE otherwise.                *
*                                                                            *
*****************************************************************************/

boolean
readinteger_sl(FILE *f, int *p)
{
    int c,ans,minus;

    GETNW(c,f);
    if (!ISDIGIT(c) && c != '-' && c != '+')
    {
        if (c != EOF) ungetc(c,f);
        return FALSE;
    }

    minus = c == '-';
    ans = (c == '-' || c == '+' ? 0 : c - '0');

    c = getc(f);
    while (ISDIGIT(c))
    {
        ans *= 10;
        ans += c - '0';
        c = getc(f);
    }

    if (c != EOF) ungetc(c,f);

    *p = (minus ? -ans : ans);
    return TRUE;
}

/*****************************************************************************
*                                                                            *
*  readstring(f,s,slen) reads a string from f. First any amount of white     *
*  space is skipped (including newlines).  If the next character is a        *
*  double-quote, everything after that before the next double-quote or       *
*  newline is put into s.  If the next character is not a double-quote,      *
*  everything before the next white space is put into s.  A nul is added,    *
*  but no more than slen characters are ever put into s.  The function       *
*  value is TRUE if a string was found and FALSE otherwise.                  *
*                                                                            *
*****************************************************************************/

boolean
readstring(FILE *f, char *s, int slen)
{
    int c;
    char *slim;

    slim = s + slen - 1;
    GETNWL(c,f);
    if (c == EOF)
    {
        *s = '\0';
        return FALSE;
    }

    if (c == '"')
    {
	c = getc(f);
	while (c != '"' && c != '\n' && c != '\r' && c != EOF)
	{
	    if (s <= slim) *s++ = (char)c;
            c = getc(f);
        }
	if (c != '"' && c != EOF) ungetc(c,f);
    }
    else
    {
        if (s <= slim) *s++ = (char)c;
        c = getc(f);
        while (c != ' ' && c != '\n' && c != '\t' && c != '\r' && c != EOF)
        {
            if (s <= slim) *s++ = (char)c;
            c = getc(f);
        }
        if (c != EOF) ungetc(c,f);
    }
    if (s <= slim) *s = '\0';
    else           *slim = '\0';

    return TRUE;
}

/*****************************************************************************
*                                                                            *
*  getint(f) reads an integer from f, optionally preceded by '='             *
*  and white space.  -1 is returned if the attempt was unsuccessful.         *
*                                                                            *
*****************************************************************************/

int
getint(FILE *f)
{
    int i,c;

    GETNWL(c,f);
    if (c != '=') ungetc(c,f);

    if (readinteger(f,&i)) return i;
    else                   return -1;
}

/*****************************************************************************
*                                                                            *
* getint_sl(f) reads an integer from f, optionally preceded by '=' and       *
* white space other than newlines.  -1 is returned if the attempt was        *
* unsuccessful.                                                              *
*****************************************************************************/

int
getint_sl(FILE *f)
{
    int i,c;

    GETNW(c,f);
    if (c != '=') ungetc(c,f);

    if (readinteger_sl(f,&i)) return i;
    else                      return -1;
}

/*****************************************************************************
*                                                                            *
*  putset(f,set1,curlenp,linelength,m,compress)   writes the set set1 to     *
*  file f using at most linelength characters per line (excluding '\n').     *
*  A value of linelength <= 0 dictates no line breaks at all.                *
*  *curlenp is the number of characters on the line so far; it is updated.   *
*  A range j1,j1+1,...,j2 for j2-j1>=2 is written as "j1:j2" if compress     *
*  is nonzero (eg. TRUE); otherwise each element is written separately.      *
*  No final '\n' is written.  labelorg is used.                              *
*                                                                            *
*  FUNCTIONS CALLED: nextelement(),itos()                                    *
*                                                                            *
*****************************************************************************/

void
putset(FILE *f, set *set1, int *curlenp, int linelength,
       int m, boolean compress)
{
    int slen,j1,j2;
    char s[40];

    j1 = -1;
    while ((j1 = nextelement(set1,m,j1)) >= 0)
    {
        j2 = j1;
        if (compress)
        {
            while (nextelement(set1,m,j2) == j2 + 1) ++j2;
            if (j2 == j1+1) j2 = j1;
        }
        slen = itos(j1+labelorg,s);
        if (j2 >= j1 + 2)
        {
            s[slen] = ':';
            slen += 1 + itos(j2+labelorg,&s[slen+1]);
        }

        if (linelength > 0 && *curlenp + slen + 1 >= linelength)
        {
            fprintf(f,"\n   ");
            *curlenp = 3;
        }
        fprintf(f," %s",s);
        *curlenp += slen + 1;
        j1 = j2;
    }
}

/*****************************************************************************
*                                                                            *
*  putset_firstbold(f,set1,curlenp,linelength,m,compress) is the same as     *
*  putset(f,set1,curlenp,linelength,m,compress) except that the first        *
*  element of the set is written bold.  This is only useful when output is   *
*  to a device that interprets ANSI control sequences.                       *
*                                                                            *
*  FUNCTIONS CALLED: nextelement(),itos()                                    *
*                                                                            *
*****************************************************************************/

void
putset_firstbold(FILE *f, set *set1, int *curlenp, int linelength,
       int m, boolean compress)
{
    int slen,slen1,j1,j2;
    char s[50],c;
    boolean bold;

    bold = TRUE;
    j1 = -1;
    while ((j1 = nextelement(set1,m,j1)) >= 0)
    {
        j2 = j1;
        if (compress)
        {
            while (nextelement(set1,m,j2) == j2 + 1) ++j2;
            if (j2 == j1+1) j2 = j1;
        }
        slen1 = slen = itos(j1+labelorg,s);
        if (j2 >= j1 + 2)
        {
            s[slen] = ':';
            slen += 1 + itos(j2+labelorg,&s[slen+1]);
        }
	c = s[slen1];

        if (linelength > 0 && *curlenp + slen + 1 >= linelength)
        {
            fprintf(f,"\n   ");
            *curlenp = 3;
        }
        if (bold)
        {
	    s[slen1] = '\0';
            fprintf(f," \033[1m%s\033[0m",s);
            s[slen1] = c;
            fprintf(f,"%s",&s[slen1]);
            bold = FALSE;
        }
        else
            fprintf(f," %s",s);

        *curlenp += slen + 1;
        j1 = j2;
    }
}

/*****************************************************************************
*                                                                            *
*  readgraph(f,g,digraph,prompt,edit,linelength,m,n) reads a graph g from f. *
*  Commands: (There is always a "current vertex" v, initially labelorg;      *
*             n is an unsigned integer.)                                     *
*  n  : add edge (v,n)                                                       *
*  -n : delete edge (v,n)                                                    *
*  n: : set v := n, and exit if v >= n.                                      *
*  ?  : type neighbours of vertex v                                          *
*  ;  : increment v, and exit if v >= n.                                     *
*  .  : exit                                                                 *
*  !  : skip rest of input line                                              *
*                                                                            *
* If digraph==FALSE, loops are illegal and (x,y) => (y,x)                    *
* If edit==FALSE, the graph is initialized to empty.                         *
* If prompt==TRUE, prompts are written to PROMPTFILE.                        *
* linelength is a limit on the number of characters per line caused by '?'   *
* A value of linelength <= 0 dictates no line breaks at all.                *
* labelorg is used.                                                          *
*                                                                            *
* FUNCTIONS CALLED : putset()                                                *
*                                                                            *
*****************************************************************************/

void
readgraph(FILE *f, graph *g, boolean digraph, boolean prompt,
      boolean edit, int linelength, int m, int n)
{
    int v,c;
    int curlen,w;
    graph *gv;
    boolean neg;

    if (!edit)
        for (v = 0, gv = g; v < n; ++v, gv += M) EMPTYSET(gv,m);

    v = 0;
    gv = g;
    neg = FALSE;

    while (TRUE)
    {
        GETNWC(c,f);
        if (ISDIGIT(c))
        {
            ungetc(c,f);
            readinteger(f,&w);
            w -= labelorg;
            if (neg)
            {
                neg = FALSE;
                if (w < 0 || w >= n || (!digraph && w == v))
                    fprintf(ERRFILE,"illegal edge (%d,%d) ignored\n\n",
                            v+labelorg,w+labelorg);
                else
                {
                    DELELEMENT(gv,w);
                    if (!digraph) DELELEMENT(GRAPHROW(g,w,M),v);
                }
            }
            else
            {
                GETNWC(c,f);
                if (c == ':')
                    if (w < 0 || w >= n)
                        fprintf(ERRFILE,
                                "illegal vertex number %d ignored\n\n",
                                w+labelorg);
                    else
                    {
                        v = w;
                        gv = GRAPHROW(g,v,M);
                    }
                else
                {
                    ungetc(c,f);
                    if (w < 0 || w >= n || (!digraph && w == v))
                        fprintf(ERRFILE,"illegal edge (%d,%d) ignored\n\n",
                                v+labelorg,w+labelorg);
                    else
                    {
                        ADDELEMENT(gv,w);
                        if (!digraph) ADDELEMENT(GRAPHROW(g,w,M),v);
                    }
                }
            }
        }
        else switch(c)
        {
            case ';':
                neg = FALSE;
                ++v;
                if (v >= n) return;
                gv = GRAPHROW(g,v,M);
                break;
            case '?':
                neg = FALSE;
                fprintf(PROMPTFILE,"%2d : ",v+labelorg);
                curlen = 5;
                putset(PROMPTFILE,gv,&curlen,linelength,M,FALSE);
                fprintf(PROMPTFILE,";\n");
                break;
            case '\n':
                neg = FALSE;
                if (prompt) fprintf(PROMPTFILE,"%2d : ",v+labelorg);
                break;
            case EOF:
            case '.':
                return;
            case '-':
                neg = TRUE;
                break;
            case '!':
                do
                    c = getc(f);
                while (c != '\n' && c != EOF);
                if (c == '\n') ungetc(c,f);
                break;
            default :
                fprintf(ERRFILE,"illegal char '%c' - use '.' to exit\n\n",
                        (char)c);
        }
    }
}

/**************************************************************************/

void
ranreg_sg(sparsegraph *sg, int degree, int n)
/* Make a random regular simple undirected graph.
 * For MAXN!=0, the maximum degree is MAXREG.
 * sg must be initialised
 */
{
    long i,k,v,w;
    boolean ok;
    int *dd,*ee;
    size_t *vv,nde,j;
 
#if MAXN
    int p[MAXREG*MAXN];
#else
    DYNALLSTAT(int,p,p_sz);

    DYNALLOC2(int,p,p_sz,degree,n,"genrang");
#endif

    nde = (size_t)n * (size_t)degree;

    SG_ALLOC(*sg,n,nde,"ranreg_sg");
    SG_VDE(sg,vv,dd,ee);
    DYNFREE(sg->w,sg->wlen);

    sg->nv = n;
    sg->nde = nde;

    for (i = j = 0; i < n; ++i)
        for (k = 0; k < degree; ++k) p[j++] = i;

    for (i = 0; i < n; ++i) vv[i] = i*(size_t)degree;

    do
    {   
        ok = TRUE;

        for (j = nde; j > 0; j -= 2)
        {   
            i = KRAN(j-1);
	    k = p[i];
	    if (k == p[j-1]) break;
	    p[i] = p[j-2]; p[j-2] = k;
        }
	if (j > 0) { ok = FALSE; continue; }

        for (i = 0; i < n; ++i) dd[i] = 0;

        for (j = nde; j > 0; )
        {   
            v = p[--j];
            w = p[--j];
            if (v != w)
	    {
                for (i = dd[w]; --i >= 0;) if (ee[vv[w]+i] == v) break;
                if (i >= 0) { ok = FALSE; break; }
            }
            ee[vv[w]+(dd[w]++)] = v;
            ee[vv[v]+(dd[v]++)] = w;
        }
    }
    while (!ok);
}

/*****************************************************************************
*                                                                            *
*  readgraph_sg(f,sg,digraph,prompt,linelength,n) reads a graph g from f.    *
*  Commands: (There is always a "current vertex" v, initially labelorg;      *
*             n is an unsigned integer.)                                     *
*  n  : add edge (v,n)                                                       *
*  -n : delete edge (v,n)                                                    *
*  n: : set v := n, and exit if v >= n.                                      *
*  ?  : type neighbours of vertex v     ** NOT IMPLEMENTED **                *
*  ;  : increment v, and exit if v >= n.                                     *
*  .  : exit                                                                 *
*  !  : skip rest of input line                                              *
*  sg must be initialised                                                    *
*                                                                            *
* If digraph==FALSE, loops are illegal and (x,y) => (y,x)                    *
* If prompt==TRUE, prompts are written to PROMPTFILE.                        *
* linelength is a limit on the number of characters per line caused by '?'   *
* A value of linelength <= 0 dictates no line breaks at all.                 *
* labelorg is used.                                                          *
*                                                                            *
*****************************************************************************/

void
readgraph_sg(FILE *f, sparsegraph *sg, boolean digraph, boolean prompt,
             int linelength, int n)
{
    int i,j,k,vv,ww,c;
    boolean neg,done;
    int *d,*e,*evi;
    echunk *ec,*ecnext,*ec_end;
    size_t ned,*v,iec,iec_end;

    sg->nv = n;
    DYNALLOC1(size_t,sg->v,sg->vlen,n,"malloc");
    DYNALLOC1(int,sg->d,sg->dlen,n,"malloc");
    DYNFREE(sg->w,sg->wlen);
    v = sg->v;
    d = sg->d;

    for (i = 0; i < n; ++i) d[i] = 0;

    ec = &first_echunk;
    iec = 0;
    vv = 0;
    neg = done = FALSE;

    while (!done)
    {
        GETNWC(c,f);
        if (ISDIGIT(c))
        {
            ungetc(c,f);
            readinteger(f,&ww);
            ww -= labelorg;
            if (neg)
            {
                neg = FALSE;
                if (ww < 0 || ww >= n || (!digraph && ww == vv))
                    fprintf(ERRFILE,"illegal edge (%d,%d) ignored\n\n",
                            vv+labelorg,ww+labelorg);
                else
                {
                    if (iec == ECHUNKSIZE)
                    {
                        if (!ec->next)
                        {
                            ecnext = (echunk*)ALLOCS(1,sizeof(echunk));
                            if (!ecnext) alloc_error("malloc");
                            ecnext->next = NULL;
                            ec->next = ecnext;
                        }
                        ec = ec->next;
                        iec = 0;
                    }
                    ec->edge[iec++] = vv;
                    ec->edge[iec++] = -1 - ww;
                    ++d[vv];
                    if (!digraph && ww != vv) ++d[ww];
                }
            }
            else
            {
                GETNWC(c,f);
                if (c == ':')
                {
                    if (ww < 0 || ww >= n)
                        fprintf(ERRFILE,
                                "illegal vertex number %d ignored\n\n",
                                ww+labelorg);
                    else
                        vv = ww;
                }
                else
                {
                    ungetc(c,f);
                    if (ww < 0 || ww >= n || (!digraph && ww == vv))
                        fprintf(ERRFILE,"illegal edge (%d,%d) ignored\n\n",
                                vv+labelorg,ww+labelorg);
                    else
                    {
                        if (iec == ECHUNKSIZE)
                        {
                            if (!ec->next)
                            {
                                ecnext = (echunk*)ALLOCS(1,sizeof(echunk));
                                if (!ecnext) alloc_error("malloc");
                                ecnext->next = NULL;
                                ec->next = ecnext;
                            }
                            ec = ec->next;
                            iec = 0;
                        }
                        ec->edge[iec++] = vv;
                        ec->edge[iec++] = ww;
                        ++d[vv];
                        if (!digraph && ww != vv) ++d[ww];
                    }
                }
            }
        }
        else switch(c)
        {
            case ';':
                neg = FALSE;
                ++vv;
                if (vv >= n) done = TRUE;
                break;
            case '?':
                neg = FALSE;
                fprintf(ERRFILE,"Command \'?\' not implemented.\n\n");
                break;
            case '\n':
                neg = FALSE;
                if (prompt) fprintf(PROMPTFILE,"%2d : ",vv+labelorg);
                break;
            case EOF:
            case '.':
                done = TRUE;
                break;
            case '-':
                neg = TRUE;
                break;
            case '!':
                do
                    c = getc(f);
                while (c != '\n' && c != EOF);
                if (c == '\n') ungetc(c,f);
                break;
            default :
                fprintf(ERRFILE,"illegal char '%c' - use '.' to exit\n\n",
                        (char)c);
        }
    }

    ned = 0;
    for (i = 0; i < n; ++i) ned += d[i];
    DYNALLOC1(int,sg->e,sg->elen,ned,"malloc");
    e = sg->e;

    v[0] = 0;
    for (i = 1; i < n; ++i) v[i] = v[i-1] + d[i-1];
    for (i = 0; i < n; ++i) d[i] = 0;

    iec_end = iec;
    ec_end = ec;

    iec = 0;
    ec = &first_echunk;

    if (ned != 0) while (TRUE)
    {
        vv = ec->edge[iec++];
        ww = ec->edge[iec++];

        if (ww >= 0)
        {
            e[v[vv]+(d[vv]++)] = ww;
            if (!digraph && ww != vv) e[v[ww] +(d[ww]++)] = vv;
        }
        else
        {
            ww = -1 - ww;
            for (i = 0; i < d[vv]; ++i)
                if (e[v[vv]+i] == ww) break;
            if (i < d[vv])
            {
                e[v[vv]+i] = e[v[vv]+d[vv]-1];
                --d[vv];
            }
            if (!digraph && ww != vv)
            {
                for (i = 0; i < d[ww]; ++i)
                    if (e[v[ww]+i] == vv) break;
                if (i < d[ww])
                {
                    e[v[ww]+i] = e[v[ww]+d[ww]-1];
                    --d[ww];
                }
            }
        }
        if (iec == iec_end && ec == ec_end) break;
        if (iec == ECHUNKSIZE) { iec = 0; ec = ec->next; }
    }

    sortlists_sg(sg);

    ned = 0;
    for (i = 0; i < n; ++i)
    {
        if (d[i] > 1)
        {
            evi = e + v[i];
            j = 1;
            for (k = 1; k < d[i]; ++k)
                if (evi[k] != evi[j-1]) evi[j++] = evi[k];
            d[i] = j;
        }
        ned += d[i];
    }
    sg->nde = ned; 
}

/*****************************************************************************
*                                                                            *
*  readgraph_swg(f,sg,digraph,prompt,linelength,n) reads a sparse weighted   *
*  graph g from f.                                                           *
*  Commands: (There is always a "current vertex" v, initially labelorg;      *
*             n is an unsigned integer, w is a weight.)                      *
*  n  : add edge (v,n)                                                       *
*  -n : delete edge (v,n)                                                    *
*  n: : set v := n, and exit if v >= n.                                      *
*  ?  : type neighbours of vertex v     ** NOT IMPLEMENTED **                *
*  ;  : increment v, and exit if v >= n.                                     *
*  .  : exit                                                                 *
*  !  : skip rest of input line                                              *
*  w# : set the weight for the next edge only                                *
*  W# : set the weight from now on                                           *
*  sg must be initialised                                                    *
*                                                                            *
* If digraph==FALSE, loops are illegal and (x,y) => (y,x)                    *
* For digraphs, an unspecified opposite edge has weight SG_MINWEIGHT         *
* If edges are inserted more than once, the largest weight counts.           *
* If prompt==TRUE, prompts are written to PROMPTFILE.                        *
* linelength is a limit on the number of characters per line caused by '?'   *
* A value of linelength <= 0 dictates no line breaks at all.                 *
* labelorg is used.                                                          *
*                                                                            *
*****************************************************************************/

void
readgraph_swg(FILE *f, sparsegraph *sg, boolean digraph, boolean prompt,
             int linelength, int n)
{
    int i,j,k,vv,ww,c;
    boolean neg,done;
    int *d,*e,*evi;
    echunkw *ec,*ecnext,*ec_end;
    size_t ned,*v,iec,iec_end;
    sg_weight *wt,currwt,defwt,*wvi;

    sg->nv = n;
    DYNALLOC1(size_t,sg->v,sg->vlen,n,"malloc");
    DYNALLOC1(int,sg->d,sg->dlen,n,"malloc");
    v = sg->v;
    d = sg->d;
    wt = sg->w;

    for (i = 0; i < n; ++i) d[i] = 0;

    defwt = currwt = 1; 
    ec = &first_echunkw;
    iec = 0;
    vv = 0;
    neg = done = FALSE;

    while (!done)
    {
        GETNWC(c,f);
        if (ISDIGIT(c))
        {
            ungetc(c,f);
            readinteger(f,&ww);
            ww -= labelorg;
            if (neg)
            {
                neg = FALSE;
                if (ww < 0 || ww >= n || (!digraph && ww == vv))
                    fprintf(ERRFILE,"illegal edge (%d,%d) ignored\n\n",
                            vv+labelorg,ww+labelorg);
                else
                {
                    if (iec == ECHUNKSIZE)
                    {
                        if (!ec->next)
                        {
                            ecnext = (echunkw*)ALLOCS(1,sizeof(echunkw));
                            if (!ecnext) alloc_error("malloc");
                            ecnext->next = NULL;
                            ec->next = ecnext;
                        }
                        ec = ec->next;
                        iec = 0;
                    }
                    ec->edge[iec].v1 = vv;
                    ec->edge[iec].v2 = -1 - ww;
                    ec->edge[iec].wt = currwt;
		    ++iec;
                    currwt = defwt;
                    ++d[vv];
                    if (ww != vv) ++d[ww];
                }
            }
            else
            {
                GETNWC(c,f);
                if (c == ':')
                {
                    if (ww < 0 || ww >= n)
                        fprintf(ERRFILE,
                                "illegal vertex number %d ignored\n\n",
                                ww+labelorg);
                    else
                        vv = ww;
                }
                else
                {
                    ungetc(c,f);
                    if (ww < 0 || ww >= n || (!digraph && ww == vv))
                        fprintf(ERRFILE,"illegal edge (%d,%d) ignored\n\n",
                                vv+labelorg,ww+labelorg);
                    else
                    {
                        if (iec == ECHUNKSIZE)
                        {
                            if (!ec->next)
                            {
                                ecnext = (echunkw*)ALLOCS(1,sizeof(echunkw));
                                if (!ecnext) alloc_error("malloc");
                                ecnext->next = NULL;
                                ec->next = ecnext;
                            }
                            ec = ec->next;
                            iec = 0;
                        }
                        ec->edge[iec].v1 = vv;
                        ec->edge[iec].v2 = ww;
                        ec->edge[iec].wt = currwt;
			++iec;
                        currwt = defwt;
                        ++d[vv];
                        if (ww != vv) ++d[ww];
                    }
                }
            }
        }
        else switch(c)
        {
            case ';':
                neg = FALSE;
                ++vv;
                if (vv >= n) done = TRUE;
                break;
            case '?':
                neg = FALSE;
                fprintf(ERRFILE,"Command \'?\' not implemented.\n\n");
                break;
            case '\n':
                neg = FALSE;
                if (prompt) fprintf(PROMPTFILE,"%2d : ",vv+labelorg);
                break;
            case EOF:
            case '.':
                done = TRUE;
                break;
            case '-':
                neg = TRUE;
                break;
            case 'w':
                readinteger(f,&currwt);
		if (currwt <= SG_MINWEIGHT)
		{
		    fprintf(ERRFILE,"Weight too small\n\n");
		    currwt = 1;
		}
		break;
            case 'W':
                readinteger(f,&currwt);
		if (currwt <= SG_MINWEIGHT)
		{
		    fprintf(ERRFILE,"Weight too small\n\n");
		    currwt = 1;
		}
		defwt = currwt;
		break;
            case '!':
                do
                    c = getc(f);
                while (c != '\n' && c != EOF);
                if (c == '\n') ungetc(c,f);
                break;
            default :
                fprintf(ERRFILE,"illegal char '%c' - use '.' to exit\n\n",
                        (char)c);
        }
    }

    ned = 0;
    for (i = 0; i < n; ++i) ned += d[i];
    DYNALLOC1(int,sg->e,sg->elen,ned,"malloc");
    DYNALLOC1(sg_weight,sg->w,sg->wlen,ned,"malloc");
    e = sg->e;
    wt = sg->w;

    v[0] = 0;
    for (i = 1; i < n; ++i) v[i] = v[i-1] + d[i-1];
    for (i = 0; i < n; ++i) d[i] = 0;

    iec_end = iec;
    ec_end = ec;

    iec = 0;
    ec = &first_echunkw;

    if (ned != 0) while (TRUE)
    {
        vv = ec->edge[iec].v1;
        ww = ec->edge[iec].v2;
        currwt = ec->edge[iec].wt;
        ++iec;

        if (ww >= 0)
        {
            e[v[vv]+d[vv]] = ww;
            wt[v[vv]+d[vv]] = currwt;
	    ++d[vv];
            if (ww != vv)
            {
		e[v[ww]+d[ww]] = vv;
		wt[v[ww]+d[ww]] = (digraph ? SG_MINWEIGHT : currwt);
		++d[ww];
	    }
        }
        else
        {
            ww = -1 - ww;
            for (i = 0; i < d[vv]; ++i)
                if (e[v[vv]+i] == ww) break;
            if (i < d[vv])
            {
                e[v[vv]+i] = e[v[vv]+d[vv]-1];
                wt[v[vv]+i] = wt[v[vv]+d[vv]-1];
                --d[vv];
            }
            if (ww != vv)
            {
                for (i = 0; i < d[ww]; ++i)
                    if (e[v[ww]+i] == vv) break;
                if (i < d[ww])
                {
                    e[v[ww]+i] = e[v[ww]+d[ww]-1];
                    wt[v[ww]+i] = wt[v[ww]+d[ww]-1];
                    --d[ww];
                }
            }
        }
        if (iec == iec_end && ec == ec_end) break;
        if (iec == ECHUNKSIZE) { iec = 0; ec = ec->next; }
    }

    sortlists_sg(sg);

    ned = 0;
    for (i = 0; i < n; ++i)
    {
        if (d[i] > 1)
        {
            evi = e + v[i];
	    wvi = wt + v[i];
            j = 1;
            for (k = 1; k < d[i]; ++k)
            {
                if (evi[k] != evi[j-1])
                {
		    evi[j] = evi[k];
		    wvi[j] = wvi[k];
		    ++j;
		}
		else if (wvi[k] > wvi[j-1])
		    wvi[j-1] = wvi[k];
            }
            d[i] = j;
        }
        ned += d[i];
    }
    sg->nde = ned; 
}

/*****************************************************************************
*                                                                            *
*  putgraph(f,g,linelength,m,n) writes a list of the edges of g to f         *
*  using at most linelength characters per line (excluding '\n').            *
*  A value of linelength <= 0 dictates no line breaks at all within the      *
*    list for each vertex.                                                   *
*  labelorg is used.                                                         *
*                                                                            *
*  FUNCTIONS CALLED: putset()                                                *
*                                                                            *
*****************************************************************************/

void
putgraph(FILE *f, graph *g, int linelength, int m, int n)
{
    int i,curlen;
    set *pg;

    for (i = 0, pg = g; i < n; ++i, pg += M)
    {
        fprintf(f,"%3d : ",i+labelorg);
        curlen = 7;
        putset(f,pg,&curlen,linelength,M,FALSE);
        fprintf(f,";\n");
    }
}

/*****************************************************************************
*                                                                            *
*  putgraph_sg(f,sg,linelength) writes a list of the edges of g to f         *
*  using at most linelength characters per line (excluding '\n').            *
*  A value of linelength <= 0 dictates no line breaks at all within the      *
*    list for each vertex.                                                   *
*  labelorg is used.                                                         *
*                                                                            *
*****************************************************************************/

void
putgraph_sg(FILE *f, sparsegraph *sg, int linelength)
{
    int i,n,curlen,slen;
    int *d,*e;
    size_t *v,j;
    sg_weight *wt;
    char s[60];

    n = sg->nv;
    SWG_VDE(sg,v,d,e,wt);

    for (i = 0; i < n; ++i)
    {
        fprintf(f,"%3d : ",i+labelorg);
        curlen = 7;

        for (j = v[i]; j < v[i]+d[i]; ++j)
        {
	    if (wt && wt[j] != 1)
	    {
		s[0] = 'w';
		if (wt[j] == SG_MINWEIGHT)
		{
		    s[1] = 'X';
		    s[2] = ' ';
		    slen = 3;
		}
		else
		{
		    slen = 2 + itos(wt[j],s+1);
		    s[slen-1] = ' ';
		}
                slen += itos(e[j]+labelorg,s+slen);
            }
            else
                slen = itos(e[j]+labelorg,s);

            if (linelength > 0 && curlen + slen + 1 > linelength)
            {
                putstring(f,"\n  ");
                curlen = 2;
            }
            PUTC(' ',f);
            putstring(f,s);
            curlen += slen + 1;
        }
        putstring(f,";\n");
    }
}

/*****************************************************************************
*                                                                            *
*  putmapping(f,lab1,org1,lab2,org2,linelength,n) writes n pairs             *
*  (lab1[i]+org1)-(lab2[i]+org2) to file f in increasing order of lab1[i].   *
*  lab1 and lab2 are assumed to be permutations.  At most linelength         *
*  characters (excluding '\n') are written per line.                         *
*  A value of linelength <= 0 dictates no line breaks at all.                *
*                                                                            *
*  FUNCTIONS CALLED: itos(),putstring()                                      *
*                                                                            *
*****************************************************************************/

void
putmapping(FILE *f, int *lab1, int org1,int *lab2, int org2,
       int linelength, int n)
{
    int i,curlen,slen;
    char s[60];

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"putmapping");
#endif

    for (i = 0; i < n; ++i) workperm[lab1[i]] = lab2[i];

    curlen = 0;
    for (i = 0; i < n; ++i)
    {
        slen = itos(i+org1,s);
        s[slen++] = '-';
        slen += itos(workperm[i]+org2,&s[slen]);
        if (linelength > 0 && curlen + slen + 1 > linelength)
        {
            putstring(f,"\n  ");
            curlen = 2;
        }
        PUTC(' ',f);
        putstring(f,s);
        curlen += slen + 1;
    }
    PUTC('\n',f);
}

/*****************************************************************************
*                                                                            *
*  putorbits(f,orbits,linelength,n) writes the orbits to f as lists          *
*  of integers separated by semicolons.  No more than linelength characters  *
*  (excluding '\n') are written per line.                                    *
*  A value of linelength <= 0 dictates no line breaks at all.                *
*  labelorg is used.                                                         *
*                                                                            *
*  FUNCTIONS CALLED: itos(),putset()                                         *
*                                                                            *
*****************************************************************************/

void
putorbits(FILE *f, int *orbits, int linelength, int n)
{
    int i,j;
    int m,curlen,sz,slen;
    char s[20];

    m = SETWORDSNEEDED(n);
#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"putorbits");
    DYNALLOC1(set,workset,workset_sz,m,"putorbits");
#endif

    for (i = n; --i >= 0;) workperm[i] = 0;
    for (i = n; --i >= 0;)
        if ((j = orbits[i]) < i)
        {
            workperm[i] = workperm[j];
            workperm[j] = i;
        }

    curlen = 0;
    for (i = 0; i < n; ++i)
        if (orbits[i] == i)
        {
	    sz = 0;
            EMPTYSET(workset,m);
            j = i;
            do
            {
                ADDELEMENT(workset,j);
                j = workperm[j];
		++sz;
            }
            while (j > 0);
            putset(f,workset,&curlen,linelength-1,m,TRUE);
	    if (sz > 1)
	    {
	        s[0] = ' ';
		s[1] = '(';
		slen = 2 + itos(sz,s+2);
		s[slen++] = ')';
		s[slen] = '\0';
		if (linelength > 0 && curlen + slen + 1 >= linelength)
        	{
            	    fprintf(f,"\n   ");
                    curlen = 3;
        	}
		fprintf(f,"%s",s);
		curlen += slen;
	    }

            PUTC(';',f);
            ++curlen;
        }
    PUTC('\n',f);
}


/*****************************************************************************
*                                                                            *
*  putorbitsplus(f,orbits,linelength,n) is the same as                       *
*  putorbits(f,orbits,linelength,n) except that the first element of each    *
*  orbit is written bold.  This only works if output is to a device that     *
*  interprets ANSI controls.                                                 *
*                                                                            *
*  FUNCTIONS CALLED: itos(),putset()                                         *
*                                                                            *
*****************************************************************************/

void
putorbitsplus(FILE *f, int *orbits, int linelength, int n)
{
    int i,j;
    int m,curlen,sz,slen;
    char s[20];

    m = SETWORDSNEEDED(n);
#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"putorbits");
    DYNALLOC1(set,workset,workset_sz,m,"putorbits");
#endif

    for (i = n; --i >= 0;) workperm[i] = 0;
    for (i = n; --i >= 0;)
        if ((j = orbits[i]) < i)
        {
            workperm[i] = workperm[j];
            workperm[j] = i;
        }

    curlen = 0;
    for (i = 0; i < n; ++i)
        if (orbits[i] == i)
        {
	    sz = 0;
            EMPTYSET(workset,m);
            j = i;
            do
            {
                ADDELEMENT(workset,j);
                j = workperm[j];
		++sz;
            }
            while (j > 0);
            putset_firstbold(f,workset,&curlen,linelength-1,m,TRUE);
	    if (sz > 1)
	    {
	        s[0] = ' ';
		s[1] = '(';
		slen = 2 + itos(sz,s+2);
		s[slen++] = ')';
		s[slen] = '\0';
		if (linelength > 0 && curlen + slen + 1 >= linelength)
        	{
            	    fprintf(f,"\n   ");
                    curlen = 3;
        	}
		fprintf(f,"%s",s);
		curlen += slen;
	    }

            PUTC(';',f);
            ++curlen;
        }
    PUTC('\n',f);
}

/*****************************************************************************
*                                                                            *
*  putquotient(f,g,lab,ptn,level,linelength,m,n) writes the quotient matrix  *
*  of g defined by the partition at the given level.  Precisely, for each    *
*  cell W, it writes the number w of the least vertex in W, then the size    *
*  of W, then the number of times w is joined FROM each cell.  A complete    *
*  join is written as "*", while an empty join is written as "-".  No more   *
*  than linelength  characters (excluding '\n') are written per line unless  *
*  linelength is very small.  A value of linelength <= 0 dictates no line    *
*  breaks at all.   labelorg is used.                                        *
*                                                                            *
*  FUNCTIONS CALLED: itos()                                                  *
*                                                                            *
*****************************************************************************/

void
putquotient(FILE *f, graph *g, int *lab, int *ptn, int level,
        int linelength, int m, int n)
{
    int i;
    char s[50];
    int ic,curlen,v,w,cell1,cell2,numcells,jc,csize,k;
    set *gw;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"putquotient");
    DYNALLOC1(set,workset,workset_sz,m,"putquotient");
#endif

    numcells = 0;
    for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
    {
        for (cell2 = cell1; ptn[cell2] > level; ++cell2) {}
        w = lab[cell1];
        for (i = cell1 + 1; i <= cell2; ++i)
            if (lab[i] < w) w = lab[i];
        workperm[numcells++] = w;
    }

    for (ic = cell1 = 0; ic < numcells; ++ic, cell1 = cell2 + 1)
    {
        for (cell2 = cell1; ptn[cell2] > level; ++cell2) {}
        EMPTYSET(workset,M);
        for (i = cell1; i <= cell2; ++i) ADDELEMENT(workset,lab[i]);
        v = workperm[ic];
        csize = cell2 - cell1 + 1;
        if (v + labelorg < 10)
        {
            s[0] = ' ';
            curlen = 1;
        }
        else
            curlen = 0;
        curlen += itos(v+labelorg,&s[curlen]);
        s[curlen++] = '[';
        curlen += itos(csize,&s[curlen]);
        fprintf(f,"%s",s);
        if (csize < 10)
        {
            fprintf(f,"]  :");
            curlen += 4;
        }
        else
        {
            fprintf(f,"] :");
            curlen += 3;
        }

        for (jc = 0; jc < numcells; ++jc)
        {
            w = workperm[jc];
            gw = GRAPHROW(g,w,m);
            k = setinter(gw,workset,M);
            if (k == 0 || k == csize)
            {
                if (linelength > 0 && curlen + 2 > linelength)
                {
                    fprintf(f,"\n    ");
                    curlen = 4;
                }
                if (k == 0) fprintf(f," -");
                else        fprintf(f," *");
                curlen += 2;
            }
            else
            {
                i = itos(k,s);
                if (linelength > 0 && curlen + i + 1 > linelength)
                {
                    fprintf(f,"\n    ");
                    curlen = 4;
                }
                fprintf(f," %s",s);
                curlen += i + 1;
            }
        }
        fprintf(f,"\n");
    }
}

/*****************************************************************************
*                                                                            *
*  putquotient_sg(f,g,lab,ptn,level,linelength) writes the quotient matrix   *
*  of g defined by the partition at the given level.  Precisely, for each    *
*  cell W, it writes the number w of the least vertex in W, then the size    *
*  of W, then the number of times w is joined FROM each cell.  A complete    *
*  join is written as "*", while an empty join is written as "-".  No more   *
*  than linelength  characters (excluding '\n') are written per line unless  *
*  linelength is very small.  A value of linelength <= 0 dictates no line    *
*  breaks at all.   labelorg is used.                                        *
*                                                                            *
*  Weughts are ignored.                                                      *
*                                                                            *
*****************************************************************************/

void
putquotient_sg(FILE *f, sparsegraph *g, int *lab, int *ptn,
						int level, int linelength)
{
    int i,m,n;
    char s[50];
    int ic,curlen,v,w,cell1,cell2,numcells,jc,csize,k;
    int *dd,*ee;
    size_t *vv,j;

    n = g->nv;
    m = SETWORDSNEEDED(n);
    SG_VDE(g,vv,dd,ee);

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"putquotient");
    DYNALLOC1(set,workset,workset_sz,m,"putquotient");
#endif

    numcells = 0;
    for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
    {
        for (cell2 = cell1; ptn[cell2] > level; ++cell2) {}
        w = lab[cell1];
        for (i = cell1 + 1; i <= cell2; ++i)
            if (lab[i] < w) w = lab[i];
        workperm[numcells++] = w;
    }

    for (ic = cell1 = 0; ic < numcells; ++ic, cell1 = cell2 + 1)
    {
        for (cell2 = cell1; ptn[cell2] > level; ++cell2) {}
        EMPTYSET(workset,M);
        for (i = cell1; i <= cell2; ++i) ADDELEMENT(workset,lab[i]);
        v = workperm[ic];
        csize = cell2 - cell1 + 1;
        if (v + labelorg < 10)
        {
            s[0] = ' ';
            curlen = 1;
        }
        else
            curlen = 0;
        curlen += itos(v+labelorg,&s[curlen]);
        s[curlen++] = '[';
        curlen += itos(csize,&s[curlen]);
        fprintf(f,"%s",s);
        if (csize < 10)
        {
            fprintf(f,"]  :");
            curlen += 4;
        }
        else
        {
            fprintf(f,"] :");
            curlen += 3;
        }

        for (jc = 0; jc < numcells; ++jc)
        {
            w = workperm[jc];
            k = 0;
	    for (j = vv[w]; j < vv[w]+dd[w]; ++j)
		if (ISELEMENT(workset,ee[j])) ++k;

            if (k == 0 || k == csize)
            {
                if (linelength > 0 && curlen + 2 > linelength)
                {
                    fprintf(f,"\n    ");
                    curlen = 4;
                }
                if (k == 0) fprintf(f," -");
                else        fprintf(f," *");
                curlen += 2;
            }
            else
            {
                i = itos(k,s);
                if (linelength > 0 && curlen + i + 1 > linelength)
                {
                    fprintf(f,"\n    ");
                    curlen = 4;
                }
                fprintf(f," %s",s);
                curlen += i + 1;
            }
        }
        fprintf(f,"\n");
    }
}

/*****************************************************************************
*                                                                            *
*  putptn(f,lab,ptn,level,linelength,n) writes the partition at the given    *
*  level as sorted lists of integers separated by semicolons.  No more than  *
*  linelength characters (excluding '\n') are written per line.              *
*  A value of linelength <= 0 dictates no line breaks at all.                *
*  labelorg is used.                                                         *
*                                                                            *
*****************************************************************************/

void
putptn(FILE *f, int *lab, int *ptn, int level, int linelength, int n)
{
    int i;
    int curlen,m;

    m = SETWORDSNEEDED(n);
#if !MAXN
    DYNALLOC1(set,workset,workset_sz,m,"putptn");
#endif

    PUTC('[',f);
    curlen = 1;
    i = 0;
    while (i < n)
    {
        EMPTYSET(workset,m);
        while (TRUE)
        {
            ADDELEMENT(workset,lab[i]);
            if (ptn[i] > level) ++i;
            else                break;
        }
        putset(f,workset,&curlen,linelength-2,m,TRUE);
        if (i < n-1)
        {
            fprintf(f," |");
            curlen += 2;
        }
        ++i;
    }
    fprintf(f," ]\n");
}

/*****************************************************************************
*                                                                            *
*  putcanon(f,canonlab,canong,linelength,m,n) writes the label canonlab      *
*  and the graph canong to f, using at most linelength characters            *
*  (excluding '\n') per line.   labelorg is used.                            *
*  A value of linelength <= 0 dictates no line breaks at all.                *
*                                                                            *
*****************************************************************************/

void
putcanon(FILE *f, int *canonlab, graph *canong, int linelength, int m, int n)
{
    int i;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"putcanon");
#endif

    for (i = 0; i < n; ++i) workperm[i] = canonlab[i];
    writeperm(f,workperm,TRUE,linelength,n);
    putgraph(f,canong,linelength,m,n);
}

/*****************************************************************************
*                                                                            *
*  putcanon_sg(f,canonlab,canong,linelength) writes the label canonlab       *
*  and the graph canong to f, using at most linelength characters            *
*  (excluding '\n') per line.   labelorg is used.                            *
*  A value of linelength <= 0 dictates no line breaks at all.                *
*                                                                            *
*****************************************************************************/

void
putcanon_sg(FILE *f, int *canonlab, sparsegraph *canong, int linelength)
{
    int i,n;

    n = canong->nv;
#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"putcanon");
#endif

    for (i = 0; i < n; ++i) workperm[i] = canonlab[i];
    writeperm(f,workperm,TRUE,linelength,n);
    putgraph_sg(f,canong,linelength);
}

/*****************************************************************************
*                                                                            *
*  readptn(f,lab,ptn,&numcells,prompt,n) reads a partition from f            *
*  and establishes it in (lab,ptn).                                          *
*  The format can be just a number, which is fixed alone, or an arbitrary    *
*  partition [...|...|...].  Ranges x:y can be used.                         *
*  labelorg is used.                                                         *
*                                                                            *
*****************************************************************************/

void
readptn(FILE *f, int *lab, int *ptn, int *numcells, boolean prompt, int n)
{
    int i,j;
    int c,v1,v2,m;

    m = SETWORDSNEEDED(n);
#if !MAXN
    DYNALLOC1(set,workset,workset_sz,m,"readptn");
#endif

    GETNW(c,f);
    if (c == '=') GETNW(c,f);
    if (ISDIGIT(c))
    {
        ungetc(c,f);
        readinteger(f,&v1);
        v1 -= labelorg;
        if (v1 >= 0 && v1 < n)
            fixit(lab,ptn,numcells,v1,n);
        else
        {
            fprintf(ERRFILE,"vertex out of range (%d), fixing nothing\n\n",
                    v1+labelorg);
            unitptn(lab,ptn,numcells,n);
        }
        return;
    }
    if (c != '[')
    {
        ungetc(c,f);
        fprintf(ERRFILE,"illegal partition, fixing nothing\n\n");
        unitptn(lab,ptn,numcells,n);
        return;
    }
    EMPTYSET(workset,m);
    *numcells = 0;
    for (i = 0; i < n; ++i) ptn[i] = NAUTY_INFINITY;
    i = 0;
    j = -1;
    while (TRUE)
    {
        GETNWC(c,f);
        if (ISDIGIT(c))
        {
            ungetc(c,f);
            readinteger(f,&v1);
            v1 -= labelorg;
            GETNWC(c,f);
            if (c == ':')
                if (!readinteger(f,&v2))
                {
                    fprintf(ERRFILE,"unfinished range\n\n");
                    v2 = v1;
                }
                else
                    v2 -= labelorg;
            else
            {
                ungetc(c,f);
                v2 = v1;
            }
            while (v1 <= v2)
            {
                if (v1 < 0 || v1 >= n || ISELEMENT(workset,v1))
                    fprintf(ERRFILE,"illegal or repeated number : %d\n\n",
                            v1+labelorg);
                else
                {
                    ADDELEMENT(workset,v1);
                    lab[++j] = v1;
                }
                ++v1;
            }
        }
        else if (c == '|' || c == ']' || c == EOF)
        {
            if (j >= i)
            {
                ++*numcells;
                ptn[j] = 0;
            }
            if (c == '|')
                i = j + 1;
            else if (j == n - 1)
                return;
            else
            {
                i = j + 1;
                ++*numcells;
                for (j = 0; j < n; ++j)
                    if (!ISELEMENT(workset,j)) lab[i++] = j;
                ptn[n-1] = 0;
                return;
            }
        }
        else if (c == '\n')
        {
            if (prompt) fprintf(PROMPTFILE,"] ");
        }
        else
            fprintf(ERRFILE,"illegal character '%c' in partition\n\n",c);
    }
}

/*****************************************************************************
*                                                                            *
*  unitptn(lab,ptn,&numcells,n) establishes the partition with one cell.     *
*                                                                            *
*****************************************************************************/

void
unitptn(int *lab,int *ptn, int *numcells, int n)
{
    int i;

    for (i = 0; i < n; ++i)
    {
        lab[i] = i;
        ptn[i] = NAUTY_INFINITY;
    }
    ptn[n-1] = 0;
    *numcells = 1;
}

/*****************************************************************************
*                                                                            *
*  individualise(lab,ptn,level,v,&pos,&numcells,n) individualises vertex v.  *
*  numcells is updated and the position of the possibly-new singleton is     *
*  returned in pos.                                                          *
*                                                                            *
*****************************************************************************/

void
individualise(int *lab,int *ptn, int level,
				int v, int *pos, int *numcells, int n)
{
    int i,j;

    for (i = 0; i < n; ++i) if (lab[i] == v) break;

    for (j = i; j > 0 && ptn[j-1] > level; --j) {};

    *pos = j;
    if (ptn[j] <= level) return;   /* individual already */

    lab[i] = lab[j];
    lab[j] = v;
    ptn[j] = level;
    ++*numcells;
}

/*****************************************************************************
*                                                                            *
*  cellstarts(ptn,level,cell,m,n) sets the set cell to contain the indices   *
*  of the starts in ptn of the partition at level level.                     *
*                                                                            *
*****************************************************************************/

void
cellstarts(int *ptn, int level, set *cell, int m, int n)
{
    int i;

    EMPTYSET(cell,m);
    i = 0;
    while (i < n)
    {
        ADDELEMENT(cell,i);
        while (ptn[i] > level) ++i;
        ++i;
    }
}

/*****************************************************************************
*                                                                            *
*  fixit(lab,ptn,&numcells,fixedvertex,n) establishes the partition          *
*  with one cell {fixedvertex} and all the other vertices (if any) in        *
*  another cell.                                                             *
*                                                                            *
*****************************************************************************/

void
fixit(int *lab, int *ptn, int *numcells, int fixedvertex, int n)
{
    int i;

    for (i = 1; i < n; ++i)
    {
        lab[i] = i;
        ptn[i] = 1;
    }

    lab[0] = fixedvertex;
    lab[fixedvertex] = 0;
    ptn[0] = 0;
    ptn[n-1] = 0;
    if (n == 1) *numcells = 1;
    else        *numcells = 2;
}

/*****************************************************************************
*                                                                            *
*  sethash(s,n,seed,key) is a function whose value depends only on the       *
*  set s, a long seed, and an integer key.  It is intended to be independent *
*  of the word size provided long ints have at least 32 bits, and also       *
*  independent of m.  n is the underlying universal set size.                *
*  31 bits of seed and 15 bits of key are significant.                       *
*  The result is in 0..2^31-1.                                               *
*                                                                            *
*****************************************************************************/

long
sethash(set *s, int n, long seed, int key)
{
    int i,j,lsh,rsh;
    unsigned long l,res,lshmask,salt;
    setword si;

    lsh = key & 0xF;
    rsh = 28 - lsh;
    salt = (key >> 4) & 0x7FFL;
    res = seed & 0x7FFFFFFFUL;
    lshmask = (1UL << lsh) - 1;

    j = 0;
    for (i = 0; ; ++i)
    {
        si = s[i];
        l = SWCHUNK0(si);
        res = (((res << lsh) ^ ((res >> rsh) & lshmask) ^ l) + salt) 
                                                            & 0x7FFFFFFFUL;
	res = FUZZ1(res);
        if ((j += 16) >= n) break;
#if WORDSIZE > 16
        l = SWCHUNK1(si);
        res = (((res << lsh) ^ ((res >> rsh) & lshmask) ^ l) + salt) 
                                                           & 0x7FFFFFFFUL;
	res = FUZZ1(res);
        if ((j += 16) >= n) break;
#if WORDSIZE == 64    
        l = SWCHUNK2(si);
        res = (((res << lsh) ^ ((res >> rsh) & lshmask) ^ l) + salt) 
                                                            & 0x7FFFFFFFUL;
	res = FUZZ1(res);
        if ((j += 16) >= n) break;
        l = SWCHUNK3(si);
        res = (((res << lsh) ^ ((res >> rsh) & lshmask) ^ l) + salt) 
                                                            & 0x7FFFFFFFUL;
	res = FUZZ1(res);
        if ((j += 16) >= n) break;
#endif
#endif
    }

    return res;
}

/*****************************************************************************
*                                                                            *
*  listhash(x,nx,key) is a function whose value depends on the SET of values *
*  in the first 'nx' entries of the array 'x', and the value of key.         *
*  Machine-independent if long ints have at least 32 bits, otherwise not.    *
*  The result is in 0..2^31-1.                                               *
*                                                                            *
*****************************************************************************/

long
listhash(int *x, int nx, long key)
{
    unsigned long lkey,val,accum;
    int i;

    lkey = (unsigned long)key & 0x7FFFFFFFUL;
    accum = nx;
    
    for (i = 0; i < nx; ++i)
    {
	val = (unsigned long)x[i] & 0x7FFFFFFFUL;
	val = (val + lkey) & 0x7FFFFFFFUL;
	accum += FUZZ1(val);
    }

    return  accum & 0x7FFFFFFFUL;
}

/*****************************************************************************
*                                                                            *
*  hashgraph_sg(sg,key) is a function whose value depends on the sparse      *
*  graph or digraph sg.                                                      *
*  Machine-independent if long ints have at least 32 bits, otherwise not.    *
*  The result is in 0..2^31-1.                                               *
*                                                                            *
*****************************************************************************/

long
hashgraph_sg(sparsegraph *sg, long key)
{
    int n,i;
    int *d,*e;
    size_t *v;
    unsigned long val,accum;

    CHECK_SWG(sg,"hashgraph_sg");
    accum = n = sg->nv;

    SG_VDE(sg,v,d,e);

    for (i = 0; i < n; ++i)
        if (d[i] == 0)
	    accum += FUZZ1(i);
        else
        {
	    accum = (accum>>7) | ((accum<<24)&0x7FFFFFFFUL);
	    val = listhash(e+v[i],d[i],key);
	    val = (val + i) & 0x7FFFFFFFUL;
	    accum += FUZZ2(val);
	}

   return (long)(accum & 0x7FFFFFFFUL);
}

/*****************************************************************************
*                                                                            *
*  hashgraph(g,m,n,key) is a function whose value depends on the             *
*  graph or digraph sg.                                                      *
*  Machine-independent if long ints have at least 32 bits, otherwise not.    *
*  The result is in 0..2^31-1.                                               *
*                                                                            *
*****************************************************************************/

long
hashgraph(graph *g, int m, int n, long key)
{
    int i;
    set *gi;
    unsigned long val,accum;

    accum = n;
    for (i = 0, gi = g; i < n; ++i, gi += m)
    {
	accum = (accum>>12) | ((accum<<19)&0x7FFFFFFFUL);
	val = sethash(gi,n,key,(key&0xFL)+i);
	val = (val + i) & 0x7FFFFFFFUL;
        accum += FUZZ2(val);
    }

    return (long)(accum & 0x7FFFFFFFUL);
}

/*****************************************************************************
*                                                                            *
*  hash(setarray,length,key) is a function whose value depends only on the   *
*  first 'length' entries of the array 'setarray', and the value of key.     *
*  key should be in the range 1..31 and preferably odd.                      *
*  This works best if long ints have at least 32 bits, but it's valid anyway.*
*  Not machine-indpendent!  Use sethash() in preference.                     *
*                                                                            *
*****************************************************************************/

long
hash(set *setarray, long length, int key)
{
    long code;
    set *sptr;

    code = length;
    sptr = setarray + length;

    while (--sptr >= setarray)
        code = (code<<key) ^ ((code>>(32-key)) + *sptr);

    return code;
}

/*****************************************************************************
*                                                                            *
*  readperm is like readvperm without the last argument.  It is provided     *
*  only for backward compatibility.                                          *
*                                                                            *
*****************************************************************************/

void
readperm(FILE *f, int *perm, boolean prompt, int n)
{
    int nv;

    readvperm(f,perm,prompt,n,&nv);
}

/*****************************************************************************
*                                                                            *
*  readvperm(f,perm,prompt,n,nv) reads a permutation of order n from         *
*  f, terminated by a semicolon.  Any repeated or illegal numbers or         *
*  characters are reported then ignored.    Missing numbers are filled in    *
*  in numerical order.  A prompt is issued for each line if prompt!=FALSE.   *
*  labelorg is used.  *nv is set equal to the number of numbers actually     *
*  given.  Ranges like v1:v2 are allowed.                                    *
*                                                                            *
*****************************************************************************/

void
readvperm(FILE *f, int *perm, boolean prompt, int n, int *nv)
{
    int i;
    int m,c,v1,v2;

    m = SETWORDSNEEDED(n);
#if !MAXN
    DYNALLOC1(set,workset,workset_sz,m,"readperm");
#endif

    EMPTYSET(workset,m);

    i = 0;

    while (TRUE)
    {
        GETNWC(c,f);
        if (c == ';' || c == EOF) break;
        if (ISDIGIT(c))
        {
            ungetc(c,f);
            readinteger(f,&v1);
            v1 -= labelorg;
            GETNWC(c,f);
            if (c == ':')
                if (!readinteger(f,&v2))
                {
                    fprintf(ERRFILE,"unfinished range\n\n");
                    v2 = v1;
                }
                else
                    v2 -= labelorg;
            else
            {
                ungetc(c,f);
                v2 = v1;
            }

            if (v1 < 0 || v1 >= n || v2 >= n || v1 > v2)
            {
                if (v1 < v2)
                    fprintf(ERRFILE,
                      "illegal range in permutation : %d:%d\n\n",
                      v1+labelorg,v2+labelorg);
                else
                    fprintf(ERRFILE,
                      "illegal number in permutation : %d\n\n",
                      v1+labelorg);
            }
            else
            for (; v1 <= v2; ++v1)
            {
                if (!ISELEMENT(workset,v1))
                {
                    perm[i++] = v1;
                    ADDELEMENT(workset,v1);
                }
                else
                    fprintf(ERRFILE,
                      "repeated number in permutation : %d\n\n",
                      v1+labelorg);
            }
        }
        else
        {
            if (c == '\n' && prompt)
                fprintf(PROMPTFILE,"+ ");
            if (c != '\n')
                fprintf(ERRFILE,"bad character '%c' in permutation\n\n",
                       (char)c);
        }
    }

    *nv = i;

    for (v1 = 0; v1 < n; ++v1)
        if (!ISELEMENT(workset,v1)) perm[i++] = v1;
}

/*****************************************************************************
*                                                                            *
*  ranperm(perm,n) creates a random permutation in perm.                     *
*                                                                            *
*****************************************************************************/

void
ranperm(int *perm, int n)
{
    int i,j,t;

    for (i = n; --i >= 0; ) perm[i] = i;

    for (i = n; --i > 0; )
    {
        j = KRAN(i+1);
        t = perm[i];
        perm[i] = perm[j];
        perm[j] = t;
    }
}

/*****************************************************************************
*                                                                            *
*  relabel(g,perm,lab,workg,m,n) replaces g by g^perm, using workg as        *
*  scratch space.  If lab!=NULL, it is taken as a labelling vector and       *
*  also permuted.                                                            *
*                                                                            *
*****************************************************************************/

void
relabel(graph *g, int *lab, int *perm, graph *workg, int m, int n)
{
    long li;
    int i;

    for (li = (long)M * (long)n; --li >= 0;) workg[li] = g[li];

    updatecan(workg,g,perm,0,M,n);
    if (lab != NULL)
    {
#if !MAXN
        DYNALLOC1(int,workperm,workperm_sz,n+2,"relabel");
#endif
        for (i = 0; i < n; ++i) workperm[perm[i]] = i;
        for (i = 0; i < n; ++i) lab[i] = workperm[lab[i]];
    }
}

/*****************************************************************************
*                                                                            *
*  relabel_sg(g,perm,lab,workg,m,n) replaces g by g^perm, using workg as     *
*  scratch space.  If lab!=NULL, it is taken as a labelling vector and       *
*  also permuted.                                                            *
*                                                                            *
*****************************************************************************/

void
relabel_sg(sparsegraph *sg, int *lab, int *perm, sparsegraph *workg)
{
    int i,n;
    sparsegraph *tempsg;
    sparsegraph tmp;

    n = sg->nv;

    if (workg)
    {
        tempsg = copy_sg(sg,workg);
        updatecan_sg((graph*)tempsg,(graph*)sg,perm,0,SETWORDSNEEDED(n),n);
    }
    else
    {
        SG_INIT(tmp);
        tempsg = copy_sg(sg,&tmp);
        updatecan_sg((graph*)tempsg,(graph*)sg,perm,0,SETWORDSNEEDED(n),n);
        SG_FREE(tmp);
    }

    if (lab != NULL)
    {
#if !MAXN
        DYNALLOC1(int,workperm,workperm_sz,n+2,"relabel_sg");
#endif
        for (i = 0; i < n; ++i) workperm[perm[i]] = i;
        for (i = 0; i < n; ++i) lab[i] = workperm[lab[i]];
    }
}

/*****************************************************************************
*                                                                            *
*  sublabel(g,perm,nperm,workg,m,n) replaces g by g^perm, using workg as     *
*  scratch space.  perm is a partial vector, of length nperm, where it is    *
*  known that the elements of perm are distinct.                             *
*                                                                            *
*****************************************************************************/

void
sublabel(graph *g, int *perm, int nperm, graph *workg, int m, int n)
{
    long li;
    int i,j,k;
    int newm;
    set *gi,*wgi;

    for (li = (long)m * (long)n; --li >= 0;) workg[li] = g[li];

    newm = SETWORDSNEEDED(nperm);

    for (li = (long)newm * (long)nperm; --li >= 0;) g[li] = 0;

    for (i = 0, gi = (set*)g; i < nperm; ++i, gi += newm)
    {
        wgi = GRAPHROW(workg,perm[i],m);
        for (j = 0; j < nperm; ++j)
        {
            k = perm[j];
            if (ISELEMENT(wgi,k)) ADDELEMENT(gi,j);
        }
    }
}

/*****************************************************************************
*                                                                            *
*  countcells(ptn,level,n) finds the number of elements of ptn[0..n-1]       *
*  that are <= level.                                                        *
*                                                                            *
*****************************************************************************/

int
countcells(int *ptn, int level, int n)
{
    int i,cnt;

    cnt = 0;
    for (i = 0; i < n; ++i) if (ptn[i] <= level) ++cnt;

    return cnt;
}

/*****************************************************************************
*                                                                            *
*  subpartion(lab,ptn,n,perm,nperm) replaces the partition (lab,ptn) of      *
*  0..n-1 by the induced partition of 0..nperm-1, using the partial          *
*  ordering of 0..n-1 given in perm[0..nperm-1].                             *
*  Return the new number of cells.                                           *
*                                                                            *
*****************************************************************************/

#define DEB(str,x) fprintf(stderr,"%s=%d\n",str,x);

int
subpartition(int *lab, int *ptn, int n, int *perm, int nperm)
{
    int i,j;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"subpartition");
#endif
    for (i = 0; i < n; ++i) workperm[i] = -1;
    for (i = 0; i < nperm; ++i) workperm[perm[i]] = i;

    j = -1;
    for (i = 0; i < n; ++i)
    {
	if (workperm[lab[i]] >= 0)
	{      
	    ++j;     
	    lab[j] = workperm[lab[i]];
	    ptn[j] = ptn[i];
	}
	else if (j >= 0 && ptn[i] < ptn[j])
	    ptn[j] = ptn[i];
    }

    return countcells(ptn,0,nperm);
}

/*****************************************************************************
*                                                                            *
*  sublabel_sg(sg,perm,nperm,workg) replaces g by g^perm, using workg as     *
*  scratch space.  perm is a partial vector, of length nperm, where it is    *
*  known that the elements of perm are distinct.                             *
*                                                                            *
*****************************************************************************/

void
sublabel_sg(sparsegraph *sg, int *perm, int nperm, sparsegraph *workg)
{
    int i,j,k,n;
    size_t newnde,kk;
    sparsegraph *tempsg;
    sparsegraph tmp;
    int *d,*e;
    int *dd,*ee;
    size_t *v,*vv;

    CHECK_SWG(sg,"sublabel_sg");
    n = sg->nv;
 
#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"relabel_sg");
#endif
    for (i = 0; i < n; ++i) workperm[i] = -1;
    for (i = 0; i < nperm; ++i) workperm[perm[i]] = i;

    newnde = 0;
    SG_VDE(sg,v,d,e);

    for (i = 0; i < nperm; ++i)
    {
        j = perm[i];
        for (k = 0; k < d[j]; ++k)
            if (workperm[e[v[j]+k]] >= 0) ++newnde;
    }
 
    if (workg)
        tempsg = workg;
    else
    {
        SG_INIT(tmp);
        tempsg = &tmp;
    }
 
    SG_ALLOC(*tempsg,nperm,newnde,"sublabel_sg");
    SG_VDE(tempsg,vv,dd,ee);

    kk = 0;
    for (i = 0; i < nperm; ++i)
    {
        j = perm[i];
        vv[i] = kk;
        dd[i] = 0;
        for (k = 0; k < d[j]; ++k)
            if (workperm[e[v[j]+k]] >= 0)
            {
                ee[vv[i]+dd[i]] = workperm[e[v[j]+k]];
                ++dd[i];
            }
	kk += dd[i];
    }
    tempsg->nv = nperm;
    tempsg->nde = newnde;

    copy_sg(tempsg,sg);

    if (!workg) SG_FREE(tmp);
}

/*****************************************************************************
*                                                                            *
*  copycomment(fin,fout,delimiter) copies fin to fout until either EOF or    *
*  the character 'delimiter' is read.  The delimiter itself isn't written.   *
*  Escape sequences \n,\t,\b,\r,\f,\\,\',\",\\n are recognised.  Otherwise,  *
*  '\' is ignored.                                                           *
*                                                                            *
*****************************************************************************/

void
copycomment(FILE *fin, FILE *fout, int delimiter)
{
    int c,backslash;

    backslash = FALSE;

    while ((c = getc(fin)) != EOF && (c != delimiter || backslash))
        if (backslash)
        {
            switch (c)
            {
            case '\n':
                break;
            case 'n':
                PUTC('\n',fout);
                break;
            case 't':
                PUTC('\t',fout);
                break;
            case 'b':
                PUTC('\b',fout);
                break;
            case 'r':
                PUTC('\r',fout);
                break;
            case 'f':
                PUTC('\f',fout);
                break;
            case '\\':
                PUTC('\\',fout);
                break;
            case '\'':
                PUTC('\'',fout);
                break;
            case '"':
                PUTC('"',fout);
                break;
            default:
                PUTC(c,fout);
            }
            backslash = FALSE;
        }
        else if (c == '\\')
            backslash = TRUE;
        else
            PUTC(c,fout);
}

/*****************************************************************************
*                                                                            *
*  converse_sg(g1,g2) performs a digraph converse operation on g1,           *
*  leaving the result in g2. g2 must exist and be initialised.               *
*  If g1 is an undirected graph, g2 will be the same.                        *
*                                                                            *
*****************************************************************************/

void
converse_sg(sparsegraph *g1, sparsegraph *g2)
{
    int *e1,*d1,*e2,*d2;
    size_t *v1,*v2,j;
    int i,k,n;

    CHECK_SWG(g1,"converse_sg");
    n = g1->nv;

    SG_ALLOC(*g2,n,g1->nde,"converse_sg");
    g2->nv = n;
    g2->nde = g1->nde;
    DYNFREE(g2->w,g2->wlen);

    SG_VDE(g1,v1,d1,e1);
    SG_VDE(g2,v2,d2,e2);

    for (i = 0; i < n; ++i) d2[i] = 0;
    for (i = 0; i < n; ++i)
        for (j = v1[i]; j < v1[i]+d1[i]; ++j) ++d2[e1[j]];

    v2[0] = 0;
    for (i = 1; i < n; ++i) v2[i] = v2[i-1] + d2[i-1];
    for (i = 0; i < n; ++i) d2[i] = 0;

    for (i = 0; i < n; ++i)
        for (j = v1[i]; j < v1[i]+d1[i]; ++j)
	{
	    k = e1[j];
	    e2[v2[k] + (d2[k]++)] = i;
	}
}

/*****************************************************************************
*                                                                            *
*  complement_sg(g1,g2) sets g2 to the complement of g1.                     *
*  If g1 has loops then the loop set is also complemented; otherwise         *
*  no loops are created.  g2 must exist and be initialised.                  *
*                                                                            *
*****************************************************************************/

void
complement_sg(sparsegraph *g1, sparsegraph *g2)
{
    int *e1,*d1,*e2,*d2;
    size_t *v1,*v2,j,ndec;
    int i,l,m,n;
    int loops;

    CHECK_SWG(g1,"complement_sg");
    n = g1->nv;
    SG_VDE(g1,v1,d1,e1);

    loops = 0;
    for (i = 0; i < n; ++i)
	for (j = v1[i]; j < v1[i] + d1[i]; ++j)
	    if (e1[j] == i) ++loops;

    if (loops > 1) ndec = n*(size_t)n - g1->nde;
    else           ndec = n*(size_t)n - n - g1->nde; 
    SG_ALLOC(*g2,n,ndec,"converse_sg");
    g2->nv = n;
    SG_VDE(g2,v2,d2,e2);

    m = SETWORDSNEEDED(n);
#if !MAXN
    DYNALLOC1(set,workset,workset_sz,m,"putorbits");
#endif

    DYNFREE(g2->w,g2->wlen);

    ndec = 0;

    for (i = 0; i < n; ++i)
    {
	EMPTYSET(workset,m);
	for (j = v1[i]; j < v1[i]+d1[i]; ++j) ADDELEMENT(workset,e1[j]);
	if (loops == 0) ADDELEMENT(workset,i);

	v2[i] = ndec;
	for (l = 0; l < n; ++l)
	    if (!ISELEMENT(workset,l)) e2[ndec++] = l;
	d2[i] = ndec - v2[i];
    }
    g2->nde = ndec;
}

/*****************************************************************************
*                                                                            *
*  mathon_sg(g1,g2) performs a Mathon doubling operation on g1,              *
*  leaving the result in g2. g2 must exist and be initialised.               *
*                                                                            *
*****************************************************************************/

void
mathon_sg(sparsegraph *g1, sparsegraph *g2)
{
    int *e1,*d1,*e2,*d2;
    size_t *v1,*v2,j;
    int i,k,m,n1,n2;

    CHECK_SWG(g1,"mathon_sg");

    n1 = g1->nv;
    n2 = 2*n1 + 2;
    SG_ALLOC(*g2,n2,n2*(size_t)n1,"mathon_sg");
    g2->nv = n2;
    g2->nde = n2*(size_t)n1;
    DYNFREE(g2->w,g2->wlen);

    SG_VDE(g1,v1,d1,e1);
    SG_VDE(g2,v2,d2,e2);

    m = SETWORDSNEEDED(n1);
#if !MAXN
    DYNALLOC1(set,workset,workset_sz,m,"mathon_sg");
#endif

    for (i = 0; i < n2; ++i)
    {
	v2[i] = i*(size_t)n1;
	d2[i] = 0;
    }

    for (i = 0; i < n1; ++i)
    {
        e2[v2[0]+(d2[0]++)] = i+1;
        e2[v2[i+1]+(d2[i+1]++)] = 0;
	e2[v2[n1+1]+(d2[n1+1]++)] = i+n1+2;
	e2[v2[i+n1+2]+(d2[i+n1+2]++)] = n1+1;
    }

    for (i = 0; i < n1; ++i)
    {
	EMPTYSET(workset,m);
	for (j = v1[i]; j < v1[i]+d1[i]; ++j)
	{
	    k = e1[j];
	    if (k == i) continue;   /* ignore loops */
	    ADDELEMENT(workset,k);
	    e2[v2[i+1]+(d2[i+1]++)] = k+1;
	    e2[v2[i+n1+2]+(d2[i+n1+2]++)] = k+n1+2;
	}
	for (k = 0; k < n1; ++k)
	    if (k != i && !ISELEMENT(workset,k))
	    {
		e2[v2[i+1]+(d2[i+1]++)] = k+n1+2;
                e2[v2[k+n1+2]+(d2[k+n1+2]++)] = i+1;
	    }
    }
}

/*****************************************************************************
*                                                                            *
*  mathon(g1,m1,n1,g2,m2,n2) performs a Mathon doubling operation on g1,     *
*  leaving the result in g2.                                                 *
*  m1,n1 and m2,n2 are the values of m,n before and after the operation.     *
*                                                                            *
*****************************************************************************/

void
mathon(graph *g1, int m1, int n1, graph *g2, int m2, int n2)
{
    int i,j,ii,jj;
    long li;
    set *rowptr,*gp;

    for (li = (long)m2 * (long)n2; --li >= 0;) g2[li] = 0;

    for (i = 1; i <= n1; ++i)
    {
        ii = i + n1 + 1;
        gp = GRAPHROW(g2,0,m2);        /* unnecessarily convoluted code */
        ADDELEMENT(gp,i);             /* needed to avoid compiler bug  */
        gp = GRAPHROW(g2,i,m2);        /* in MSDOS version */
        ADDELEMENT(gp,0);
        gp = GRAPHROW(g2,n1+1,m2);
        ADDELEMENT(gp,ii);
        gp = GRAPHROW(g2,ii,m2);
        ADDELEMENT(gp,n1+1);
    }

    for (i = 0, rowptr = g1; i < n1; ++i, rowptr += m1)
        for (j = 0; j < n1; ++j)
            if (j != i)
            {
                ii = i + n1 + 2;
                jj = j + n1 + 2;
                if (ISELEMENT(rowptr,j))
                {
                    gp = GRAPHROW(g2,i+1,m2);
                    ADDELEMENT(gp,j+1);
                    gp = GRAPHROW(g2,ii,m2);
                    ADDELEMENT(gp,jj);
                }
                else
                {
                    gp = GRAPHROW(g2,i+1,m2);
                    ADDELEMENT(gp,jj);
                    gp = GRAPHROW(g2,ii,m2);
                    ADDELEMENT(gp,j+1);
                }
            }
}

/*****************************************************************************
*                                                                            *
*  rangraph(g,digraph,invprob,m,n) makes a random graph (or digraph if       *
*  digraph!=FALSE) with edge probability 1/invprob.                          *
*                                                                            *
*****************************************************************************/

void
rangraph(graph *g, boolean digraph, int invprob, int m, int n)
{
    int i,j;
    long li;
    set *row,*col;

    for (li = (long)m * (long)n; --li >= 0;) g[li] = 0;

    for (i = 0, row = g; i < n; ++i, row += m)
        if (digraph)
        {
            for (j = 0; j < n; ++j)
                if (KRAN(invprob) == 0) ADDELEMENT(row,j);
        }
        else
        {
            for (j = i + 1, col = GRAPHROW(g,j,m); j < n; ++j, col += m)
                if (KRAN(invprob) == 0)
                {
                    ADDELEMENT(row,j);
                    ADDELEMENT(col,i);
                }
        }
}


/*****************************************************************************
*                                                                            *
*  rangraph2(g,digraph,p1,p2,m,n) makes a random graph (or digraph if        *
*  digraph!=FALSE) with edge probability p1/p2.                              *
*                                                                            *
*****************************************************************************/

void
rangraph2(graph *g, boolean digraph, int p1, int p2, int m, int n)
{
    int i,j;
    long li;
    set *row,*col;

    for (li = (long)m * (long)n; --li >= 0;) g[li] = 0;

    for (i = 0, row = g; i < n; ++i, row += m)
        if (digraph)
        {
            for (j = 0; j < n; ++j)
                if (KRAN(p2) < p1) ADDELEMENT(row,j);
        }
        else
            for (j = i + 1, col = GRAPHROW(g,j,m); j < n; ++j, col += m)
                if (KRAN(p2) < p1)
                {
                    ADDELEMENT(row,j);
                    ADDELEMENT(col,i);
                }
}

/*****************************************************************************
*                                                                            *
*  rangraph2_sg(sg,digraph,p1,p2,n) makes a random graph (or digraph if    *
*  digraph!=FALSE) with edge probability p1/p2. sg must be initialised.      *
*                                                                            *
*****************************************************************************/

void
rangraph2_sg(sparsegraph *sg, boolean digraph, int p1, int p2, int n)
{
    int i,j,k;
    int *dd,*ee;
    double rn,expec,var,sd;
    int ldeg;
    size_t *vv,inc,nde;
 
    sg->nv = n;

    rn = n;
    expec = (rn*rn-rn)*(double)p1/(double)p2;
    var = expec*(double)(p2-p1)/(double)p2;
    if (!digraph) var *= 2.0;
    sd = 1.0;
    if (var > 1.0)
        for (i = 0; i < 19; ++i) sd = (sd + var/sd) / 2.0;
    inc = sd + 20;
 
    SG_ALLOC(*sg,n,(size_t)expec+4*inc,"rangraph2_sg");
    SG_VDE(sg,vv,dd,ee);
    DYNFREE(sg->w,sg->wlen);

    for (i = 0; i < n; ++i) dd[i] = 0;
    vv[0] = 0;
    nde = 0;

    if (!digraph)
    {
	for (i = 0; i < n; ++i)
	{
	    ldeg = 0;
	    for (j = i+1; j < n; ++j)
	    if (KRAN(p2) < p1)
	    {
		nde += 2;
		if (nde > sg->elen)
		{
		    DYNREALLOC(int,sg->e,sg->elen,sg->elen+inc,
						"rangraph2_sg realloc");
		    ee = sg->e;
		}
		ee[vv[i]+ldeg++] = j;
		++dd[j];
	    }
	    if (i < n-1) vv[i+1] = vv[i] + dd[i] + ldeg;
	    dd[i] = ldeg;
	}
        for (i = 0; i < n; ++i)
	    for (k = 0; k < dd[i]; ++k)
	    {
		j = ee[vv[i]+k];
		if (j > i) ee[vv[j]+dd[j]++] = i;
	    }
	sg->nde = nde;	
    }
    else
    {
	for (i = 0; i < n; ++i)
        {
            ldeg = 0;
            for (j = 0; j < n; ++j)
            if (j != i && KRAN(p2) < p1)
            {
                ++nde;
                if (nde > sg->elen)
                {
                    DYNREALLOC(int,sg->e,sg->elen,sg->elen+inc,
                                             "rangraph2_sg realloc");
                    ee = sg->e;
                }
                ee[vv[i]+ldeg++] = j;
            }
	    if (i < n-1) vv[i+1] = vv[i] + ldeg;
            dd[i] = ldeg;
        }
	sg->nde = nde;	
    }
}

/****************************************************************************/

static void
putsequence(FILE *f, int *x, int linelength, int n)
/* Write n integers to f with equal values collapsed.
 * labelorg is used. */
{
    char s[60];
    int j,v1,v2,xval,curlen;

    curlen = 0;
    v1 = 0;
    while (v1 < n)
    {
        xval = x[v1];

        for (v2 = v1; v2 < n - 1 && x[v2+1] == xval; ++v2) {}
        j = itos(v1+labelorg,s);
        if (v2 > v1)
        {
            s[j++] = '-';
            j += itos(v2+labelorg,&s[j]);
        }
        s[j++] = ':';
        j += itos(xval,&s[j]);
        s[j] = ' ';
        s[j+1] = '\0';
        if (linelength > 0 && curlen + j >= linelength)
        {
            PUTC('\n',f);
            curlen = 0;
        }
        curlen += j + 1;
        putstring(f,s);
        v1 = v2 + 1;
    }
    PUTC('\n',f);
}

/****************************************************************************/

static void
putnumbers(FILE *f, int *x, int linelength, int n)
/* Write n integers to f with equal values combined as multiplicities.
 * labelorg is NOT used. */
{
    char s[60];
    int j,v1,v2,xval,curlen;

    curlen = 0;
    v1 = 0;
    while (v1 < n)
    {
        xval = x[v1];

        for (v2 = v1; v2 < n - 1 && x[v2+1] == xval; ++v2) {}
	if (v2 > v1)
	{
	    j = itos(v2-v1+1,s);
	    s[j++] = '*';
	}
	else
	    j = 0;

        j += itos(xval,&s[j]);
        s[j] = ' ';
        s[j+1] = '\0';
        if (linelength > 0 && curlen + j >= linelength)
        {
            PUTC('\n',f);
            curlen = 0;
        }
        curlen += j + 1;
        putstring(f,s);
        v1 = v2 + 1;
    }
    PUTC('\n',f);
}

/*****************************************************************************
*                                                                            *
*  putdegs(f,g,linelength,m,n)  writes the degree of each vertex of g to     *
*  file f, using at most linelength characters per line (excluding '\n').    *
*  A value of linelength <= 0 dictates no line breaks at all.                *
*  labelorg is used.                                                         *
*                                                                            *
*  FUNCTIONS CALLED : itos(),putstring(),setsize()                           *
*                                                                            *
*****************************************************************************/

void
putdegs(FILE *f, graph *g, int linelength, int m, int n)
{
    int i;
    graph *gp;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"putdegs");
#endif

    for (i = 0, gp = g; i < n; ++i, gp += M)
        workperm[i] = setsize(gp,m);

    putsequence(f,workperm,linelength,n);
}

/*****************************************************************************
*                                                                            *
*  putdegseq(f,g,linelength,m,n)  writes the sorted degree sequence of g     *
*  file f, using at most linelength characters per line (excluding '\n').    *
*  A value of linelength <= 0 dictates no line breaks at all.                *
*                                                                            *
*****************************************************************************/

void
putdegseq(FILE *f, graph *g, int linelength, int m, int n)
{
    int i;
    graph *gp;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n,"putdegs");
#endif

    for (i = 0, gp = g; i < n; ++i, gp += M)
        workperm[i] = setsize(gp,m);

    sort1int(workperm,n);
    putnumbers(f,workperm,linelength,n);
}

/*****************************************************************************
*                                                                            *
*  putdegs_sg(f,sg,linelength)  writes the degree of each vertex of sg to    *
*  file f, using at most linelength characters per line (excluding '\n').    *
*  A value of linelength <= 0 dictates no line breaks at all.                *
*  labelorg is used.                                                         *
*                                                                            *
*  FUNCTIONS CALLED : itos(),putstring(),                                    *
*                                                                            *
*****************************************************************************/

void
putdegs_sg(FILE *f, sparsegraph *sg, int linelength)
{
    putsequence(f,sg->d,linelength,sg->nv);
}

/*****************************************************************************
*                                                                            *
*  putdegseq_sg(f,sg,linelength) writes the sorted degree sequence of sg to  *
*  file f, using at most linelength characters per line (excluding '\n').    *
*  A value of linelength <= 0 dictates no line breaks at all.                *
*                                                                            *
*****************************************************************************/

void
putdegseq_sg(FILE *f, sparsegraph *sg, int linelength)
{
    int i;
#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,sg->nv,"putdegs");
#endif

    for (i = 0; i < sg->nv; ++i)
        workperm[i] = sg->d[i];

    sort1int(workperm,sg->nv);
    putnumbers(f,workperm,linelength,sg->nv);
}

/*****************************************************************************
*                                                                            *
*  complement(g,m,n) replaces the graph g by its complement                  *
*  No loops are created unless there are loops present, in which case the    *
*  loops are also complemented.                                              *
*                                                                            *
*****************************************************************************/

void
complement(graph *g, int m, int n)
{
    boolean loops;
    int i,j;
    graph *gp;

#if !MAXN
    DYNALLOC1(set,workset,workset_sz,m,"complement");
#endif

    loops = FALSE;
    for (i = 0, gp = g; i < n && !loops; ++i, gp += M)
        if (ISELEMENT(gp,i)) loops = TRUE;

    EMPTYSET(workset,m);
    for (i = 0; i < n; ++ i) ADDELEMENT(workset,i);

    for (i = 0, gp = g; i < n; ++i, gp += M)
    {
        for (j = 0; j < M; ++j) gp[j] = workset[j] & ~gp[j];
        if (!loops) DELELEMENT(gp,i);
    }
}

/*****************************************************************************
*                                                                            *
*  converse(g,m,n) replaces the digraph g by its converse.                   *
*  There is no effect on an undirected graph.                                *
*                                                                            *
*****************************************************************************/

void
converse(graph *g, int m, int n)
{
    int i,j;
    graph *gi,*gj;

    for (i = 0, gi = g; i < n; ++i, gi += M)
        for (j = i+1, gj = gi+M; j < n; ++j, gj += M)
            if ((ISELEMENT(gi,j)!=0) + (ISELEMENT(gj,i)!=0) == 1)
            {
                FLIPELEMENT(gi,j);
                FLIPELEMENT(gj,i);
            }
}

/*****************************************************************************
*                                                                            *
*  naututil_check() checks that this file is compiled compatibly with the    *
*  given parameters.   If not, call exit(1).                                 *
*                                                                            *
*****************************************************************************/

void
naututil_check(int wordsize, int m, int n, int version)
{
    if (wordsize != WORDSIZE)
    {
        fprintf(ERRFILE,"Error: WORDSIZE mismatch in naututil.c\n");
        exit(1);
    }

#if MAXN
    if (m > MAXM)
    {
        fprintf(ERRFILE,"Error: MAXM inadequate in naututil.c\n");
        exit(1);
    }

    if (n > MAXN)
    {
        fprintf(ERRFILE,"Error: MAXN inadequate in naututil.c\n");
        exit(1);
    }
#endif

    if (version < NAUTYREQUIRED)
    {
        fprintf(ERRFILE,"Error: naututil.c version mismatch\n");
        exit(1);
    }
}

/*****************************************************************************
*                                                                            *
*  naututil_freedyn() - free the dynamic memory in this module               *
*                                                                            *
*****************************************************************************/

void
naututil_freedyn(void)
{
    echunk *ec1,*ec2;

#if !MAXN
    DYNFREE(workperm,workperm_sz);
    DYNFREE(workset,workset_sz);
#endif
    ec1 = first_echunk.next;

    while (ec1)
    {
        ec2 = ec1->next;
        FREES(ec1);
        ec1 = ec2;
    }
}
