/*****************************************************************************
*                                                                            *
* Auxiliary procedures for use with nauty 2.5.                               *
* None of these procedures are needed by nauty or by dreadnaut.              *
*                                                                            *
*   Copyright (1984-2013) Brendan McKay.  All rights reserved.               *
*   Subject to waivers and disclaimers in nauty.h.                           *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       26-Apr-89 : initial creation for Version 1.5.                        *
*       14-Oct-90 : renamed to version 1.6 (no changes to this file)         *
*        5-Jun-93 : renamed to version 1.7+ (no changes to this file)        *
*       18-Aug-93 : renamed to version 1.8 (no changes to this file)         *
*       17-Sep-93 : renamed to version 1.9 (no changes to this file)         *
*       24-Jan-00 : renamed to version 2.0 (no changes to this file)         *
*       16-Nov-00 : made changes listed in nauty.h                           *
*        8-Aug-02 : updated for version 2.2 (dynamic storage)                *
*        3-Nov-04 : fixed names of nautaux_freedyn() and nautaux_check()     *
*       10-Dec-06 : removed BIGNAUTY                                         *
*       10-Nov-09 : removed types permutation and shortish                   *
*       15-Jan-12 : add TLS_ATTR attributes                                  *
*                                                                            *
*****************************************************************************/

#define ONE_WORD_SETS
#include "naututil.h"          /* which includes "nauty.h" and <stdio.h> */
#include "nautaux.h"

#if  MAXM==1
#define M 1
#else
#define M m
#endif

#if !MAXN
DYNALLSTAT(set,workset,workset_sz);
DYNALLSTAT(int,workperm,workperm_sz);
#else
static TLS_ATTR set workset[MAXM];   /* used for scratch work */
static TLS_ATTR int workperm[MAXN+2];
#endif

/*****************************************************************************
*                                                                            *
*  ptncode(g,lab,ptn,level,m,n) returns a long integer invariant which       *
*  depends on the (assumed equitable) partition at the stated level, and     *
*  the number of edges betwen the various cells.                             *
*  Neither nauty nor dreadnaut use this.                                     *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>,setinter()                                       *
*                                                                            *
*****************************************************************************/

long
ptncode(graph *g, int *lab, int *ptn, int level, int m, int n)
{
    int i;
    long code;
    int v1,v2,nc,inter,cellend;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"testcanlab");
    DYNALLOC1(set,workset,workset_sz,m,"testcanlab");
#endif

    /* find all cells: put starts in workperm[0..n] */

    i = nc = 0;
    code = 0;

    while (i < n)
    {
        workperm[nc++] = i;
        code = ((code << 13) ^ (code >> 19)) + i;
        while (ptn[i] > level) ++i;
        ++i;
    }
    workperm[nc] = n;

    for (v2 = 0; v2 < nc; ++v2)
    {
        EMPTYSET(workset,m);
        for (i = workperm[v2], cellend = workperm[v2+1] - 1;
                                                      i <= cellend; ++i)
            ADDELEMENT(workset,lab[i]);
        for (v1 = 0; v1 < nc; ++v1)
        {
            i = workperm[v1];
            cellend = workperm[v1+1] - 1;
            inter = setinter(workset,GRAPHROW(g,lab[i],M),M);
            code = ((code << 13) ^ (code >> 19)) + inter;
        }
    }

    return code;
}

/*****************************************************************************
*                                                                            *
*  equitable(g,lab,ptn,level,m,n) checks that the partition at the given     *
*  level is equitable.  Neither nauty nor dreadnaut use this.                *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>,setinter()                                       *
*                                                                            *
*****************************************************************************/

boolean
equitable(graph *g, int *lab, int *ptn, int level, int m, int n)
{
    int i;
    int v1,v2,nc,inter,cellend;
    boolean ok;

    /* find all cells: put starts in workperm[0..n] */

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"testcanlab");
    DYNALLOC1(set,workset,workset_sz,m,"testcanlab");
#endif

    i = nc = 0;

    while (i < n)
    {
        workperm[nc++] = i;
        while (ptn[i] > level)
            ++i;
        ++i;
    }
    workperm[nc] = n;

    ok = TRUE;
    for (v2 = 0; v2 < nc && ok; ++v2)
    {
        EMPTYSET(workset,m);
        for (i = workperm[v2], cellend = workperm[v2+1] - 1;
                                                      i <= cellend; ++i)
            ADDELEMENT(workset,lab[i]);
        for (v1 = 0; v1 < nc; ++v1)
        {
            i = workperm[v1];
            cellend = workperm[v1+1] - 1;
            if (i == cellend)
                continue;
            inter = setinter(workset,GRAPHROW(g,lab[i],M),M);
            while (++i <= cellend)
                 if (setinter(workset,GRAPHROW(g,lab[i],M),M) != inter)
                     ok = FALSE;
        }
    }

    return ok;
}

/*****************************************************************************
*                                                                            *
*  component(g,v,c,m,n) determines the set of all vertices that can be       *
*  reached along directed paths starting at vertex v, including v itself.    *
*  This set is returned as c, unless c is null.  The size of the set is      *
*  returned as the function value.                                           *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>,setinter(),nextelement()                         *
*                                                                            *
*****************************************************************************/

int
component(graph *g, int v, set *cmpt, int m, int n)
{
    int i,z;
    set newverts[MAXM],*gx;
    int head,tail,x;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n+2,"testcanlab");
#endif

    EMPTYSET(workset,m);
    ADDELEMENT(workset,v);
    head = 0;
    tail = 1;
    workperm[head] = v;

    while (head < tail && tail < n)
    {
        x = workperm[head++];
        gx = GRAPHROW(g,x,m);
        for (i = m; --i >= 0;)
        {
            newverts[i] = gx[i] & ~workset[i];
            workset[i] |= gx[i];
        }
        for (z = -1; (z = nextelement(newverts,m,z)) >= 0; )
            workperm[tail++] = z;
    }

    if (cmpt != NULL)
        for (i = m; --i >= 0;) cmpt[i] = workset[i];

    return tail;
}

/*****************************************************************************
*                                                                            *
*  nautaux_check() checks that this file is compiled compatibly with the    *
*  given parameters.   If not, call exit(1).                                 *
*                                                                            *
*****************************************************************************/

void
nautaux_check(int wordsize, int m, int n, int version)
{
    if (wordsize != WORDSIZE)
    {
        fprintf(ERRFILE,"Error: WORDSIZE mismatch in nautaux.c\n");
        exit(1);
    }

#if MAXN
    if (m > MAXM)
    {
        fprintf(ERRFILE,"Error: MAXM inadequate in nautaux.c\n");
        exit(1);
    }

    if (n > MAXN)
    {
        fprintf(ERRFILE,"Error: MAXN inadequate in nautaux.c\n");
        exit(1);
    }
#endif

    if (version < NAUTYREQUIRED)
    {
        fprintf(ERRFILE,"Error: nautaux.c version mismatch\n");
        exit(1);
    }
}

/*****************************************************************************
*                                                                            *
*  nautaux_freedyn() - free the dynamic memory in this module               *
*                                                                            *
*****************************************************************************/

void
nautaux_freedyn(void)
{
#if !MAXN
    DYNFREE(workset,workset_sz);
    DYNFREE(workperm,workperm_sz);
#endif
}
