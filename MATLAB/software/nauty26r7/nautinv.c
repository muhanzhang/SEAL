/*****************************************************************************
*                                                                            *
*  Vertex-invariants source file for nauty 2.6.                              *
*                                                                            *
*   Copyright (1989-2013) Brendan McKay.  All rights reserved.               *
*   Subject to waivers and disclaimers in nauty.h.                           *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       13-Mar-90 : initial release for version 1.5                          *
*       10-Nov-90 : changes for version 1.6 :                                *
*                 - added dummy routine nautinv_null()                       *
*       27-Aug-92 : renamed to version 1.7, no changes to this file          *
*        5-Jun-93 : renamed to version 1.7+, no changes to this file         *
*       18-Aug-93 : renamed to version 1.8, no changes to this file          *
*       17-Sep-93 : changes for version 1.9 :                                *
*                 - added invariant routine adjacencies()                    *
*       20-Jun-96 : changes for version 2.0 :                                *
*                 - added invariants cellfano() and cellfano2()              *
*       11-Jul-96 - added dynamic allocation                                 *
*       21-Oct-98 - use shortish in place of short for BIGNAUTY              *
*        9-Jan-00 - added nautinv_check()                                    *
*       12-Feb-00 - minor code formatting                                    *
*       16-Nov-00 - made changes listed in nauty.h                           *
*       22-Apr-01 : changes for version 2.1 :                                *
*                 - made all large dynamic memory external to routines       *
*                 - added nautinv_freedyn() to free all such memory          *
*                 - include nautinv.h rather than naututil.h                 *
*                 - removed nautinv_null()                                   *
*                 - added code to facilitate compilation into Magma          *
*                 - removed EXTDEFS                                          *
*       12-Jul-01 - use invararg in distances()                              *
*                 - fixed comments in ind and cliq routines                  *
*       21-Nov-01 : use NAUTYREQUIRED in nautinv_check()                     *
*       10-Dec-06 : remove BIGNAUTY                                          *
*       10-Nov-09 : remove types shortish and permutation                    *
*       23-Nov-09 : add refinvar()                                           *
*       12-Jun-10 : fixed identical errors in cellcliq() and cellind()       *
*       15-Jan-12 : add TLS_ATTR attributes                                  *
*       23-Aug-12 : fix getbigcells(), thanks to Fatih Demirkale             *
*       23-Jan-13 : add some parens to satisfy icc                           *
*                                                                            *
*****************************************************************************/

#define ONE_WORD_SETS
#include "nautinv.h"

#if  MAXM==1
#define M 1
#else
#define M m
#endif

#define MASH(l,i) ((((l) ^ 056345) + (i)) & 077777)
    /* : expression whose long value depends only on long l and int/long i.
	 Anything goes, preferably non-commutative. */

#define CLEANUP(l) ((int)((l) % 077777))
    /* : expression whose value depends on long l and is less than 077777
	 when converted to int then short.  Anything goes. */

#define ACCUM(x,y)   x = (((x) + (y)) & 077777)
    /* : must be commutative. */

static const int fuzz1[] = {037541,061532,005257,026416};
static const int fuzz2[] = {006532,070236,035523,062437};

#define FUZZ1(x) ((x) ^ fuzz1[(x)&3])
#define FUZZ2(x) ((x) ^ fuzz2[(x)&3])

#define MAXCLIQUE 10    /* max clique size for cliques() and maxindset() */

#if MAXN
static TLS_ATTR int workshort[MAXN+2];
static TLS_ATTR int vv[MAXN],ww[MAXN];
static TLS_ATTR int workperm[MAXN];
static TLS_ATTR int bucket[MAXN+2];
static TLS_ATTR int count[MAXN];
static TLS_ATTR set workset[MAXM];
static TLS_ATTR set w01[MAXM],w02[MAXM],w03[MAXM],w12[MAXM],w13[MAXM],w23[MAXM];
static TLS_ATTR set pt0[MAXM],pt1[MAXM],pt2[MAXM];
static TLS_ATTR set wss[MAXCLIQUE-1][MAXM];
static TLS_ATTR set ws1[MAXM],ws2[MAXM];
#else
DYNALLSTAT(int,workshort,workshort_sz);
DYNALLSTAT(int,vv,vv_sz);
DYNALLSTAT(int,ww,ww_sz);
DYNALLSTAT(int,workperm,workperm_sz);
DYNALLSTAT(int,bucket,bucket_sz);
DYNALLSTAT(int,count,count_sz);
DYNALLSTAT(set,ws1,ws1_sz);
DYNALLSTAT(set,ws2,ws2_sz);
DYNALLSTAT(set,workset,workset_sz);
DYNALLSTAT(set,w01,w01_sz);
DYNALLSTAT(set,w02,w02_sz);
DYNALLSTAT(set,w03,w03_sz);
DYNALLSTAT(set,w12,w12_sz);
DYNALLSTAT(set,w13,w13_sz);
DYNALLSTAT(set,w23,w23_sz);
DYNALLSTAT(set,pt0,pt0_sz);
DYNALLSTAT(set,pt1,pt1_sz);
DYNALLSTAT(set,pt2,pt2_sz);
DYNALLSTAT(set,wss,wss_sz);
#endif

/* aproto: header new_nauty_protos.h */

/*****************************************************************************
*                                                                            *
*  This file contains a number of procedures which compute vertex-invariants *
*  for stronger partition refinement.   Since entirely different             *
*  vertex-invariants seem to work better for different types of graph, we    *
*  cannot do more than give a small collection of representative examples.   *
*  Any serious computations with difficult graphs may well need to use       *
*  specially-written vertex-invariants.  The use of vertex-invariants        *
*  procedures is supported by nauty from version 1.5 onwards, via the        *
*  options userinvarproc, mininvarlevel, maxinvarlevel and invararg.         *
*  The meaning of these fields in detail are as follows:                     *
*     userinvarproc  is the address of the vertex-invariant procedure.  If   *
*                    no vertex-invariants is required, this field should     *
*                    have the value NULL.                                    *
*     maxinvarlevel  The absolute value of this is the maximum level in the  *
*                    search tree at which the vertex-invariant will be       *
*                    computed.  The root of the tree is at level 1, so the   *
*                    vertex-invariant will not be invoked at all if          *
*                    maxinvarlevel==0.  Negative values of maxinvarlevel     *
*                    request nauty to not compute the vertex-invariant at    *
*                    a level greater than that of the earliest node (if any) *
*                    on the path to the first leaf of the search tree at     *
*                    which the vertex-invariant refines the partition.       *
*     mininvarlevel  The absolute value of this is the minimum level in the  *
*                    search tree at which the vertex-invariant will be       *
*                    computed.  The root of the tree is at level 1, so there *
*                    is no effective limit if mininvarlevel is -1, 0 or 1.   *
*                    Negative values of mininvarlevel request nauty to not   *
*                    compute the vertex-invariant at a level less than       *
*                    that of the earliest node (if any) on the path to the   *
*                    first leaf of the search tree at which the              *
*                    vertex-invariant refines the partition.                 *
*     invararg       is passed to the vertex-invariant procedure via the     *
*                    argument of the same name.  It can be used by the       *
*                    procedure for any purpose.                              *
*  Note that negative values of maxinvarlevel and mininvarlevel make the     *
*  canonical labelling invalid, but can speed automorphism group finding.    *
*  Nauty already knows this and takes their absolute values.                 *
*                                                                            *
*  A vertex-invariant must be declared thus:                                 *
*  void invarproc(g,lab,ptn,level,numcells,tvpos,invar,invararg,digraph,m,n) *
*  All of these arguments must be treated as read-only except for invar.     *
*  g        : the graph, exactly as passed to nauty()                        *
*  lab,ptn  : the current partition nest (see nauty.h for the format)        *
*  level    : the level of this node in the search tree.                     *
*  numcells : the number of cells in the partition at this node.             *
*  tvpos    : the index in (lab,ptn) of one cell in the partition.           *
*             If level <= 1, the cell will be the first fragment of the      *
*             first active cell (as provided by the initial call to nauty),  *
*             or the first cell, if there were no active cells.              *
*             If level > 1, the cell will be the singleton cell which was    *
*             created to make this node of the search tree from its parent.  *
*  invararg : a copy of options.invararg                                     *
*  digraph  : a copy of options.digraph                                      *
*  m,n      : size parameters as passed to nauty()                           *
*  invar    : an array to return the answer in.   The procedure must put in  *
*             each invar[i]  (0 <= i < n)  an invariant of the 6-tuple       *
*             (<vertex i>,g,<the partition nest to this level>,level,        *
*               invararg,digraph)                                            *
*             Note that invar[] is declared as an int array.  Since the      *
*             absolute value of the invariant is irrelevant, only the        *
*             comparative values, any short, int or long value can be        *
*             assigned to the entries of invar[] without fear.  However,     *  
*             you should assign a value less than 077777 to ensure machine-  *
*             independence of the canonical labelling.                       *
*                                                                            *
*  The refinement procedure has already been called before the invariant     *
*  procedure is called.  That means that the partition is equitable if       *
*  digraph==FALSE.                                                           *
*                                                                            *
*****************************************************************************/

/*****************************************************************************
*                                                                            *
*  twopaths() assigns to each vertex v the sum of the weights of each vertex *
*  which can be reached from v along a walk of length two (including itself  *
*  usually).  The weight of each vertex w is defined as the ordinal number   *
*  of the cell containing w, starting at 1 for the first cell.               *
*                                                                            *
*****************************************************************************/

void
twopaths(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
         int *invar, int invararg, boolean digraph, int m, int n)
{
        int i,v,w;
        int wt;
        set *gv,*gw;

#if !MAXN
	DYNALLOC1(set,workset,workset_sz,m,"twopaths");
	DYNALLOC1(int,workshort,workshort_sz,n+2,"twopaths");
#endif

        wt = 1;
        for (i = 0; i < n; ++i)
        {
            workshort[lab[i]] = wt;
            if (ptn[i] <= level) ++wt;
        }

        for (v = 0, gv = (set*)g; v < n; ++v, gv += M)
        {
            EMPTYSET(workset,m);
            w = -1;
            while ((w = nextelement(gv,M,w)) >= 0)
            {
                gw = GRAPHROW(g,w,m);
                for (i = M; --i >= 0;) UNION(workset[i],gw[i]);
            }
            wt = 0;
            w = -1;
            while ((w = nextelement(workset,M,w)) >= 0) ACCUM(wt,workshort[w]);
            invar[v] = wt;
        }
}

/*****************************************************************************
*                                                                            *
*  quadruples() assigns to each vertex v a value depending on the set of     *
*  weights w(v,v1,v2,v3), where w(v,v1,v2,v3) depends on the number of       *
*  vertices adjacent to an odd number of {v,v1,v2,v3}, and to the cells      *
*  that v,v1,v2,v3 belong to.  {v,v1,v2,v3} are permitted to range over all  *
*  distinct 4-tuples which contain at least one member in the cell tvpos.    *
*                                                                            *
*****************************************************************************/

void
quadruples(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
           int *invar, int invararg, boolean digraph, int m, int n)
{
        int i,pc;
        setword sw;
        set *gw;
        int wt;
        int v,iv,v1,v2,v3;
        set *gv;
        long wv,wv1,wv2,wv3;

#if !MAXN
        DYNALLOC1(int,workshort,workshort_sz,n+2,"quadruples");
	DYNALLOC1(set,ws1,ws1_sz,m,"quadruples");
	DYNALLOC1(set,workset,workset_sz,m,"quadruples");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        wt = 1;
        for (i = 0; i < n; ++i)
        {
            workshort[lab[i]] = FUZZ2(wt);
            if (ptn[i] <= level) ++wt;
        }

        iv = tvpos - 1;
        do
        {
             v = lab[++iv];
             gv = GRAPHROW(g,v,m);
             wv = workshort[v];
             for (v1 = 0; v1 < n-2; ++v1)
             {
                wv1 = workshort[v1];
                if (wv1 == wv && v1 <= v) continue;
                wv1 += wv;
                gw = GRAPHROW(g,v1,m);
                for (i = M; --i >= 0;) workset[i] = gv[i] ^ gw[i];
                for (v2 = v1+1; v2 < n-1; ++v2)
                {
                    wv2 = workshort[v2];
                    if (wv2 == wv && v2 <= v) continue;
                    wv2 += wv1;
                    gw = GRAPHROW(g,v2,m);
                    for (i = M; --i >= 0;) ws1[i] = workset[i] ^ gw[i];
                    for (v3 = v2+1; v3 < n; ++v3)
                    {
                        wv3 = workshort[v3];
                        if (wv3 == wv && v3 <= v) continue;
                        wv3 += wv2;
                        gw = GRAPHROW(g,v3,m);
                        pc = 0;
                        for (i = M; --i >= 0;)
                            if ((sw = ws1[i] ^ gw[i]) != 0) pc += POPCOUNT(sw);
                        wt = (FUZZ1(pc)+wv3) & 077777;
                        wt = FUZZ2(wt);
                        ACCUM(invar[v],wt);
                        ACCUM(invar[v1],wt);
                        ACCUM(invar[v2],wt);
                        ACCUM(invar[v3],wt);
                    }
                }
            }
        }
        while (ptn[iv] > level);
}

/*****************************************************************************
*                                                                            *
*  triples() assigns to each vertex v a value depending on the set of        *
*  weights w(v,v1,v2), where w(v,v1,v2) depends on the number of vertices    *
*  adjacent to an odd number of {v,v1,v2}, and to the cells that             *
*  v,v1,v2 belong to.  {v,v1,v2} are permitted to range over all distinct    *
*  triples which contain at least one member in the cell tvpos.              *
*                                                                            *
*****************************************************************************/

void
triples(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
        int *invar, int invararg, boolean digraph, int m, int n)
{
        int i,pc;
        setword sw;
        set *gw;
        int wt;
        int v,iv,v1,v2;
        set *gv;
        long wv,wv1,wv2;

#if !MAXN
        DYNALLOC1(set,workset,workset_sz,m,"triples");
        DYNALLOC1(int,workshort,workshort_sz,n+2,"triples");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        wt = 1;
        for (i = 0; i < n; ++i)
        {
            workshort[lab[i]] = FUZZ1(wt);
            if (ptn[i] <= level) ++wt;
        }

        iv = tvpos - 1;
        do
        {
             v = lab[++iv];
             gv = GRAPHROW(g,v,m);
             wv = workshort[v];
             for (v1 = 0; v1 < n-1; ++v1)
             {
                wv1 = workshort[v1];
                if (wv1 == wv && v1 <= v) continue;
                wv1 += wv;
                gw = GRAPHROW(g,v1,m);
                for (i = M; --i >= 0;) workset[i] = gv[i] ^ gw[i];
                for (v2 = v1+1; v2 < n; ++v2)
                {
                    wv2 = workshort[v2];
                    if (wv2 == wv && v2 <= v) continue;
                    wv2 += wv1;
                    gw = GRAPHROW(g,v2,m);
                    pc = 0;
                    for (i = M; --i >= 0;)
                        if ((sw = workset[i] ^ gw[i]) != 0) pc += POPCOUNT(sw);
                    wt = (FUZZ1(pc)+wv2) & 077777;
                    wt = FUZZ2(wt);
                    ACCUM(invar[v],wt);
                    ACCUM(invar[v1],wt);
                    ACCUM(invar[v2],wt);
                }
            }
        }
        while (ptn[iv] > level);
}

/*****************************************************************************
*                                                                            *
*  adjtriang() assigns to each vertex v a value depending on the numbers     *
*  of common neighbours between each pair {v1,v2} of neighbours of v, and    *
*  which cells v1 and v2 lie in.  The vertices v1 and v2 must be adjacent    *
*  if invararg == 0 and not adjacent if invararg == 1.                       *
*                                                                            *
*****************************************************************************/

void
adjtriang(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
          int *invar, int invararg, boolean digraph, int m, int n)
{
        int j,pc;
        setword sw;
        set *gi;
        int wt;
        int i,v1,v2;
        boolean v1v2;
        set *gv1,*gv2;

#if !MAXN
        DYNALLOC1(set,workset,workset_sz,m,"adjtriang");
        DYNALLOC1(int,workshort,workshort_sz,n+2,"adjtriang");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        wt = 1;
        for (i = 0; i < n; ++i)
        {
            workshort[lab[i]] = FUZZ1(wt);
            if (ptn[i] <= level) ++wt;
        }

        for (v1 = 0, gv1 = g; v1 < n; ++v1, gv1 += M)
        {
            for (v2 = (digraph ? 0 : v1+1); v2 < n; ++v2)
            {
                if (v2 == v1) continue;
                v1v2 = (ISELEMENT(gv1,v2) != 0);
                if ((invararg == 0 && !v1v2)
			 || (invararg == 1 && v1v2)) continue;
                wt = workshort[v1];
                ACCUM(wt,workshort[v2]);
                ACCUM(wt,v1v2);

                gv2 = GRAPHROW(g,v2,m);
                for (i = M; --i >= 0;) workset[i] = gv1[i] & gv2[i];
                i = -1;
                while ((i = nextelement(workset,M,i)) >= 0)
                {
                    pc = 0;
                    gi = GRAPHROW(g,i,m);
                    for (j = M; --j >= 0;)
                        if ((sw = workset[j] & gi[j]) != 0) pc += POPCOUNT(sw);
                    pc = (pc + wt) & 077777;
                    ACCUM(invar[i],pc);
                }
            }
        }
}

/*****************************************************************************
*                                                                            *
*  getbigcells(ptn,level,minsize,bigcells,cellstart,cellsize,n) is an        *
*  auxiliary procedure to make a list of all the large cells in the current  *
*  partition.  On entry, ptn, level and n have their usual meanings,         *
*  while minsize is the smallest size of an interesting cell.  On return,    *
*  bigcells is the number of cells of size at least minsize, cellstart[0...] *
*  contains their starting positions in ptn, and cellsize[0...] contain      *
*  their sizes.  These two arrays are in increasing order of cell size,      *
*  then position.                                                            *
*                                                                            *
*****************************************************************************/

void
getbigcells(int *ptn, int level, int minsize, int *bigcells,
            int *cellstart, int *cellsize, int n)
{
        int cell1,cell2,j;
        int si,st;
        int bc,i,h;

        bc = 0;
        for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
        {
            for (cell2 = cell1; ptn[cell2] > level; ++cell2) {}

            if (cell2 >= cell1 + minsize - 1)
            {
                cellstart[bc] = cell1;
                cellsize[bc] = cell2 - cell1 + 1;
                ++bc;
            }
        }
        *bigcells = bc;

        j = bc / 3;
        h = 1;
        do
            h = 3 * h + 1;
        while (h < j);

        do                      /* shell sort */
        {
            for (i = h; i < bc; ++i)
            {
                st = cellstart[i];
                si = cellsize[i];
                for (j = i; cellsize[j-h] > si ||
                            (cellsize[j-h] == si && cellstart[j-h] > st); )
                {
                    cellsize[j] = cellsize[j-h];
                    cellstart[j] = cellstart[j-h];
                    if ((j -= h) < h) break;
                }
                cellsize[j] = si;
                cellstart[j] = st;
            }
            h /= 3;
        }
        while (h > 0);
}

/*****************************************************************************
*                                                                            *
*  celltrips() assigns to each vertex v a value depending on the set of      *
*  weights w(v,v1,v2), where w(v,v1,v2) depends on the number of vertices    *
*  adjacent to an odd number of {v,v1,v2}.  {v,v1,v2} are  constrained to    *
*  belong to the same cell.  We try the cells in increasing order of size,   *
*  and stop as soon as any cell splits.                                      *
*                                                                            *
*****************************************************************************/

void
celltrips(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
          int *invar, int invararg, boolean digraph, int m, int n)
{
        int i,pc;
        setword sw;
        set *gw;
        int wt;
        int v,iv,v1,iv1,v2,iv2;
        int icell,bigcells,cell1,cell2;
        int *cellstart,*cellsize;
        set *gv;

#if !MAXN
        DYNALLOC1(set,workset,workset_sz,m,"celltrips");
        DYNALLOC1(int,workshort,workshort_sz,n+2,"celltrips");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        cellstart = workshort;
        cellsize = workshort + (n/2);
        getbigcells(ptn,level,3,&bigcells,cellstart,cellsize,n);

        for (icell = 0; icell < bigcells; ++icell)
        {
            cell1 = cellstart[icell];
            cell2 = cell1 + cellsize[icell] - 1;
            for (iv = cell1; iv <= cell2 - 2; ++iv)
            {
                v = lab[iv];
                gv = GRAPHROW(g,v,m);
                for (iv1 = iv + 1; iv1 <= cell2 - 1; ++iv1)
                {
                    v1 = lab[iv1];
                    gw = GRAPHROW(g,v1,m);
                    for (i = M; --i >= 0;) workset[i] = gv[i] ^ gw[i];
                    for (iv2 = iv1 + 1; iv2 <= cell2; ++iv2)
                    {
                        v2 = lab[iv2];
                        gw = GRAPHROW(g,v2,m);
                        pc = 0;
                        for (i = M; --i >= 0;)
                            if ((sw = workset[i] ^ gw[i]) != 0)
                                pc += POPCOUNT(sw);
                        wt = FUZZ1(pc);
                        ACCUM(invar[v],wt);
                        ACCUM(invar[v1],wt);
                        ACCUM(invar[v2],wt);
                    }
                }
            }
            wt = invar[lab[cell1]];
            for (i = cell1 + 1; i <= cell2; ++i)
                if (invar[lab[i]] != wt) return;
        }
}

/*****************************************************************************
*                                                                            *
*  cellquads() assigns to each vertex v a value depending on the set of      *
*  weights w(v,v1,v2,v3), where w(v,v1,v2,v3) depends on the number of       *
*  vertices adjacent to an odd number of {v,v1,v2,v3}.  {v,v1,v2,v3} are     *
*  constrained to belong to the same cell.  We try the cells in increasing   *
*  order of size, and stop as soon as any cell splits.                       *
*                                                                            *
*****************************************************************************/

void
cellquads(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
          int *invar, int invararg, boolean digraph, int m, int n)
{
        int i,pc;
        setword sw;
        set *gw;
        int wt;
        int v,iv,v1,iv1,v2,iv2,v3,iv3;
        int icell,bigcells,cell1,cell2;
        int *cellstart,*cellsize;
        set *gv;

#if !MAXN
        DYNALLOC1(set,workset,workset_sz,m,"cellquads");
        DYNALLOC1(int,workshort,workshort_sz,n+2,"cellquads");
	DYNALLOC1(set,ws1,ws1_sz,m,"cellquads");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        cellstart = workshort;
        cellsize = workshort + (n/2);
        getbigcells(ptn,level,4,&bigcells,cellstart,cellsize,n);

        for (icell = 0; icell < bigcells; ++icell)
        {
            cell1 = cellstart[icell];
            cell2 = cell1 + cellsize[icell] - 1;
            for (iv = cell1; iv <= cell2 - 3; ++iv)
            {
                v = lab[iv];
                gv = GRAPHROW(g,v,m);
                for (iv1 = iv + 1; iv1 <= cell2 - 2; ++iv1)
                {
                    v1 = lab[iv1];
                    gw = GRAPHROW(g,v1,m);
                    for (i = M; --i >= 0;) workset[i] = gv[i] ^ gw[i];
                    for (iv2 = iv1 + 1; iv2 <= cell2 - 1; ++iv2)
                    {
                        v2 = lab[iv2];
                        gw = GRAPHROW(g,v2,m);
                        for (i = M; --i >= 0;) ws1[i] = workset[i] ^ gw[i];
                        for (iv3 = iv2 + 1; iv3 <= cell2; ++iv3)
                        {
                            v3 = lab[iv3];
                            gw = GRAPHROW(g,v3,m);
                            pc = 0;
                            for (i = M; --i >= 0;)
                                if ((sw = ws1[i] ^ gw[i]) != 0)
                                    pc += POPCOUNT(sw);
                            wt = FUZZ1(pc);
                            ACCUM(invar[v],wt);
                            ACCUM(invar[v1],wt);
                            ACCUM(invar[v2],wt);
                            ACCUM(invar[v3],wt);
                        }
                    }
                }
            }
            wt = invar[lab[cell1]];
            for (i = cell1 + 1; i <= cell2; ++i)
                if (invar[lab[i]] != wt) return;
        }
}

/*****************************************************************************
*                                                                            *
*  cellquins() assigns to each vertex v a value depending on the set of      *
*  weights w(v,v1,v2,v3,v4), where w(v,v1,v2,v3,v4) depends on the number    *
*  of vertices adjacent to an odd number of {v,v1,v2,v3,v4}.                 *
*  {v,v1,v2,v3,v4} are constrained to belong to the same cell.  We try the   *
*  cells in increasing order of size, and stop as soon as any cell splits.   *
*                                                                            *
*****************************************************************************/

void
cellquins(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
          int *invar, int invararg, boolean digraph, int m, int n)
{
        int i,pc;
        setword sw;
        set *gw;
        int wt;
        int v,iv,v1,iv1,v2,iv2,v3,iv3,v4,iv4;
        int icell,bigcells,cell1,cell2;
        int *cellstart,*cellsize;
        set *gv;

#if !MAXN
        DYNALLOC1(set,workset,workset_sz,m,"cellquins");
        DYNALLOC1(int,workshort,workshort_sz,n+2,"cellquins");
	DYNALLOC1(set,ws1,ws1_sz,m,"cellquins");
	DYNALLOC1(set,ws2,ws2_sz,m,"cellquins");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        cellstart = workshort;
        cellsize = workshort + (n/2);
        getbigcells(ptn,level,5,&bigcells,cellstart,cellsize,n);

        for (icell = 0; icell < bigcells; ++icell)
        {
            cell1 = cellstart[icell];
            cell2 = cell1 + cellsize[icell] - 1;
            for (iv = cell1; iv <= cell2 - 4; ++iv)
            {
                v = lab[iv];
                gv = GRAPHROW(g,v,m);
                for (iv1 = iv + 1; iv1 <= cell2 - 3; ++iv1)
                {
                    v1 = lab[iv1];
                    gw = GRAPHROW(g,v1,m);
                    for (i = M; --i >= 0;) workset[i] = gv[i] ^ gw[i];
                    for (iv2 = iv1 + 1; iv2 <= cell2 - 2; ++iv2)
                    {
                        v2 = lab[iv2];
                        gw = GRAPHROW(g,v2,m);
                        for (i = M; --i >= 0;) ws1[i] = workset[i] ^ gw[i];
                        for (iv3 = iv2 + 1; iv3 <= cell2 - 1; ++iv3)
                        {
                            v3 = lab[iv3];
                            gw = GRAPHROW(g,v3,m);
                            for (i = M; --i >= 0;) ws2[i] = ws1[i] ^ gw[i];
                            for (iv4 = iv3 + 1; iv4 <= cell2; ++iv4)
                            {
                                v4 = lab[iv4];
                                gw = GRAPHROW(g,v4,m);
                                pc = 0;
                                for (i = M; --i >= 0;)
                                    if ((sw = ws2[i] ^ gw[i]) != 0)
                                        pc += POPCOUNT(sw);
                                wt = FUZZ1(pc);
                                ACCUM(invar[v],wt);
                                ACCUM(invar[v1],wt);
                                ACCUM(invar[v2],wt);
                                ACCUM(invar[v3],wt);
                                ACCUM(invar[v4],wt);
                            }
                        }
                    }
                }
            }
            wt = invar[lab[cell1]];
            for (i = cell1 + 1; i <= cell2; ++i)
                if (invar[lab[i]] != wt) return;
        }
}

/*****************************************************************************
*                                                                            *
*  uniqinter(s1,s2,m) returns the number in both sets if it is unique,       *
*  or -1 if there is none or it is not unique.                               *
*****************************************************************************/

static int
uniqinter(set *s1, set *s2, int m)
{
	int i,j;
	setword w;

	for (i = 0; i < M; ++i)
	{
	    if ((w = s1[i] & s2[i]) != 0)
	    {
		j = FIRSTBITNZ(w);
		if (w != BITT[j]) return -1;
		j += TIMESWORDSIZE(i);
		for (++i; i < M; ++i)
		    if (s1[i] & s2[i]) return -1;
		return j;
	    }
	}
	return -1;
}

/*****************************************************************************
*                                                                            *
*  cellfano2() assigns to each vertex v a value depending on the set of      *
*  weights w(v,v1,v2,v3), where w(v,v1,v2,v3) depends on the number of       *
*  fano-plane analogues containing {v,v1,v2,v3}.  {v,v1,v2,v3} are           *
*  constrained to belong to the same cell and being independent and          *
*  non-collinear.  We try the cells in increasing order of size, and stop    *
*  as soon as any cell splits.                                               *
*                                                                            *
*****************************************************************************/

void 
cellfano2(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
          int *invar, int invararg, boolean digraph, int m, int n)
{
        int i,pc;
        setword sw;
        int wt;
        int v0,v1,v2,v3,iv0,iv1,iv2,iv3;
        int icell,bigcells,cell1,cell2;
        int *cellstart,*cellsize;
	int nw,x01,x02,x03,x12,x13,x23;
	int pnt0,pnt1,pnt2;
	set *gv0,*gv1,*gv2,*gv3;
	set *gp0,*gp1,*gp2;

#if !MAXN
        DYNALLOC1(int,workshort,workshort_sz,n+2,"cellfano2");
	DYNALLOC1(int,vv,vv_sz,n,"cellfano2");
	DYNALLOC1(int,ww,ww_sz,n,"cellfano2");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        cellstart = workshort;
        cellsize = workshort + (n/2);
        getbigcells(ptn,level,4,&bigcells,cellstart,cellsize,n);

        for (icell = 0; icell < bigcells; ++icell)
        {
            cell1 = cellstart[icell];
            cell2 = cell1 + cellsize[icell] - 1;
            for (iv0 = cell1; iv0 <= cell2 - 3; ++iv0)
            {
                v0 = lab[iv0];
                gv0 = GRAPHROW(g,v0,m);
		nw = 0;
		for (iv1 = iv0 + 1; iv1 <= cell2; ++iv1)
		{
		    v1 = lab[iv1];
                    if (ISELEMENT(gv0,v1)) continue;
                    if ((x01 = uniqinter(gv0,GRAPHROW(g,v1,m),m)) < 0) continue;
		    vv[nw] = v1;
		    ww[nw] = x01;
		    ++nw;
		}	

                for (iv1 = 0; iv1 < nw-2; ++iv1)
                {
                    v1 = vv[iv1];
                    gv1 = GRAPHROW(g,v1,m);
		    x01 = ww[iv1];

                    for (iv2 = iv1 + 1; iv2 < nw-1; ++iv2)
                    {
			x02 = ww[iv2];
			if (x02 == x01) continue;
                        v2 = vv[iv2];
			if (ISELEMENT(gv1,v2)) continue;
                        gv2 = GRAPHROW(g,v2,m);
			if ((x12 = uniqinter(gv1,gv2,m)) < 0) continue;

                        for (iv3 = iv2 + 1; iv3 < nw; ++iv3)
                        {
			    x03 = ww[iv3];
			    if (x03 == x01 || x03 == x02) continue;
                            v3 = vv[iv3];
			    if (ISELEMENT(gv1,v3) || ISELEMENT(gv2,v3))
				continue;
                            gv3 = GRAPHROW(g,v3,m);
			    if ((x13 = uniqinter(gv1,gv3,m)) < 0) continue;
			    if ((x23 = uniqinter(gv2,gv3,m)) < 0
                                                   || x23 == x13) continue;

			    if ((pnt0 = uniqinter(GRAPHROW(g,x01,m),
						 GRAPHROW(g,x23,m),m)) < 0)
				continue;
			    if ((pnt1 = uniqinter(GRAPHROW(g,x02,m),
                                                 GRAPHROW(g,x13,m),m)) < 0)
                                continue;
                            if ((pnt2 = uniqinter(GRAPHROW(g,x03,m),
                                                 GRAPHROW(g,x12,m),m)) < 0)
                                continue;

			    gp0 = GRAPHROW(g,pnt0,m);
			    gp1 = GRAPHROW(g,pnt1,m);
			    gp2 = GRAPHROW(g,pnt2,m);

			    pc = 0;
			    for (i = M; --i >= 0;)
			    {
				sw = gp0[i] & gp1[i] & gp2[i];
				if (sw) pc += POPCOUNT(sw);
			    }
			    wt = FUZZ1(pc);
			    ACCUM(invar[v0],wt);
			    ACCUM(invar[v1],wt);
                            ACCUM(invar[v2],wt);
                            ACCUM(invar[v3],wt);
                        }
                    }
                }
            }
            wt = invar[lab[cell1]];
            for (i = cell1 + 1; i <= cell2; ++i)
                if (invar[lab[i]] != wt) return;
        }
}

/*****************************************************************************
*                                                                            *
*  setnbhd(g,m,n,w,wn) is an auxiliary routine that sets wn to the union     *
*  of the neighbours of the vertices in w.                                   *
*                                                                            *
*****************************************************************************/

void
setnbhd(graph *g, int m, int n, set *w, set *wn)
{
	int i,j;
	set *gi;

	i = nextelement(w,M,-1);
	if (i < 0)
	{
	    EMPTYSET(wn,M);
	    return;
	}

	gi = GRAPHROW(g,i,M);
	for (j = M; --j >= 0;) wn[j] = gi[j];

	while ((i = nextelement(w,M,i)) >= 0)
	{
	    gi = GRAPHROW(g,i,M);
            for (j = M; --j >= 0;) wn[j] |= gi[j];
	}
}    

/*****************************************************************************
*                                                                            *
*  cellfano() assigns to each vertex v a value depending on the set of       *
*  weights w(v,v1,v2,v3), where w(v,v1,v2,v3) depends on the number of       *
*  fano-plane analogues containing {v,v1,v2,v3}.  {v,v1,v2,v3} are           *
*  constrained to belong to the same cell and being independent.  We try     *
*  the cells in increasing order of size, and stop as soon as any cell       *
*  splits.                                                                   *
*                                                                            *
*****************************************************************************/

void 
cellfano(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
         int *invar, int invararg, boolean digraph, int m, int n)
{
        int i,pc;
        setword sw;
        int wt;
        int v0,v1,v2,v3,iv0,iv1,iv2,iv3;
        int icell,bigcells,cell1,cell2;
        int *cellstart,*cellsize;
	set *gv0,*gv1,*gv2,*gv3;

#if !MAXN
        DYNALLOC1(int,workshort,workshort_sz,n+2,"cellfano");
	DYNALLOC1(set,w01,w01_sz,m,"cellfano");
	DYNALLOC1(set,w02,w02_sz,m,"cellfano");
	DYNALLOC1(set,w03,w03_sz,m,"cellfano");
	DYNALLOC1(set,w12,w12_sz,m,"cellfano");
	DYNALLOC1(set,w13,w13_sz,m,"cellfano");
	DYNALLOC1(set,w23,w23_sz,m,"cellfano");
	DYNALLOC1(set,pt0,pt0_sz,m,"cellfano");
	DYNALLOC1(set,pt1,pt1_sz,m,"cellfano");
	DYNALLOC1(set,pt2,pt2_sz,m,"cellfano");
	DYNALLOC1(set,workset,workset_sz,m,"cellfano");
#else
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        cellstart = workshort;
        cellsize = workshort + (n/2);
        getbigcells(ptn,level,4,&bigcells,cellstart,cellsize,n);

        for (icell = 0; icell < bigcells; ++icell)
        {
            cell1 = cellstart[icell];
            cell2 = cell1 + cellsize[icell] - 1;
            for (iv0 = cell1; iv0 <= cell2 - 3; ++iv0)
            {
                v0 = lab[iv0];
                gv0 = GRAPHROW(g,v0,m);
                for (iv1 = iv0 + 1; iv1 <= cell2 - 2; ++iv1)
                {
                    v1 = lab[iv1];
		    if (ISELEMENT(gv0,v1)) continue;
                    gv1 = GRAPHROW(g,v1,m);
                    for (i = M; --i >= 0;) workset[i] = gv0[i] & gv1[i];
		    setnbhd(g,m,n,workset,w01);

                    for (iv2 = iv1 + 1; iv2 <= cell2 - 1; ++iv2)
                    {
                        v2 = lab[iv2];
			if (ISELEMENT(gv0,v2) || ISELEMENT(gv1,v2))
			    continue;
                        gv2 = GRAPHROW(g,v2,m);
			for (i = M; --i >= 0;) workset[i] = gv0[i] & gv2[i];
                        setnbhd(g,m,n,workset,w02);
                        for (i = M; --i >= 0;) workset[i] = gv1[i] & gv2[i];
                        setnbhd(g,m,n,workset,w12);

                        for (iv3 = iv2 + 1; iv3 <= cell2; ++iv3)
                        {
                            v3 = lab[iv3];
			    if (ISELEMENT(gv0,v3) || ISELEMENT(gv1,v3) ||
					ISELEMENT(gv2,v3))
				continue;
                            gv3 = GRAPHROW(g,v3,m);
                            for (i = M; --i >= 0;) workset[i] = gv0[i] & gv3[i];
                            setnbhd(g,m,n,workset,w03);
                            for (i = M; --i >= 0;) workset[i] = gv1[i] & gv3[i];
                            setnbhd(g,m,n,workset,w13);
                            for (i = M; --i >= 0;) workset[i] = gv2[i] & gv3[i];
                            setnbhd(g,m,n,workset,w23);
			
			    for (i = M; --i >= 0;) workset[i] = w01[i] & w23[i];
			    setnbhd(g,m,n,workset,pt0);
                            for (i = M; --i >= 0;) workset[i] = w03[i] & w12[i];
                            setnbhd(g,m,n,workset,pt1);
                            for (i = M; --i >= 0;) workset[i] = w02[i] & w13[i];
                            setnbhd(g,m,n,workset,pt2);
			    pc = 0;
			    for (i = M; --i >= 0;)
			    {
				sw = pt0[i] & pt1[i] & pt2[i];
				if (sw) pc += POPCOUNT(sw);
			    }
			    wt = FUZZ1(pc);
			    ACCUM(invar[v0],wt);
			    ACCUM(invar[v1],wt);
                            ACCUM(invar[v2],wt);
                            ACCUM(invar[v3],wt);
                        }
                    }
                }
            }
            wt = invar[lab[cell1]];
            for (i = cell1 + 1; i <= cell2; ++i)
                if (invar[lab[i]] != wt) return;
        }
}

/*****************************************************************************
*                                                                            *
*  distances() assigns to each vertex v a value depending on the number of   *
*  vertices at each distance from v, and what cells they lie in.             *
*  If we find any cell which is split in this manner, we don't try any       *
*  further cells.                                                            *
*                                                                            *
*****************************************************************************/

void
distances(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
          int *invar, int invararg, boolean digraph, int m, int n)
{
        int i;
        set *gw;
        int wt;
        int d,dlim,cell1,cell2,iv,v,w;
        boolean success;

#if !MAXN
        DYNALLOC1(set,workset,workset_sz,m,"distances");
        DYNALLOC1(int,workshort,workshort_sz,n+2,"distances");
	DYNALLOC1(set,ws1,ws1_sz,m,"distances");
	DYNALLOC1(set,ws2,ws2_sz,m,"distances"); 
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        wt = 1;
        for (i = 0; i < n; ++i)
        {
            workshort[lab[i]] = FUZZ1(wt);
            if (ptn[i] <= level) ++wt;
        }

	if (invararg > n || invararg == 0) dlim = n;
	else                               dlim = invararg+1;

        success = FALSE;
        for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
        {
            for (cell2 = cell1; ptn[cell2] > level; ++cell2) {}
            if (cell2 == cell1) continue;

            for (iv = cell1; iv <= cell2; ++iv)
            {
                v = lab[iv];
                EMPTYSET(ws1,m);
                ADDELEMENT(ws1,v);
                EMPTYSET(ws2,m);
                ADDELEMENT(ws2,v);
                for (d = 1; d < dlim; ++d)
                {
                    EMPTYSET(workset,m);
                    wt = 0;
                    w = -1;
                    while ((w = nextelement(ws2,M,w)) >= 0)
                    {
                        gw = GRAPHROW(g,w,m);
                        ACCUM(wt,workshort[w]);
                        for (i = M; --i >= 0;) workset[i] |= gw[i];
                    }
                    if (wt == 0) break;
                    ACCUM(wt,d);
                    wt = FUZZ2(wt);
                    ACCUM(invar[v],wt);
                    for (i = M; --i >= 0;)
                    {
                        ws2[i] = workset[i] & ~ws1[i];
                        ws1[i] |= ws2[i];
                    }
                }
                if (invar[v] != invar[lab[cell1]]) success = TRUE;
            }
            if (success) break;
        }
}

/*****************************************************************************
*                                                                            *
*  indsets() assigns to each vertex v a value depending on which cells the   *
*  vertices which join v in an independent set lie in.  The size of the      *
*  independent sets which are used is the smallest of invararg and MAXCLIQUE.*
*                                                                            *
*****************************************************************************/

void
indsets(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
        int *invar, int invararg, boolean digraph, int m, int n)
{
        int i;
        int wt;
        set *gv;
        int ss,setsize;
        int v[MAXCLIQUE];
        long wv[MAXCLIQUE];
        set *s0,*s1;

#if !MAXN
        DYNALLOC1(int,workshort,workshort_sz,n+2,"indsets");
	DYNALLOC2(set,wss,wss_sz,m,MAXCLIQUE-1,"indsets");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        if (invararg <= 1 || digraph) return;

        if (invararg > MAXCLIQUE) setsize = MAXCLIQUE;
        else                      setsize = invararg;

        wt = 1;
        for (i = 0; i < n; ++i)
        {
            workshort[lab[i]] = FUZZ2(wt);
            if (ptn[i] <= level) ++wt;
        }

        for (v[0] = 0; v[0] < n; ++v[0])
        {
            wv[0] = workshort[v[0]];
            s0 = (set*)wss;
            EMPTYSET(s0,m);
            for (i = v[0]+1; i < n; ++i) ADDELEMENT(s0,i);
            gv = GRAPHROW(g,v[0],m);
            for (i = M; --i >= 0;) s0[i] &= ~gv[i];
            ss = 1;
            v[1] = v[0];
            while (ss > 0)
            {
                if (ss == setsize)
                {
                    wt = FUZZ1(wv[ss-1]);
                    for (i = ss; --i >= 0;) ACCUM(invar[v[i]],wt);
                    --ss;
                }
                else if ((v[ss] = nextelement((set*)wss+M*(ss-1),M,v[ss])) < 0)
                    --ss;
                else
                {
                    wv[ss] = wv[ss-1] + workshort[v[ss]];
                    ++ss;
                    if (ss < setsize)
                    {
                        gv = GRAPHROW(g,v[ss-1],m);
			s1 = (set*)wss + M*(ss-2);
                        for (i = M; --i >= 0;) s1[i+M] = s1[i] & ~gv[i];
                        v[ss] = v[ss-1];
                    }
                }
            }
        }
}

/*****************************************************************************
*                                                                            *
*  cliques() assigns to each vertex v a value depending on which cells the   *
*  vertices which join v in a clique lie in.  The size of the cliques used   *
*  is the smallest of invararg and MAXCLIQUE.                                *
*                                                                            *
*****************************************************************************/

void
cliques(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
        int *invar, int invararg, boolean digraph, int m, int n)
{
        int i;
        int wt;
        set *gv;
        int ss,setsize;
        int v[MAXCLIQUE];
        long wv[MAXCLIQUE];
	set *ns;

#if !MAXN
        DYNALLOC1(int,workshort,workshort_sz,n+2,"cliques");
	DYNALLOC2(set,wss,wss_sz,m,MAXCLIQUE-1,"cliques");
#else
	set wss[MAXCLIQUE-1][MAXM];
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        if (invararg <= 1 || digraph) return;

        if (invararg > MAXCLIQUE) setsize = MAXCLIQUE;
        else                      setsize = invararg;

        wt = 1;
        for (i = 0; i < n; ++i)
        {
            workshort[lab[i]] = FUZZ2(wt);
            if (ptn[i] <= level) ++wt;
        }

        for (v[0] = 0; v[0] < n; ++v[0])
        {
            wv[0] = workshort[v[0]];
            gv = GRAPHROW(g,v[0],m);
	    ns = (set*)wss;
            for (i = M; --i >= 0;) ns[i] = gv[i];
            ss = 1;
            v[1] = v[0];
            while (ss > 0)
            {
                if (ss == setsize)
                {
                    wt = FUZZ1(wv[ss-1]);
                    for (i = ss; --i >= 0;) ACCUM(invar[v[i]],wt);
                    --ss;
                }
                else if ((v[ss] = nextelement((set*)wss+M*(ss-1),M,v[ss])) < 0)
                    --ss;
                else
                {
                    wv[ss] = wv[ss-1] + workshort[v[ss]];
                    ++ss;
                    if (ss < setsize)
                    {
                        gv = GRAPHROW(g,v[ss-1],m);
			ns = (set*)wss + M*(ss-2);
                        for (i = M; --i >= 0;) ns[i+M] = ns[i] & gv[i];
                        v[ss] = v[ss-1];
                    }
                }
            }
        }
}

/*****************************************************************************
*                                                                            *
*  cellcliq() assigns to each vertex v a value depending on the number of    *
*  cliques which v lies in and which lie in the same cell as v.  The size    *
*  of clique counted is the smallest of invararg and MAXCLIQUE.  We try the  *
*  cells in increasing order of size and stop as soon as any cell splits.    *
*                                                                            *
*****************************************************************************/

void
cellcliq(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
         int *invar, int invararg, boolean digraph, int m, int n)
{
        int i;
        int wt;
        set *gv;
        int ss,setsize;
        int v[MAXCLIQUE];
        set *ns;
        int *cellstart,*cellsize;
        int iv,icell,bigcells,cell1,cell2;
        int pc;
        setword sw;

#if !MAXN
        DYNALLOC1(set,workset,workset_sz,m,"cellcliq");
        DYNALLOC1(int,workshort,workshort_sz,n+2,"cellcliq");
	DYNALLOC2(set,wss,wss_sz,m,MAXCLIQUE-1,"cellcliq");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        if (invararg <= 1 || digraph) return;

        if (invararg > MAXCLIQUE) setsize = MAXCLIQUE;
        else                      setsize = invararg;

        cellstart = workshort;
        cellsize = workshort + (n/2);
        getbigcells(ptn,level,setsize > 6 ? setsize : 6,&bigcells,
                    cellstart,cellsize,n);

        for (icell = 0; icell < bigcells; ++icell)
        {
            cell1 = cellstart[icell];
            cell2 = cell1 + cellsize[icell] - 1;

            EMPTYSET(workset,m);
            for (iv = cell1; iv <= cell2; ++iv) ADDELEMENT(workset,lab[iv]);

            for (iv = cell1; iv <= cell2; ++iv)
            {
                v[0] = lab[iv];
                gv = GRAPHROW(g,v[0],m);
		ns = (set*)wss;
                pc = 0;

                for (i = M; --i >= 0;)
                {
                    ns[i] = gv[i] & workset[i];
                    if ((sw = ns[i]) != 0) pc += POPCOUNT(sw);
                }
                if (pc <= 1 || pc >= cellsize[icell] - 2) continue;

                ss = 1;
                v[1] = v[0];
                while (ss > 0)
                {
                    if (ss == setsize)
                    {
                        for (i = ss; --i >= 0;) ++invar[v[i]];
                        --ss;
                    }
                    else if ((v[ss] 
				= nextelement((set*)wss+M*(ss-1),M,v[ss])) < 0)
                        --ss;
                    else
                    {
                        ++ss;
                        if (ss < setsize)
                        {
                            gv = GRAPHROW(g,v[ss-1],m);
			    ns = (set*)wss + M*(ss-2);
                            for (i = M; --i >= 0;) ns[i+M] = ns[i] & gv[i];
                            v[ss] = v[ss-1];
                        }
                    }
                }
            }
            wt = invar[lab[cell1]];
            for (iv = cell1 + 1; iv <= cell2; ++iv)
                if (invar[lab[iv]] != wt) return;
        }
}

/*****************************************************************************
*                                                                            *
*  cellind() assigns to each vertex v a value depending on the number of     *
*  independent sets which v lies in and which lie in the same cell as v.     *
*  The size of clique counted is the smallest of invararg and MAXCLIQUE.     *
*  We try the cells in increasing order of size and stop as soon as any      *
*  cell splits.                                                              *
*                                                                            *
*****************************************************************************/

void
cellind(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
        int *invar, int invararg, boolean digraph, int m, int n)
{
        int i;
        int wt;
        set *gv;
        int ss,setsize;
        int v[MAXCLIQUE];
        set *ns;
        int *cellstart,*cellsize;
        int iv,icell,bigcells,cell1,cell2;
        int pc;
        setword sw;

#if !MAXN
        DYNALLOC1(set,workset,workset_sz,m,"cellind");
        DYNALLOC1(int,workshort,workshort_sz,n+2,"cellind");
	DYNALLOC2(set,wss,wss_sz,m,MAXCLIQUE-1,"cellind");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        if (invararg <= 1 || digraph) return;

        if (invararg > MAXCLIQUE) setsize = MAXCLIQUE;
        else                      setsize = invararg;

        cellstart = workshort;
        cellsize = workshort + (n/2);
        getbigcells(ptn,level,setsize > 6 ? setsize : 6,&bigcells,
                    cellstart,cellsize,n);

        for (icell = 0; icell < bigcells; ++icell)
        {
            cell1 = cellstart[icell];
            cell2 = cell1 + cellsize[icell] - 1;

            EMPTYSET(workset,m);
            for (iv = cell1; iv <= cell2; ++iv) ADDELEMENT(workset,lab[iv]);

            for (iv = cell1; iv <= cell2; ++iv)
            {
                v[0] = lab[iv];
                gv = GRAPHROW(g,v[0],m);
	 	ns = (set*)wss;
                pc = 0;

                for (i = M; --i >= 0;)
                {
                    ns[i] = ~gv[i] & workset[i];
                    if ((sw = ns[i]) != 0) pc += POPCOUNT(sw);
                }
                if (pc <= 1 || pc >= cellsize[icell] - 2) continue;

                ss = 1;
                v[1] = v[0];
                while (ss > 0)
                {
                    if (ss == setsize)
                    {
                        for (i = ss; --i >= 0;) ++invar[v[i]];
                        --ss;
                    }
                    else if ((v[ss] 
			   = nextelement((set*)wss+M*(ss-1),M,v[ss])) < 0)
                        --ss;
                    else
                    {
                        ++ss;
                        if (ss < setsize)
                        {
                            gv = GRAPHROW(g,v[ss-1],m);
			    ns = (set*)wss + M*(ss-2);
                            for (i = M; --i >= 0;) ns[i+M] = ns[i] & ~gv[i];
                            v[ss] = v[ss-1];
                        }
                    }
                }
            }
            wt = invar[lab[cell1]];
            for (iv = cell1 + 1; iv <= cell2; ++iv)
                if (invar[lab[iv]] != wt) return;
        }
}

/*****************************************************************************
*                                                                            *
*  adjacencies() assigns to each vertex v a code depending on which cells    *
*  it is joined to and from, and how many times.  It is intended to provide  *
*  better partitioning that the normal refinement routine for digraphs.      *
*  It will not help with undirected graphs in nauty at all.                  *
*                                                                            *
*****************************************************************************/

void
adjacencies(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
            int *invar, int invararg, boolean digraph, int m, int n)
{
        int i,v,w;
        int vwt,wwt;
        set *gv;

#if !MAXN
        DYNALLOC1(int,workshort,workshort_sz,n+2,"adjacencies");
#endif

        vwt = 1;
        for (i = 0; i < n; ++i)
        {
            workshort[lab[i]] = vwt;
            if (ptn[i] <= level) ++vwt;
            invar[i] = 0;
        }

        for (v = 0, gv = (set*)g; v < n; ++v, gv += M)
        {
            vwt = FUZZ1(workshort[v]);
            wwt = 0;
            w = -1;
            while ((w = nextelement(gv,M,w)) >= 0)
            {
                ACCUM(wwt,FUZZ2(workshort[w]));
                ACCUM(invar[w],vwt);
            }
            ACCUM(invar[v],wwt);
        }
}

/*****************************************************************************
*                                                                            *
*  nautinv_check() checks that this file is compiled compatibly with the     *
*  given parameters.   If not, call exit(1).                                 *
*                                                                            *
*****************************************************************************/

void
nautinv_check(int wordsize, int m, int n, int version)
{
        if (wordsize != WORDSIZE)
        {
            fprintf(ERRFILE,"Error: WORDSIZE mismatch in nautinv.c\n");
            exit(1);
        }

#if MAXN
        if (m > MAXM)
        {
            fprintf(ERRFILE,"Error: MAXM inadequate in nautinv.c\n");
            exit(1);
        }

        if (n > MAXN)
        {
            fprintf(ERRFILE,"Error: MAXN inadequate in nautinv.c\n");
            exit(1);
        }
#endif

        if (version < NAUTYREQUIRED)
        {
            fprintf(ERRFILE,"Error: nautinv.c version mismatch\n");
            exit(1);
        }
}

/*****************************************************************************
*                                                                            *
*  nautinv_freedyn() - free the dynamic memory in this module                *
*                                                                            *
*****************************************************************************/
 
void
nautinv_freedyn(void)
{
#if !MAXN
	DYNFREE(workset,workset_sz);
	DYNFREE(workshort,workshort_sz);
	DYNFREE(ws1,ws1_sz);
	DYNFREE(ws2,ws2_sz);
	DYNFREE(vv,vv_sz);
	DYNFREE(ww,ww_sz);
	DYNFREE(w01,w01_sz);
	DYNFREE(w02,w02_sz);
	DYNFREE(w03,w03_sz);
	DYNFREE(w12,w12_sz);
	DYNFREE(w13,w13_sz);
	DYNFREE(w23,w23_sz);
	DYNFREE(pt0,pt0_sz);
	DYNFREE(pt1,pt1_sz);
	DYNFREE(pt2,pt2_sz);
	DYNFREE(wss,wss_sz);
#endif
}

/*===================================================================*/


/*****************************************************************************
*                                                                            *
*  semirefine(g,lab,ptn,level,numcells,strength,active,m,n) performs a       *
*  refinement operation on the partition at the specified level of the       *
*  partition nest (lab,ptn).  *numcells is assumed to contain the number of  *
*  cells on input, and is updated.  The initial set of active cells (alpha   *
*  in the paper) is specified in the set active.  Precisely, x is in active  *
*  iff the cell starting at index x in lab is active.                        *
*  *code is set to a value which depends on the fine detail of the           *
*  algorithm, but which is independent of the labelling of the graph.        *
*                                                                            *
*****************************************************************************/

static int
semirefine(graph *g, int *lab, int *ptn, int level, int *numcells,
           int strength, set *active, int m, int n)
{
	int i,c1,c2,labc1;
	setword x;
	set *set1,*set2;
	int split1,split2,cell1,cell2;
	int cnt,bmin,bmax;
	long longcode;
	set *gptr;
	int maxcell,maxpos,hint;

#if !MAXN
	DYNALLOC1(int,workperm,workperm_sz,n,"refine");
	DYNALLOC1(set,workset,workset_sz,m,"refine");
	DYNALLOC1(int,bucket,bucket_sz,n+2,"refine");
	DYNALLOC1(int,count,count_sz,n,"refine");
#endif

	longcode = *numcells;
	split1 = -1;
	hint = 0;
	while (*numcells < n && ((split1 = hint, ISELEMENT(active,split1))
                             || (split1 = nextelement(active,M,split1)) >= 0
	                     || (split1 = nextelement(active,M,-1)) >= 0))
	{
	    DELELEMENT(active,split1);
	    for (split2 = split1; ptn[split2] > level; ++split2) {}
	    longcode = MASH(longcode,split1+split2);
	    if (split1 == split2)       /* trivial splitting cell */
	    {
	        gptr = GRAPHROW(g,lab[split1],M);
	        for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
	        {
	            for (cell2 = cell1; ptn[cell2] > level; ++cell2) {}
	            if (cell1 == cell2) continue;
	            c1 = cell1;
	            c2 = cell2;
	            while (c1 <= c2)
	            {
	                labc1 = lab[c1];
	                if (ISELEMENT(gptr,labc1))
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
	                ptn[c2] = level;
	                longcode = MASH(longcode,FUZZ1(c2));
	                ++*numcells;
			if (ISELEMENT(active,cell1) || c2-cell1 >= cell2-c1)
			{
     			    ADDELEMENT(active,c1);
			    if (c1 == cell2) hint = c1;
			}
			else
			{
     			    ADDELEMENT(active,cell1);
			    if (c2 == cell1) hint = cell1;
			}
	            }
	        }
	    }
	    else        /* nontrivial splitting cell */
	    {
	        EMPTYSET(workset,m);
	        for (i = split1; i <= split2; ++i)
	            ADDELEMENT(workset,lab[i]);
	        longcode = MASH(longcode,FUZZ2(split2-split1+1));

	        for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
	        {
	            for (cell2 = cell1; ptn[cell2] > level; ++cell2) {}
	            if (cell1 == cell2) continue;
	            i = cell1;
	            set1 = workset;
	            set2 = GRAPHROW(g,lab[i],m);
	            cnt = 0;
	            for (c1 = m; --c1 >= 0;)
	                if ((x = ((*set1++) & (*set2++))) != 0)
	                    cnt += POPCOUNT(x);

	            count[i] = bmin = bmax = cnt;
	            bucket[cnt] = 1;
	            while (++i <= cell2)
	            {
	                set1 = workset;
	                set2 = GRAPHROW(g,lab[i],m);
	                cnt = 0;
	                for (c1 = m; --c1 >= 0;)
	                    if ((x = ((*set1++) & (*set2++))) != 0)
	                        cnt += POPCOUNT(x);

	                while (bmin > cnt) bucket[--bmin] = 0;
	                while (bmax < cnt) bucket[++bmax] = 0;
	                ++bucket[cnt];
	                count[i] = cnt;
	            }
	            if (bmin == bmax)
	            {
	                longcode = MASH(longcode,FUZZ1(bmin+cell1));
	                continue;
	            }
	            c1 = cell1;
		    maxcell = -1;
	            for (i = bmin; i <= bmax; ++i)
	                if (bucket[i])
	                {
	                    c2 = c1 + bucket[i];
	                    bucket[i] = c1;
	                    longcode = MASH(longcode,i+c1);
			    if (c2-c1 > maxcell)
			    {
				maxcell = c2-c1;
				maxpos = c1;
			    }
	                    if (c1 != cell1)
	                    {
	                        ADDELEMENT(active,c1);
			        if (c2-c1 == 1) hint = c1;
	                        ++*numcells;
	                    }
	                    if (c2 <= cell2) ptn[c2-1] = level;
	                    c1 = c2;
	                }
	            for (i = cell1; i <= cell2; ++i)
	                workperm[bucket[count[i]]++] = lab[i];
	            for (i = cell1; i <= cell2; ++i) lab[i] = workperm[i];
		    if (!ISELEMENT(active,cell1))
		    {
     			ADDELEMENT(active,cell1);
     			DELELEMENT(active,maxpos);     /* check maxpos is alwas defined */
		    }
	        }
	    }
	    if (--strength == 0) break;   /* negative is fine! */
	}

	longcode = MASH(longcode,FUZZ2(*numcells));
	return CLEANUP(longcode);
}

void 
refinvar(graph *g, int *lab, int *ptn, int level, int numcells, int tvpos,
          int *invar, int invararg, boolean digraph, int m, int n)
{
        int i,j;
        int wt;
        int icell,bigcells,cell1,cell2;
        int *cellstart,*cellsize;
	int newnumcells;

#if !MAXN
        DYNALLOC1(int,workshort,workshort_sz,n+2,"refinvar");
        DYNALLOC1(int,vv,vv_sz,n,"refinvar");
        DYNALLOC1(int,ww,ww_sz,n,"refinvar");
        DYNALLOC1(set,ws1,ws1_sz,n,"refinvar");
#endif

        for (i = n; --i >= 0;) invar[i] = 0;

        cellstart = workshort;
        cellsize = workshort + (n/2);
        getbigcells(ptn,level,2,&bigcells,cellstart,cellsize,n);

        for (icell = 0; icell < bigcells; ++icell)
        {
            cell1 = cellstart[icell];
            cell2 = cell1 + cellsize[icell] - 1;
            for (i = cell1; i <= cell2; ++i)
            {
		for (j = 0; j < n; ++j)
		{
		    vv[j] = lab[j];
		    ww[j] = ptn[j];
		}
		newnumcells = numcells + 1;
		ww[cell1] = level;
		EMPTYSET(ws1,m);
		ADDELEMENT(ws1,cell1);
		vv[i] = lab[cell1];
		vv[cell1] = lab[i];
		invar[lab[i]] = semirefine(g,vv,ww,level,&newnumcells,
						invararg,ws1,m,n);
            }
            wt = invar[lab[cell1]];
            for (i = cell1 + 1; i <= cell2; ++i)
                if (invar[lab[i]] != wt) return;
        }
}
