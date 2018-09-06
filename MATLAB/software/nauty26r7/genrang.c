/* genrang.c  version 2.1; B D McKay, Feb 15, 2016 */
/* TODO:  Check allocs for no edges */

#define USAGE \
"genrang [-P#|-P#/#|-e#|-r#|-R#|-d#] [-l#] [-m#] [-t] [-T] [-a] \n" \
"         [-s|-g|-z] [-S#] [-q] n|n1,n2 num [outfile]"

#define HELPTEXT \
" Generate random graphs.\n\
     n  : number of vertices\n\
     n1,n2 : number of vertices (bipartite graph)\n\
    num : number of graphs\n\
\n\
    A bipartite variant is only available if specified below.\n\
\n\
    -s  : Write in sparse6 format (default)\n\
    -g  : Write in graph6 format\n\
    -z  : Make random digraphs and write in digraph6 format\n\
    -P#/# : Give edge probability; -P# means -P1/#.\n\
          Bipartite version available.\n\
    -e# : Give the number of edges\n\
          Bipartite version available.\n\
    -r# : Make regular of specified degree\n\
    -d# : Make regular of specified degree (pseudorandom)\n\
          Bipartite version: this is the degree on the first side\n\
    -R# : Make regular of specified degree but output\n\
          as vertex count, edge count, then list of edges\n\
    -l# : Maximum loop multiplicity (default 0)\n\
    -m# : Maximum multiplicity of non-loop edge (default and minimum 1)\n\
    -t  : Make a random tree\n\
    -T  : Make a random tournament (implies -z)\n\
    -a  : Make invariant under a random permutation\n\
    -S# : Specify random generator seed (default nondeterministic)\n\
    -q  : suppress auxiliary output\n"

#define MAXLREG 88   /* Max value for -r or -R switch (multigraphs) */
/* This is also the limit for -r if MAXN > 0. */
/* The limit for simple graphs is in naututil-h.in */

/*************************************************************************

   Oct 27, 2004 : corrected handling of -P values
**************************************************************************/

#include "gtools.h"

static long seed;

/*************************************************************************/

static void
perminvar(graph *g, int *perm, int m, int n)
/* Add to g the least number of edges needed to make perm
   an automorphism. */
{
    int i,j,ii,jj;
    set *gi,*gpi;

    for (i = 0, gi = (set*)g; i < n; ++i, gi += m)
    {
        gpi = g + m * 1L * perm[i];
        for (j = -1; (j = nextelement(gi,m,j)) >= 0; )
            if (!ISELEMENT(gpi,perm[j]))
            {
                ii = perm[i];
                jj = perm[j];
                while (ii != i || jj != j)
                {
                    ADDELEMENT(g+m*1L*ii,jj);
                    ii = perm[ii];
                    jj = perm[jj];
                }
            }
    }
}

/**************************************************************************/

static void
gcomplement(graph *g, boolean loopstoo, int m, int n)
/* Replace g by its complement */
{
    int i,j;
    graph *gp;
#if MAXN
    set mask[MAXM];
#else
    DYNALLSTAT(set,mask,mask_sz);
    DYNALLOC1(set,mask,mask_sz,m,"complement");
#endif

    EMPTYSET(mask,m);
    for (i = 0; i < n; ++i) ADDELEMENT(mask,i);

    if (loopstoo)
    {
        for (i = 0, gp = g; i < n; ++i, gp += m)
        {
            for (j = 0; j < m; ++j) gp[j] ^= mask[j];
        }
    }
    else
    {
        for (i = 0, gp = g; i < n; ++i, gp += m)
        {
	    DELELEMENT(mask,i);
            for (j = 0; j < m; ++j) gp[j] ^= mask[j];
	    ADDELEMENT(mask,i);
        }
    }
}

/**************************************************************************/

static void
gcomplement_bip(graph *g, int m, int n1, int n2)
/* Replace g by its bipartite complement */
{
    int i,j,n;
    graph *gp;
#if MAXN
    set mask[MAXM];
#else
    DYNALLSTAT(set,mask,mask_sz);
    DYNALLOC1(set,mask,mask_sz,m,"complement");
#endif

    n = n1 + n2;

    EMPTYSET(mask,m);
    for (i = n1; i < n; ++i) ADDELEMENT(mask,i);

    for (i = 0, gp = g; i < n1; ++i, gp += m)
        for (j = 0; j < m; ++j) gp[j] ^= mask[j];

    EMPTYSET(mask,m);
    for (i = 0; i < n1; ++i) ADDELEMENT(mask,i);

    for (i = n1, gp = GRAPHROW(g,n1,m); i < n; ++i, gp += m)
        for (j = 0; j < m; ++j) gp[j] ^= mask[j];
}

/**************************************************************************/

static void
ranedges(long e, boolean loopsok, graph *g, int m, int n)
/* Random graph with n vertices and e edges */
{
    unsigned long ln,li,nc2,ned,sofar;
    set *gi,*gj;
    int i,j;

    ln = n;
    nc2 = (ln&1) ? ln*((ln-1)/2) : (ln/2)*(ln-1);
    if (loopsok) nc2 += ln;

    if (e + e > nc2) ned = nc2 - e;
    else             ned = e;
    sofar = 0;

    for (li = m*ln; --li != 0;) g[li] = 0; g[0] = 0;

    while (sofar < ned)
    {
        i = KRAN(n);
        if (loopsok) j = KRAN(n);
        else do { j = KRAN(n); } while (i == j);
        gi = GRAPHROW(g,i,m);
        if (!ISELEMENT(gi,j))
        {
            ADDELEMENT(gi,j);
            gj = GRAPHROW(g,j,m);
            ADDELEMENT(gj,i);
            ++sofar;
        }
    }

    if (ned != e) gcomplement(g,loopsok,m,n);
}

/**************************************************************************/

static void
ranedges_bip(long e, graph *g, int m, int n1, int n2)
/* Random bipartite graph with n1+n2 vertices and e edges */
{
    size_t ln,li,nc2,ned,sofar;
    set *gi,*gj;
    int i,j,n;

    n = n1 + n2;
    nc2 = (size_t)n1*n2;

    if (e + e > nc2) ned = nc2 - e;
    else             ned = e;
    sofar = 0;

    for (li = m*(size_t)n; --li != 0;) g[li] = 0; g[0] = 0;

    while (sofar < ned)
    {
        i = KRAN(n1);
        j = n1 + KRAN(n2);
        gi = GRAPHROW(g,i,m);
        if (!ISELEMENT(gi,j))
        {
            ADDELEMENT(gi,j);
            gj = GRAPHROW(g,j,m);
            ADDELEMENT(gj,i);
            ++sofar;
        }
    }

    if (ned != e) gcomplement_bip(g,m,n1,n2);
}

/**************************************************************************/

static void
grandtourn(graph *g, int m, int n)
/* Random tournament */
{
    int i,j;
    long li;
    set *row,*col;

    for (li = (long)m * (long)n; --li >= 0;) g[li] = 0;

    for (i = 0, row = g; i < n; ++i, row += m)
    {
        for (j = i+1, col = GRAPHROW(g,i+1,m); j < n; ++j, col += m)
            if (KRAN(2) < 1) ADDELEMENT(row,j);
	    else             ADDELEMENT(col,i);
    }
}

/**************************************************************************/

static void
grandtourn_bip(graph *g, int m, int n1, int n2)
/* Random bipartite tournament */
{
    int i,j,n;
    long li;
    set *row,*col;

    n = n1 + n2;
    for (li = (long)m * (long)n; --li >= 0;) g[li] = 0;

    for (i = 0, row = g; i < n1; ++i, row += m)
    {
        for (j = n1, col = GRAPHROW(g,n1,m); j < n; ++j, col += m)
            if (KRAN(2) < 1) ADDELEMENT(row,j);
	    else             ADDELEMENT(col,i);
    }
}

/**************************************************************************/

static void
grandgraph(graph *g, boolean digraph, boolean loopsok,
                                         int p1, int p2, int m, int n)
/* Random graph with probability p1/p2 */
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
	    if (!loopsok) DELELEMENT(row,i);
        }
        else
        {
            for (j = i + (!loopsok), col = GRAPHROW(g,j,m);
                                            j < n; ++j, col += m)
                if (KRAN(p2) < p1)
                {
                    ADDELEMENT(row,j);
                    ADDELEMENT(col,i);
                }
        }
}

/**************************************************************************/

static void
grandgraph_bip(graph *g, boolean digraph,
                         int p1, int p2, int m, int n1, int n2)
/* Random bipartite graph with probability p1/p2 */
{
    int i,j,n;
    long li;
    set *row,*col;

    n = n1 + n2;

    for (li = (long)m * (long)n; --li >= 0;) g[li] = 0;

    for (i = 0, row = g; i < n1; ++i, row += m)
        if (digraph)
        {
            for (j = n1; j < n; ++j)
                if (KRAN(p2) < p1) ADDELEMENT(row,j);
        }
        else
        {
            for (j = n1, col = GRAPHROW(g,j,m);
                                            j < n; ++j, col += m)
                if (KRAN(p2) < p1)
                {
                    ADDELEMENT(row,j);
                    ADDELEMENT(col,i);
                }
        }
}

/**************************************************************************/

static void
ranarcs(long e, boolean loopsok, graph *g, int m, int n)
/* Random digraph graph with n vertices and e edges */
{
    unsigned long ln,li,nn,ned,sofar;
    set *gi;
    int i,j;

    ln = n;
    nn = (loopsok ? n*n : n*(n-1));

    if (e + e > nn) ned = nn - e;
    else            ned = e;
    sofar = 0;

    for (li = m*ln; --li != 0;) g[li] = 0;
    g[0] = 0;

    while (sofar < ned)
    {
        i = KRAN(n);
	if (loopsok) j = KRAN(n);
        else do { j = KRAN(n); } while (i == j);
        gi = GRAPHROW(g,i,m);
        if (!ISELEMENT(gi,j))
        {
            ADDELEMENT(gi,j);
            ++sofar;
        }
    }

    if (ned != e) gcomplement(g,loopsok,m,n);
}

/**************************************************************************/

static void
makeranreg(int *cub, int degree, int multmax, int loopmax, int n)
/* Make a random regular graph in cub[].  Each consecutive degree
   entries of cub[] is set to the neighbours of one vertex.
   The length of cub had better be at least degree*n  */
{
    long i,j,k,v,w,nn,mult;
    boolean ok;
#if MAXN
    int deg[MAXN],p[MAXLREG*MAXN];
#else
    DYNALLSTAT(int,deg,deg_sz);
    DYNALLSTAT(int,p,p_sz);
    DYNALLSTAT(int,loops,loops_sz);

    DYNALLOC1(int,deg,deg_sz,n,"genrang");
    DYNALLOC2(int,p,p_sz,degree,n,"genrang");
    DYNALLOC1(int,loops,loops_sz,n,"genrang");
#endif

    nn = n;

    for (i = j = 0; i < nn; ++i)
        for (k = 0; k < degree; ++k)
           p[j++] = i;

    do
    {
        ok = TRUE;

        for (j = degree*nn-1; j >= 1; j -= 2)
        {
            i = KRAN(j);
            k = p[j-1];
            p[j-1] = p[i];
            p[i] = k;
        }
        for (i = 0; i < nn; ++i) deg[i] = loops[i] = 0;

        for (j = degree*nn-1; j >= 1;)
        {
            v = p[j--];
            w = p[j--];
            if (v == w && ++loops[v] > loopmax)
            {
                ok = FALSE;
                break;
            }
            if (v != w && multmax < degree)
            {
                mult = 0;
                for (i = deg[w]; --i >= 0;)
                    if (cub[degree*w+i] == v && ++mult >= multmax) break;
                if (i >= 0)
                {
                    ok = FALSE;
                    break;
                }
            }
            cub[degree*w+deg[w]++] = v;
            cub[degree*v+deg[v]++] = w;
        }
    }
    while (!ok);
}

/**************************************************************************/

static void
rundmodel(int *cub, int degree, int n)
/* Make a random-ish regular graph in cub[] using the d-model.
   Each consecutive degree entries of cub[] is set to the neighbours
   of one vertex.  The length of cub had better be at least degree*n  */
{
    long iters,fails;
    size_t i,j,navail;
    int *cubi,*cubj,vi,vj,k;
    boolean ok;
#if MAXN
    int deg[MAXN],avail[MAXN*MAXLREG];
#else
    DYNALLSTAT(int,deg,deg_sz);
    DYNALLSTAT(int,avail,avail_sz);

    DYNALLOC1(int,deg,deg_sz,n,"genrang");
    DYNALLOC2(int,avail,avail_sz,n,degree,"genrang");
#endif

    iters = 0;
    do
    {
	ok = TRUE;
	++iters;

	k = 0;
        for (i = 0; i < n; ++i)
        {
	    deg[i] = 0;
	    for (j = 0; j < degree; ++j) avail[k++] = i;
        }
        navail = n*degree;

	while (navail >= 2 && ok)
	{
	    for (fails = 100 + navail; --fails >= 0;)
	    {
		i = KRAN(navail);
		do { j = KRAN(navail); } while (j == i);
		vi = avail[i];
		vj = avail[j];
                if (vi == vj) continue;
		cubi = cub + vi*degree;
		cubj = cub + vj*degree;
		for (k = deg[vi]; --k >= 0; ) if (cubi[k] == vj) break;
		if (k < 0) break;
	    }

	    if (fails >= 0)
	    {
		cubi[deg[vi]++] = vj;
		cubj[deg[vj]++] = vi;

		avail[i] = avail[navail-1];
		--navail;
		if (avail[i] == vj) j = i;

		avail[j] = avail[navail-1];
		--navail;
	    }
	    else
		ok = FALSE;
	}
	if (navail > 0) ok = FALSE;
    } while (!ok);

    /* fprintf(stderr,">C %ld iters\n",iters); */
}

/**************************************************************************/

static void
ranregR(FILE *f, int degree, int multmax, int loopmax, int n)
/* Make a random regular graph of order n and degree d and write
   it in f, as number of vertices, number of edges, list of edges */
{
    long i,j,k,l;
    int loops;
#if MAXN
    int cub[MAXLREG*MAXN];
#else
    DYNALLSTAT(int,cub,cub_sz);
    DYNALLOC2(int,cub,cub_sz,degree,n,"genrang");
#endif

    makeranreg(cub,degree,multmax,loopmax,n);

    fprintf(f,"%d %ld\n",n,n*(long)degree/2);
    l = j = 0;
    for (i = 0; i < n; ++i)
    {
        loops = 0;
        for (k = 0; k < degree; ++k, ++j)
            if (i < cub[j] || (i == cub[j] && (++loops & 1) == 0))
            {
                if (l > 0 && l % 5 == 0) fprintf(f,"\n");
                fprintf(f," %ld %d",i,cub[j]);
                ++l;
            }
    }
    fprintf(f,"\n");
    if (ferror(f)) gt_abort(">E genrang output error\n");
}
 
/**************************************************************************/

static void
ranreg(int degree, graph *g, int m, int n)
/* Make a random simple regular graph of order n and degree d and return
   it in g. */
{
    int i,j,k;
    set *gi;
#if MAXN
    int cub[MAXLREG*MAXN];
#else
    DYNALLSTAT(int,cub,cub_sz);
    DYNALLOC1(int,cub,cub_sz,degree*n,"genrang");
#endif

    makeranreg(cub,degree,1,0,n);

    j = 0;
    for (i = 0, gi = (set*)g; i < n; ++i, gi += m)
    {
        EMPTYSET(gi,m);
        for (k = 0; k < degree; ++k)
        {
            ADDELEMENT(gi,cub[j]);
            j++;
        }
    }
}

/**************************************************************************/

static void
ranreglm_sg(int degree, sparsegraph *sg, int multmax, int loopmax, int n)
/* Make a sparse random regular graph of order n and degree d
 * and return it in sg. */
{
    int i,j,k,deg,loops;
    long nde,k0;
#if MAXN
    int cub[MAXLREG*MAXN];
#else
    DYNALLSTAT(int,cub,cub_sz);
    DYNALLOC1(int,cub,cub_sz,degree*n,"genrang");
#endif

    makeranreg(cub,degree,multmax,loopmax,n);

    SG_ALLOC(*sg,n,degree*n,"genrang");

    sg->nv = n;
    j = nde = 0;
    for (i = 0; i < n; ++i)
    {
        sg->v[i] = k0 = i*degree;
        loops = deg = 0;
        for (k = 0; k < degree; ++k, ++j)
        {
            if (cub[j] == i)
            {
                /* Loops are in cub twice but sg once */
                ++loops;
                if ((loops&1)) sg->e[k0+deg++] = i;
            }
            else
                sg->e[k0+deg++] = cub[j];
        }
        sg->d[i] = deg;
        nde += deg;
    }
    sg->nde = nde;
}

/**************************************************************************/

static void
dmodel_sg(int degree, sparsegraph *sg, int n)
/* Make a sparse random-ish regular graph of order n and degree d
 * and return it in sg. */
{
    int i,j,k,deg,comdeg;
    long k0,nde;
#if MAXN
    int cub[MAXLREG*MAXN];
    boolean adj[MAXN];
#else
    DYNALLSTAT(int,cub,cub_sz);
    DYNALLSTAT(boolean,adj,adj_sz);
#endif

    SG_ALLOC(*sg,n,degree*n,"dmodel_sg");

    if (degree <= n-degree-1)
    {
#if !MAXN
        DYNALLOC1(int,cub,cub_sz,degree*n,"dmodel_sg");
#endif
        rundmodel(cub,degree,n);

        sg->nv = n;
        j = nde = 0;
        for (i = 0; i < n; ++i)
        {
            sg->v[i] = k0 = i*(long)degree;
            deg = 0;
            for (k = 0; k < degree; ++k, ++j)
                sg->e[k0+deg++] = cub[j];
            sg->d[i] = deg;
            nde += deg;
        }
        sg->nde = nde;
    }
    else
    {
	comdeg = n - degree - 1;
#if !MAXN
        DYNALLOC1(int,cub,cub_sz,comdeg*n,"dmodel_sg");
        DYNALLOC1(boolean,adj,adj_sz,n,"dmodel_sg");
#endif
        rundmodel(cub,comdeg,n);

        sg->nv = n;
        j = nde = 0;
        for (i = 0; i < n; ++i)
        {
            sg->v[i] = k0 = i*(long)degree;
            deg = 0;
	    for (k = 0; k < n; ++k) adj[k] = TRUE;
	    adj[i] = FALSE;
	    for (k = 0; k < comdeg; ++k, ++j) adj[cub[j]] = FALSE;
            for (k = 0; k < n; ++k)
		if (adj[k]) sg->e[k0+deg++] = k;
            sg->d[i] = deg;
            nde += deg;
        }
        sg->nde = nde;
    }
}

/**************************************************************************/

static void
rundmodel_bip(int *cub, int deg1, int deg2, int n1, int n2)
/* Make a random-ish semiregular bipartite graph in cub[] using the d-model.
   Each consecutive deg1/deg2 entries of cub[] is set to the neighbours
   of one vertex.  The length of cub had better be deg1*n1  */
{
    long iters,fails;
    size_t i,j,k,navail,ne;
    int *cubi,*cubj,vi,vj,n;
    boolean ok;
#if MAXN
    int deg[MAXN],avail[MAXN*MAXLREG];
#else
    DYNALLSTAT(int,deg,deg_sz);
    DYNALLSTAT(int,avail,avail_sz);

    n = n1 + n2;
    DYNALLOC1(int,deg,deg_sz,n,"genrang");
    DYNALLOC2(int,avail,avail_sz,n,2*deg1,"genrang");
#endif

    ne = n1*(size_t)deg1;

    iters = 0;
    do
    {
	ok = TRUE;
	++iters;

	k = 0;
        for (i = 0; i < n1; ++i)
        {
	    deg[i] = 0;
	    for (j = 0; j < deg1; ++j) avail[k++] = i;
        }
        for (i = n1; i < n; ++i)
        {
	    deg[i] = 0;
	    for (j = 0; j < deg2; ++j) avail[k++] = i;
        }
        navail = ne;

	while (navail >= 1 && ok)
	{
	    for (fails = 100 + navail; --fails >= 0;)
	    {
		i = KRAN(navail);
		j = ne + KRAN(navail);
		vi = avail[i];
		vj = avail[j];
		cubi = cub + vi*deg1;
		cubj = cub + ne + (vj-n1)*deg2;
		for (k = 0; k < deg[vi]; ++k) if (cubi[k] == vj) break;
		if (k == deg[vi]) break;
	    }

	    if (fails >= 0)
	    {
		cubi[deg[vi]++] = vj;
		cubj[deg[vj]++] = vi;

		avail[i] = avail[navail-1];
		avail[j] = avail[ne+navail-1];
		--navail;
	    }
	    else
		ok = FALSE;
	}
	if (navail > 0) ok = FALSE;
    } while (!ok);

    /* fprintf(stderr,">C %ld iters\n",iters); */
}

/**************************************************************************/

static void
dmodel_bip_sg(int deg1, sparsegraph *sg, int n1, int n2)
/* Make a sparse random-ish semiregular bipartite graph of order n1+n2
   and degree deg1 on the left and return it in sg. */
{
    int i,k,deg,comdeg1,comdeg2,n,deg2;
    size_t j,k0,nde,ne,comne;
#if MAXN
    int cub[MAXLREG*MAXN];
    boolean adj[MAXN];
#else
    DYNALLSTAT(int,cub,cub_sz);
    DYNALLSTAT(boolean,adj,adj_sz);
#endif

    n = n1 + n2;
    ne = n1*(size_t)deg1;
    deg2 = ne / n2;
    if (deg2*(size_t)n2 != ne || deg1 > n2 || deg2 > n1)
	gt_abort(">E genrang: impossible bipartite degrees\n");

    SG_ALLOC(*sg,n,2*ne,"dmodel_bip_sg");

    if (deg1 <= n2-deg1)
    {
#if !MAXN
        DYNALLOC1(int,cub,cub_sz,2*ne,"dmodel_bip_sg");
#endif
        rundmodel_bip(cub,deg1,deg2,n1,n2);

        sg->nv = n;
        j = nde = 0;
        for (i = 0; i < n1; ++i)
        {
            sg->v[i] = k0 = i*(size_t)deg1;
            deg = 0;
            for (k = 0; k < deg1; ++k, ++j)
                sg->e[k0+deg++] = cub[j];
            sg->d[i] = deg;
            nde += deg;
        }
        for (i = n1; i < n; ++i)
        {
            sg->v[i] = k0 = ne + (i-n1)*(size_t)deg2;
            deg = 0;
            for (k = 0; k < deg2; ++k, ++j)
                sg->e[k0+deg++] = cub[j];
            sg->d[i] = deg;
            nde += deg;
        }
        sg->nde = nde;
    }
    else
    {
	comdeg1 = n2 - deg1;
	comdeg2 = n1 - deg2;
	comne = n1*(size_t)comdeg1;
#if !MAXN
        DYNALLOC1(int,cub,cub_sz,2*comne,"dmodel_bip_sg");
        DYNALLOC1(boolean,adj,adj_sz,n,"dmodel_bip_sg");
#endif
        rundmodel_bip(cub,comdeg1,comdeg2,n1,n2);

        sg->nv = n;
        j = nde = 0;
        for (i = 0; i < n1; ++i)
        {
            sg->v[i] = k0 = i*(long)deg1;
            deg = 0;
	    for (k = n1; k < n; ++k) adj[k] = TRUE;
	    for (k = 0; k < comdeg1; ++k, ++j) adj[cub[j]] = FALSE;
            for (k = n1; k < n; ++k)
		if (adj[k]) sg->e[k0+deg++] = k;
            sg->d[i] = deg;
            nde += deg;
        }
        for (i = n1; i < n; ++i)
        {
            sg->v[i] = k0 = ne + (i-n1)*(long)deg2;
            deg = 0;
	    for (k = 0; k < n1; ++k) adj[k] = TRUE;
	    for (k = 0; k < comdeg2; ++k, ++j) adj[cub[j]] = FALSE;
            for (k = 0; k < n1; ++k)
		if (adj[k]) sg->e[k0+deg++] = k;
            sg->d[i] = deg;
            nde += deg;
        }
        sg->nde = nde;
    }
}

/**************************************************************************/

static void
randomtree(sparsegraph *sg, int n)
/* Make a random tree with n vertices */
{
    int i,v0,v1,ne,k;
#if MAXN
    int ed[2*MAXN];
#else
    DYNALLSTAT(int,ed,ed_sz);
    DYNALLOC1(int,ed,ed_sz,2*n,"randomtree");
#endif

    SG_ALLOC(*sg,n,2*(n-1),"randomtree");
    sg->nv = n;
    sg->nde = 2*(n-1);
    sg->w = NULL;

    for (i = 0; i < n; ++i) sg->d[i] = 0;

    v0 = KRAN(n);
    ne = k = 0;
    while (ne < n-1)
    {
	do { v1 = KRAN(n); } while (v1 == v0);
	if (sg->d[v1] == 0)
	{
	    ed[k++] = v0;
	    ed[k++] = v1;
	    ++ne;
	    ++sg->d[v0];
	    ++sg->d[v1];
	}
	v0 = v1;
    }

    sg->v[0] = 0;
    for (i = 1; i < n; ++i) sg->v[i] = sg->v[i-1] + sg->d[i-1];

    for (i = 0; i < n; ++i) sg->d[i] = 0;

    for (k = 0; k < 2*(n-1); )
    {
        v0 = ed[k++];
        v1 = ed[k++];
	sg->e[sg->v[v0]+(sg->d[v0])++] = v1;
	sg->e[sg->v[v1]+(sg->d[v1])++] = v0;
    }
}

/**************************************************************************/
/**************************************************************************/

#define NOBIP if (bipartite) { fprintf(stderr, \
 ">E genrang: This option is not available for bipartite graphs." \
 " Feel free to request it.\n"); exit(1); }

int
main(int argc, char *argv[])
{
    int m,n,n1,n2,codetype;
    int argnum,j;
    char *arg,sw;
    boolean badargs;
    boolean gswitch,sswitch,qswitch,Sswitch,Rswitch,lswitch,tswitch;
    boolean aswitch,P1switch,P2switch,eswitch,rswitch,mswitch,dswitch;
    boolean Tswitch;
    long numgraphs,nout,P1value,P2value,evalue,rvalue;
    nauty_counter ln,nc2;
    int Svalue,loopmax,multmax;
    static FILE *outfile;
    char *outfilename;
    sparsegraph sg;
    boolean usesparse,digraph,bipartite;

#if MAXN
    graph g[MAXM*1L*MAXN];
    int perm[MAXN];
#else
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,perm,perm_sz);
#endif

    HELP; PUTVERSION;

    gswitch = sswitch = qswitch = Sswitch = Rswitch = FALSE;
    aswitch = P1switch = P2switch = eswitch = rswitch = FALSE;
    digraph = dswitch = tswitch = lswitch = mswitch = FALSE;
    Tswitch = FALSE;
    outfilename = NULL;

    argnum = 0;
    badargs = FALSE;
    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                     SWBOOLEAN('g',gswitch)
                else SWBOOLEAN('s',sswitch)
                else SWBOOLEAN('z',digraph)
                else SWBOOLEAN('a',aswitch)
                else SWBOOLEAN('t',tswitch)
                else SWBOOLEAN('T',Tswitch)
                else SWBOOLEAN('q',qswitch)
                else SWLONG('P',P1switch,P1value,"genrang -P")
                else SWLONG('/',P2switch,P2value,"genrang -P")
                else SWLONG('e',eswitch,evalue,"genrang -e")
                else SWLONG('d',dswitch,rvalue,"genrang -d")
                else SWLONG('r',rswitch,rvalue,"genrang -r")
                else SWLONG('R',Rswitch,rvalue,"genrang -R")
                else SWINT('S',Sswitch,Svalue,"genrang -S")
                else SWINT('l',lswitch,loopmax,"genrang -l")
                else SWINT('m',mswitch,multmax,"genrang -m")
                else badargs = TRUE;
            }
        }
        else
        {
            ++argnum;
            if      (argnum == 1)
            {
		if (sscanf(arg,"%d,%d",&n1,&n2) == 2)
		{
		    bipartite = TRUE;
		    if (n1 < 1 || n2 < 1) badargs = TRUE;
		    n = n1 + n2;
		}
		else if (sscanf(arg,"%d",&n) == 1)
		{
		    bipartite = FALSE;
		    if (n < 1) badargs = TRUE;
		}
		else
		    badargs = TRUE;
            }
            else if (argnum == 2)
            {
                if (sscanf(arg,"%ld",&numgraphs) != 1 || numgraphs < 1)
                    badargs = TRUE;
            }
            else if (argnum == 3) outfilename = arg;
            else                  badargs = TRUE;
        }
    }

    if (Tswitch) digraph = TRUE;

    if ((gswitch!=0) + (sswitch!=0) + (digraph!=0) > 1)
        gt_abort(">E genrang: -gsz are incompatible\n");

    if (gswitch)      codetype = GRAPH6;
    else if (digraph) codetype = DIGRAPH6;
    else              codetype = SPARSE6;

    if (P1switch && !P2switch)
    {
        P2value = P1value;
        P1value = 1;
    }
    else if (P2switch && !P1switch)
    {
        P1value = 1;
        P1switch = TRUE;
    }

    if (P1switch && (P1value < 0 || P2value <= 0 || P1value > P2value))
        gt_abort(">E genrang: bad value for -P switch\n");

    if ((P1switch!=0) + (eswitch!=0) + (rswitch!=0) + (dswitch!=0) 
           + (Tswitch!=0) + (Rswitch!=0) + (tswitch!=0) > 1)
        gt_abort(">E genrang: -PerRdtT are incompatible\n");

    if ((sswitch!=0) + (gswitch!=0) + (Rswitch!=0) > 1)  /* REVISE */
        gt_abort(">E genrang: -sgR are incompatible\n");

    if ((aswitch!=0) + (Rswitch!=0) > 1)
        gt_abort(">E genrang: -aR are incompatible\n");

    if (!lswitch) loopmax = 0;
    if (!mswitch) multmax = 1;
    if (digraph &&
           (Rswitch || multmax>1 || loopmax>1 || dswitch || tswitch))
        gt_abort(">E genrang: -z is only compatible with -T, -P, -e and -l1\n");

    if (!digraph && (loopmax>1 || multmax>1)
                 && !(Rswitch || (rswitch && !gswitch)))
        gt_abort(">E genrang: -l>1,-m>1 need -R or -r without -g\n");

    if (multmax < 1 || loopmax < 0)
        gt_abort(">E genrang: bad value for -l or -m\n");

    if (argnum < 2 || argnum > 3) badargs = TRUE;

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!Sswitch)
    {
#ifdef INITSEED
        INITSEED;
        ran_init(seed);
#endif
    }
    else
        ran_init(Svalue);

    if (!outfilename || outfilename[0] == '-')
    {
        outfilename = "stdout";
        outfile = stdout;
    }
    else if ((outfile = fopen(outfilename,"w")) == NULL)
    {
        fprintf(stderr,"Can't open output file %s\n",outfilename);
        gt_abort(NULL);
    }

    m = (n + WORDSIZE + 1) / WORDSIZE;
    usesparse = tswitch || dswitch ||
                 (rswitch && !aswitch && codetype==SPARSE6);
#if !MAXN
    if (!Rswitch && !usesparse)
    {
        DYNALLOC2(graph,g,g_sz,n,m,"genrang");
        if (aswitch) DYNALLOC1(int,perm,perm_sz,n,"genrang");
    }
#endif

    rswitch = rswitch || Rswitch || dswitch;

#if MAXN
    if (rswitch && rvalue > MAXLREG)
    {
        fprintf(stderr,
                ">E -r/-R is only implemented for degree <= %d\n",MAXLREG);
        exit(1);
    }
#endif

    ln = n;
    if (bipartite)
	nc2 = (unsigned long)n1*n2;
    else
        nc2 = ln*loopmax + (1+(digraph!=0))*ln*(ln-1)/2*multmax;

    if (eswitch && evalue > nc2)
    {   
        fprintf(stderr,
             ">E There are no graphs of order %d and %ld edges\n",
             n,evalue);
        exit(1);
    }

    if (rswitch && (((n&1) != 0 && (rvalue&1) != 0)
        || rvalue > (n-1)*multmax+2*loopmax))
    {
        fprintf(stderr,     
             ">E There are no such graphs of order %d and degree %ld\n",
             n,rvalue);
        exit(1);
    }

    if (!P1switch)
    {
        P1value = 1;
        P2value = 2;
    }

    SG_INIT(sg);

    for (nout = 1; nout <= numgraphs; ++nout)
    {
        if (eswitch)
        {
	    if (digraph)
            {
		NOBIP;
		ranarcs(evalue,loopmax>0,g,m,n);
	    }
	    else 
	    { 
		if (bipartite)
	            ranedges_bip(evalue,g,m,n1,n2);
		else
	            ranedges(evalue,loopmax>0,g,m,n);
	    }
	}
	else if (dswitch)
        {
	    if (bipartite)
		dmodel_bip_sg(rvalue,&sg,n1,n2);
	    else
		dmodel_sg(rvalue,&sg,n);
	}
        else if (Rswitch)
	{
	    NOBIP;
	    ranregR(outfile,rvalue,multmax,loopmax,n);
	}
        else if (rswitch && usesparse)
	{
	    NOBIP;
            ranreglm_sg(rvalue,&sg,multmax,loopmax,n);
	}
        else if (rswitch && !usesparse)
	{
	    NOBIP;
	    ranreg(rvalue,g,m,n);
	}
	else if (tswitch)
	{
	    NOBIP;
	    randomtree(&sg,n);
	}
	else if (Tswitch)
	{
	    if (bipartite)
		grandtourn_bip(g,m,n1,n2);
	    else
		grandtourn(g,m,n);
	}
        else
        {
	    if (bipartite)
	        grandgraph_bip(g,digraph,P1value,P2value,m,n1,n2);
	    else
	        grandgraph(g,digraph,loopmax>0,P1value,P2value,m,n);
	}

        if (Rswitch) continue;

        if (aswitch && !usesparse)
        {
            ranperm(perm,n);
            perminvar(g,perm,m,n);
        }
        if (codetype == SPARSE6)
	{
            if (usesparse)
            {
                sortlists_sg(&sg);
                writes6_sg(outfile,&sg);
            }
            else
                writes6(outfile,g,m,n);
	}
        else if (codetype == DIGRAPH6)
	{
	    if (usesparse)
		writed6_sg(outfile,&sg);
            else
                writed6(outfile,g,m,n);
	}
        else 
	{
	    if (usesparse)
		writeg6_sg(outfile,&sg);
            else
                writeg6(outfile,g,m,n);
	}
    }

    exit(0);
}
