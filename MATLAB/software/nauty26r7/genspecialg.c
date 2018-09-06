/* genspecialg.c  version 1.1; B D McKay, Feb 12, 2016 */

#define USAGE "genspecialg \n\
[-s|-g|-z|-d] [-q] \
[-p#|-c#|-e#|-k#|-b#,#|-Q#|-f#|-J#,#|-P#,#|C#,#...|G#,#...|T#,#...] [outfile]"

#define HELPTEXT \
" Generate one particular graph.\n\
     #  : size parameter called n in the descriptions\n\
\n\
    -s : Write in sparse6 format (default)\n\
    -g : Write in graph6 format\n\
    -z : Make digraph versions and write in digraph6 format\n\
    -d : Write in dreadnaut format (can be used with -z)\n\
    -q : Suppress summary\n\
\n\
    If defined, the digraph version is shown in parentheses:\n\
    -p#   : path (directed path) on n vertices.\n\
    -c#   : cycle (directed cycle) on n vertices.\n\
    -e#   : empty graph (digraph with loops only) on n vertices.\n\
    -k#   : complete graph (with loops) on n vertices\n\
    -b#,# : complete bipartite graph (directed l->r) on n vertices\n\
    -f#   : flower snark on 4*# vertices\n\
    -P#,# : generalized Petersen graph; usual one is -P5,2\n\
    -Q#   : hypercube on 2^n vertices and degree n.\n\
    -J#,# : Johnson graph J(n,k), args are n and k.\n\
    -C#,#... : circulant (di)graph.\n\
    -T#,#... : theta (di)graph Theta(#,#,...), give path lengths.\n\
    -G#,#... : (directed) grid, use negative values for open directions\n"

/* Ideas: multipartite, knesser, full trees */

#include "gtools.h"

#define MAXARGS 1000  /* Maximum argument list for multi-argument parameters */
#define SWAP(x,y) {int w = x; x = y; y = w;}

static long args[MAXARGS];

static short vmark_val = 32000;
DYNALLSTAT(short,vmark,vmark_sz);
#define MARK(i) vmark[i] = vmark_val
#define UNMARK(i) vmark[i] = 0
#define ISMARKED(i) (vmark[i] == vmark_val)
#define ISNOTMARKED(i) (vmark[i] != vmark_val)
#define RESETMARKS {if (vmark_val++ >= 32000) \
    {size_t ij; for (ij=0;ij<vmark_sz;++ij) vmark[ij]=0; vmark_val=1;}}

static void
preparemarks(size_t nn)
{
    size_t oldsize;
    short *oldpos;

    oldsize = vmark_sz;
    oldpos = vmark;
    DYNALLOC1(short,vmark,vmark_sz,nn,"preparemarks");
    if (vmark_sz != oldsize || vmark != oldpos) vmark_val = 32000;
}

/**************************************************************************/

static void
writedread(FILE *f, sparsegraph *sg, boolean digraph)
/* Write in dreadnaut format */
{
    size_t *v;
    int *d,*e,n,i,j,k;

    SG_VDE(sg,v,d,e);
    n = sg->nv;

    if (digraph) fprintf(f,"n=%d $=0 dg\n",n);
    else         fprintf(f,"n=%d $=0 g\n",n);

    for (i = 0; i < n; ++i)
    {
	for (j = 0; j < d[i]; ++j)
        {
	    k = e[v[i]+j];
	    if (k >= i || digraph) fprintf(f," %d",k);
	}
	if (i == n-1) fprintf(f,".\n$$\n");
	else          fprintf(f,";\n");
    }
} 

/**************************************************************************/

static int binom[32][16];   /* Cached binomial coefficients */

static int
binomial(int n, int k)
/* Value of binomial(n,k), error if too big for int */
{
    int i,nki,ans;
    nauty_counter work;

    if (k > n/2) k = n - k;
    if (k < 0) return 0;

    if (n < 32  && binom[n][k] > 0) return binom[n][k];

    work = 1;
    for (i = 1; i <= k; ++i)
    {
        nki = n-k+i;
        work = (work/i) * nki + (work%i) * nki / i;
        if ((int)work != work) { fprintf(stderr,"Overflow\n"); exit(1); }
    }

    ans = (int)work;
    if (n < 32) binom[n][k] = ans;

    return ans;
}

/**************************************************************************/

static void
unrank(int r, int k, int *a)
/* r-th k-set in colex order (r=0,1,...) */
{
    int i,p;

    for (i = k; i > 0; --i)
    {
        p = i - 1;
        do ++p; while (binomial(p,i) <= r);
        r -= binomial(p-1,i);
        a[i-1] = p-1;
    }
}

static int
rank(int k, int *a)
/* Rank of a[0..k-1] in colex order */
{
    int i,r;

    r = 0;
    for (i = 0; i < k; ++i)
        r += binomial(a[i],i+1);

    return r;
}

/**************************************************************************/

static int
vnumber(long *dimen, int *index, int ndimen)
{
    int i,v;
    
    v = 0;
    for (i = 0; i < ndimen; ++i)
        v = v*dimen[i] + index[i];

    return v;
}

/**************************************************************************/

static void
makepath(long n, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i;
    size_t *v,k;

    if (n < 1 || n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: bad argument for -p\n");

    if (digraph) SG_ALLOC(*sg,n,n-1,"genspecialg");
    else         SG_ALLOC(*sg,n,2UL*n-2,"genspecialg");

    SG_VDE(sg,v,d,e);

    if (digraph || n == 1)
    {
	sg->nv = n;
	sg->nde = n-1;

	for (i = 0; i < n-1; ++i)
        {
	    d[i] = 1;
	    v[i] = i;
	    e[i] = i+1;
        }
	d[n-1] = 0;
	v[n-1] = 0;
    }
    else
    {
	sg->nv = n;
	sg->nde = 2*n-2;

	d[0] = 1;
	v[0] = 0;
	e[0] = 1;
	for (i = 1, k = 1; i < n-1; ++i, k += 2)
	{
	    d[i] = 2;
	    v[i] = k;
	    e[k] = i-1;
	    e[k+1] = i+1;
	}
	d[n-1] = 1;
	v[n-1] = k;
	e[k] = n-2;
    }
}

/**************************************************************************/

static void
makecycle(long n, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i;
    size_t *v,k;

    if (!digraph && (n < 1 || n == 2 || n > NAUTY_INFINITY-2))
        gt_abort(">E genspecialg: bad argument for -c\n");
    if (digraph && (n < 1 || n > NAUTY_INFINITY-2))
        gt_abort(">E genspecialg: bad argument for -zc\n");

    if (digraph) SG_ALLOC(*sg,n,n,"genspecialg");
    else         SG_ALLOC(*sg,n,2UL*n,"genspecialg");

    SG_VDE(sg,v,d,e);

    if (digraph || n == 1)
    {
	sg->nv = n;
	sg->nde = n;

	for (i = 0; i < n-1; ++i)
        {
	    d[i] = 1;
	    v[i] = i;
	    e[i] = i+1;
        }
	d[n-1] = 1;
	v[n-1] = n-1;
	e[n-1] = 0;
    }
    else
    {
	sg->nv = n;
	sg->nde = 2UL*n;

	d[0] = 2;
	v[0] = 0;
	e[0] = 1;
	e[1] = n-1;

	for (i = 1; i < n-1; ++i)
        {
	    d[i] = 2;
	    v[i] = 2UL*i;
	    e[2UL*i] = i-1;
	    e[2UL*i+1] = i+1;
        }
	d[n-1] = 2;
	v[n-1] = 2UL*n-2;
	e[2UL*n-2] = 0;
	e[2UL*n-1] = n-2;
    }
}

/**************************************************************************/

static void
makeflowersnark(long k, boolean digraph, sparsegraph *sg)
/* Flower snark on 4*k vertices, no digraph variant 
*
 * The flower snark Jn can be constructed with the following process :
 * Build n copies of the star graph on 4 vertices. Denote the
 * central vertex of each star Ai and the outer vertices Bi, Ci and
 * Di. This results in a disconnected graph on 4n vertices with 3n
 * edges (Ai-Bi, Ai-Ci and Ai-Di for 1?i?n). Construct the n-cycle
 * (B1... Bn). This adds n edges. Finally construct the 2n-cycle
 * (C1... CnD1... Dn). This adds 2n edges. By construction, the
 * Flower snark Jn is a cubic graph with 4n vertices and 6n edges.
*/

#define FSA(i) (4*(i))
#define FSB(i) (4*(i)+1)
#define FSC(i) (4*(i)+2)
#define FSD(i) (4*(i)+3)

{
    int n,*d,*e,i,j;
    size_t *v,nde;

    if (k < 3 || k > (NAUTY_INFINITY-2)/4)
        gt_abort(">E genspecialg: bad argument for -f\n");

    n = 4*k;
    nde = 12*(size_t)k;

    SG_ALLOC(*sg,n,nde,"genspecialg");

    SG_VDE(sg,v,d,e);
    sg->nv = n;
    sg->nde = nde;

    for (i = 0; i < n; ++i)
    {
	d[i] = 0;
        v[i] = 3*(size_t)i;
    }

    for (i = 0; i < k; ++i)
    {
        e[v[FSA(i)]+d[FSA(i)]++] = FSB(i);
        e[v[FSB(i)]+d[FSB(i)]++] = FSA(i);
        e[v[FSA(i)]+d[FSA(i)]++] = FSC(i);
        e[v[FSC(i)]+d[FSC(i)]++] = FSA(i);
        e[v[FSA(i)]+d[FSA(i)]++] = FSD(i);
        e[v[FSD(i)]+d[FSD(i)]++] = FSA(i);
    }

    for (i = 0; i < k; ++i)
    {
        j = FSB((i+1)%k);
	e[v[FSB(i)]+d[FSB(i)]++] = j;
	e[v[j]+d[j]++] = FSB(i);
    }

    for (i = 0; i < k-1; ++i)
    {
	e[v[FSC(i)]+d[FSC(i)]++] = FSC(i+1);
	e[v[FSC(i+1)]+d[FSC(i+1)]++] = FSC(i);
    }

    for (i = 0; i < k-1; ++i)
    {
	e[v[FSD(i)]+d[FSD(i)]++] = FSD(i+1);
	e[v[FSD(i+1)]+d[FSD(i+1)]++] = FSD(i);
    }

    e[v[FSD(0)]+d[FSD(0)]++] = FSC(k-1);
    e[v[FSC(k-1)]+d[FSC(k-1)]++] = FSD(0);
    e[v[FSC(0)]+d[FSC(0)]++] = FSD(k-1);
    e[v[FSD(k-1)]+d[FSD(k-1)]++] = FSC(0);
}

/**************************************************************************/

static void
makeJohnson(long n, long k, boolean digraph, sparsegraph *sg)
{
    size_t *v;
    int *d,*e,*ep,nv,deg,i,j,s,t,u;
    DYNALLSTAT(int,a,a_sz);
    DYNALLSTAT(int,b,b_sz);

    if (k > n/2) k = n - k;
    if (k < 0) gt_abort(">E genspecialg: bad parameters for -J\n");

    nv = binomial(n,k);
    if (nv > NAUTY_INFINITY-2) gt_abort(">E genspecialg: too big -J\n");
    deg = k*(n-k);

    SG_ALLOC(*sg,nv,nv*(size_t)deg,"genspecialg");
    sg->nv = nv;
    sg->nde = nv*(size_t)deg;
    SG_VDE(sg,v,d,e);

    DYNALLOC1(int,a,a_sz,k,"genspecialg");
    DYNALLOC1(int,b,b_sz,k,"genspecialg");
    preparemarks(n);

    for (i = 0; i < nv; ++i)
    {
	v[i] = i*(size_t)deg;
	d[i] = deg;
	ep = e + v[i];
	unrank(i,k,a);
//{int x;for(x=0;x<k;++x)printf(" %d",a[x]);printf("\n");}
	RESETMARKS;
	for (j = 0; j < k; ++j) MARK(a[j]);

	for (j = 0; j < n; ++j)
	if (ISNOTMARKED(j))
	{
	    for (s = 0; s < k; ++s)
	    {
		for (t = 0; t < k; ++t) b[t] = a[t];
		u = s;
		while (u > 0 && b[u-1] > j)
		{
		    b[u] = b[u-1];
		    --u;
		}
		while (u < k-1 && b[u+1] < j)
		{
		    b[u] = b[u+1];
		    ++u;
		}
		b[u] = j;
//{int x;printf("-");for(x=0;x<k;++x)printf(" %d",b[x]);printf(" (%d)\n",rank(k,b));}
		*(ep++) = rank(k,b);
	    }
	}
    }
}

/**************************************************************************/

static void
makecomplete(long n, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j;
    size_t *v,k;

    if (n < 1 || n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: bad argument for -k\n");

    if (digraph) SG_ALLOC(*sg,n,n*(size_t)n,"genspecialg");
    else         SG_ALLOC(*sg,n,n*(size_t)(n-1),"genspecialg");

    SG_VDE(sg,v,d,e);

    if (digraph)
    {
        sg->nv = n;
	sg->nde = n*(size_t)n;

	for (i = 0, k = 0; i < n; ++i, k += n)
	{
	    d[i] = n;
	    v[i] = k;
	    for (j = 0; j < n; ++j) e[k+j] = j;
        }
    }
    else
    {
        sg->nv = n;
	sg->nde = n*(size_t)(n-1);

	for (i = 0, k = 0; i < n; ++i)
	{
	    d[i] = n-1;
	    v[i] = k;
	    for (j = 0; j < n; ++j)
               if (j != i) e[k++] = j;
        }
    }
}

/**************************************************************************/

static void
makeempty(long n, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i;
    size_t *v;

    if (n < 1 || n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: bad argument for -e\n");

    if (digraph) SG_ALLOC(*sg,n,n,"genspecialg");
    else         SG_ALLOC(*sg,n,0,"genspecialg");

    SG_VDE(sg,v,d,e);

    if (digraph)
    {
        sg->nv = n;
	sg->nde = n;

	for (i = 0; i < n; ++i)
	{
	    d[i] = 1;
	    v[i] = i;
	    e[i] = i;
        }
    }
    else
    {
        sg->nv = n;
	sg->nde = 0;

	for (i = 0; i < n; ++i)
	{
	    d[i] = 0;
	    v[i] = 0;
        }
    }
}

/**************************************************************************/

static void
makehypercube(long deg, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j;
    size_t *v,k,nv;

    if (deg < 0 || deg > 30)
        gt_abort(">E genspecialg: bad argument for -q\n");
    if (digraph)
        gt_abort(">E genspecialg: -zq is not implemented\n");

    nv = 1UL << deg;
    SG_ALLOC(*sg,nv,deg*nv,"genspecialg");

    SG_VDE(sg,v,d,e);

    sg->nv = nv;
    sg->nde = deg*nv;

    for (i = 0, k = 0; i < nv; ++i, k += deg)
    {
	d[i] = deg;
	v[i] = k;
	for (j = 0; j < deg; ++j) e[k+j] = i ^ (1<<j);
    }
}

/**************************************************************************/

static void
maketheta(long *len, int npaths, boolean digraph, sparsegraph *sg)
{
    int i,j,k,n,ntemp,*d,*e;
    size_t *v,ne,etemp;
    boolean hasone;

    hasone = FALSE;
    n = 2;
    ne = 0;
    for (i = 0; i < npaths; ++i)
    {
	if (len[i] < 1)
	    gt_abort(">E genspecialg: -T paths must be at least length 1\n");
	if (len[i] == 1)
	{
	    if (hasone) gt_abort(
                  ">E genspecialg: -T only one path of length 1 allowed\n");
	    hasone = TRUE;
	}
	ntemp = n;
	n += len[i]-1;
	if (n < ntemp)
	    gt_abort(">E genspecialg: -T too many vertices\n");
	etemp = ne;
	ne += len[i];
	if (ne < etemp) gt_abort(">E genspecialg: -T too many edges\n");
    }

    if (n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: -T size is too big\n");

    if (!digraph)
    {
	etemp = ne;
	ne *= 2;
	if (ne < etemp) gt_abort(">E genspecialg: -T too many edges\n");
    }

    SG_ALLOC(*sg,n,ne,"genspecialg");
    SG_VDE(sg,v,d,e);
    sg->nv = n;
    sg->nde = ne;    

    v[0] = 0;
    v[1] = npaths;
    if (digraph)
    {
	v[2] = v[1];
	for (i = 3; i < n; ++i) v[i] = v[i-1] + 1;
    }
    else
    {
	v[2] = v[1] + npaths;
	for (i = 3; i < n; ++i) v[i] = v[i-1] + 2;
    }

    for (i = 0; i < n; ++i) d[i] = 0;

    if (hasone)
    {
	e[v[0]+(d[0]++)] = 1;
	if (!digraph) e[v[1]+(d[1]++)] = 0;
    }

    k = 2;
    for (i = 0; i < npaths; ++i)
    {
	if (len[i] == 1) continue;

	e[v[0]+(d[0]++)] = k;
	if (!digraph) e[v[k]+(d[k]++)] = 0;
	
	for (j = 0; j < len[i]-2; ++j)
	{
	    e[v[k]+(d[k]++)] = k+1;
	    if (!digraph) e[v[k+1]+(d[k+1]++)] = k;
	    ++k;
	}
	e[v[k]+(d[k]++)] = 1;
	if (!digraph) e[v[1]+(d[1]++)] = k;
        ++k;
    }
}

/**************************************************************************/

static void
makegrid(long *dim, int ndim, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j,deg,n,oldn;
    size_t *v,k;
    boolean closed[30];
    int index[30];

    n = 1;
    deg = 0;
    for (i = 0; i < ndim; ++i)
    {
        if (dim[i] >= -1 && dim[i] <= 1)
            gt_abort(">E genspecialg: -G dimensions must be at least 2\n");
	if (dim[i] == 2 && !digraph)
            gt_abort(">E genspecialg: -G dimen 2 is only ok for digraphs\n");

	closed[i] = (dim[i] > 0);
	if (dim[i] < 0) dim[i] = -dim[i];

	oldn = n;
        n *= dim[i];
	if (n < 0 || n / dim[i] != oldn)
	    gt_abort(">E genspecialg: -G size is too big\n");

	if (digraph || dim[i] == 2) ++deg;
        else                        deg += 2;

        index[i] = 0;
    }

    if (n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: -G size is too big\n");

    SG_ALLOC(*sg,n,deg*(size_t)n,"genspecialg");

    SG_VDE(sg,v,d,e);

    sg->nv = n;
    sg->nde = deg*(size_t)n;

    k = 0;
    for (i = 0; i < n; ++i)
    {
	v[i] = k;
	for (j = 0; j < ndim; ++j)
	{
	    if (index[j] < dim[j]-1)
	    {
		++index[j];
		e[k++] = vnumber(dim,index,ndim);
		--index[j];
	    }
	    if (!digraph && index[j] > 0)
	    {
		--index[j];
		e[k++] = vnumber(dim,index,ndim);
		++index[j];
	    }
	    if (closed[j] && index[j] == dim[j]-1)
	    {
		index[j] = 0;
		e[k++] = vnumber(dim,index,ndim);
		index[j] = dim[j]-1;
	    }
	    if (closed[j] && !digraph && index[j] == 0)
	    {
		index[j] = dim[j]-1;
		e[k++] = vnumber(dim,index,ndim);
		index[j] = 0;
	    }
	}

        d[i] = k - v[i];

	for (j = ndim; --j >= 0;)
	{
	    if (index[j] != dim[j]-1)
	    {
		++index[j];
		break;
	    }
	    else
		index[j] = 0;
	}
    }
}

/**************************************************************************/

static void
makecirculant(long n, long *conn, int nconn, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j,deg;
    size_t *v,k;

    if (nconn > 0 && conn[0] <= 0)
        gt_abort(">E genspecialg: -C connections must be nonzero\n");

    for (i = 1; i < nconn; ++i)
	if (conn[i] <= conn[i-1])
	    gt_abort(">E genspecialg: -C connections must be increasing\n");

    if (nconn == 0)
	deg = 0;
    else
    {
        if (digraph)
	{
	    if (conn[nconn-1] >= n) gt_abort(
                 ">E genspecialg: -C connections must be 1..n-1\n");
	    deg = nconn;
	}
	else
        {
            if (conn[nconn-1] > n/2) gt_abort(
                 ">E genspecialg: -C connections must be 1..n/2\n");
	    deg = 2*nconn - (2*conn[nconn-1]==n);
	}
    }

    SG_ALLOC(*sg,n,deg*n,"genspecialg");

    SG_VDE(sg,v,d,e);
    sg->nv = n;
    sg->nde = deg*n;

    for (i = 0; i < n; ++i)
    {
        d[i] = deg;
	v[i] = deg*(size_t)i;
    }
 
    for (i = 0; i < n; ++i)
    {
	k = v[i];
	for (j = 0; j < nconn; ++j)
	{
	    e[k++] = (i + conn[j]) % n;
	    if (!digraph && 2*conn[j] != n)
		e[k++] = (i - conn[j] + n) % n;
	}
    }
}

/**************************************************************************/

static void
makegenpetersen(long n1, long n2, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j,n;
    size_t *v,k;

    if (digraph) gt_abort(">E no digraph version of -P is implemented\n");

    n = 2*n1;
    if (n < 1 || n1 > NAUTY_INFINITY/2-1 || n2 < 1 || 2*n2 >= n1)
	gt_abort(">E -Pm,k needs m>0,0<k<m/2; or m too large\n");

    SG_ALLOC(*sg,n,3UL*n,"genspecialg");

    SG_VDE(sg,v,d,e);
    sg->nv = n;
    sg->nde = 3UL*n;

    for (i = 0; i < n; ++i)
    {
        d[i] = 3;
	v[i] = 3UL*i;
    }

    for (i = 0; i < n1; ++i)
    {
	k = v[i];
	e[k] = (i + 1) % n1;
	e[k+1] = (i + n1 - 1) % n1;
	e[k+2] = i + n1;
    }
    
    for (i = 0; i < n1; ++i)
    {
	k = v[n1+i];
	e[k] = n1 + (i + n2) % n1;
        e[k+1] = n1 + (i - n2 + n1) % n1;
	e[k+2] = i;
    }
} 

/**************************************************************************/

static void
makecompletebipartite(long n1, long n2, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j,n;
    size_t *v,k;

    n = n1 + n2;

    if (n1 < 1 || n2 < 1 || n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: bad argument for -b\n");

    if (digraph) SG_ALLOC(*sg,n,n1*n2,"genspecialg");
    else         SG_ALLOC(*sg,n,2*n1*n2,"genspecialg");

    SG_VDE(sg,v,d,e);

    if (digraph)
    {
	sg->nv = n;
	sg->nde = n1*n2;

	for (i = 0, k = 0; i < n1; ++i)
	{
	    d[i] = n2;
	    v[i] = k;
	    for (j = n1; j < n; ++j) e[k++] = j;
	}
	for (i = n1; i < n; ++i)
	{
	    d[i] = 0;
	    v[i] = 0;
        }
    }
    else
    {
	sg->nv = n;
	sg->nde = 2*n1*n2;

	for (i = 0, k = 0; i < n1; ++i)
	{
	    d[i] = n2;
	    v[i] = k;
	    for (j = n1; j < n; ++j) e[k++] = j;
	}
	for (i = n1; i < n; ++i)
	{
	    d[i] = n1;
	    v[i] = k;
	    for (j = 0; j < n1; ++j) e[k++] = j;
        }
    }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    int n,codetype;
    int argnum,i,j,nreq;
    char *arg,sw;
    boolean badargs,quiet;
    boolean Cswitch,Pswitch,gswitch,sswitch,zswitch,Jswitch,dswitch;
    boolean pswitch,cswitch,eswitch,kswitch,bswitch,Qswitch,Gswitch;
    boolean fswitch,Tswitch;
    long size;
    static FILE *outfile;
    char *outfilename;
    sparsegraph sg;
    boolean usesparse,digraph;
    long Pargs[2],bargs[2],Jargs[2];
    int nPargs,nbargs,nCargs,nGargs,nJargs,nTargs;

    HELP; PUTVERSION;

    gswitch = sswitch = zswitch = Pswitch = FALSE;
    pswitch = cswitch = eswitch = kswitch = FALSE;
    Gswitch = Cswitch = bswitch = Qswitch = FALSE;
    dswitch = Jswitch = fswitch = Tswitch = quiet = FALSE;

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
                else SWBOOLEAN('z',zswitch)
                else SWBOOLEAN('d',dswitch)
                else SWBOOLEAN('q',quiet)
                else SWLONG('p',pswitch,size,"genspecialg -p")
                else SWLONG('c',cswitch,size,"genspecialg -c")
                else SWLONG('e',eswitch,size,"genspecialg -e")
                else SWLONG('k',kswitch,size,"genspecialg -k")
                else SWLONG('f',fswitch,size,"genspecialg -f")
                else SWLONG('Q',Qswitch,size,"genspecialg -Q")
                else SWSEQUENCE('b',",",bswitch,bargs,2,
                                nbargs,"genspecialg -b")
                else SWSEQUENCE('J',",",Jswitch,Jargs,2,
                                nJargs,"genspecialg -J")
                else SWSEQUENCE('P',",",Pswitch,Pargs,2,
                                nPargs,"genspecialg -P")
                else SWSEQUENCE('C',",",Cswitch,args,MAXARGS,
                                nCargs,"genspecialg -C")
                else SWSEQUENCE('G',",",Gswitch,args,30,
                                nGargs,"genspecialg -G")
                else SWSEQUENCE('T',",",Tswitch,args,MAXARGS,
                                nTargs,"genspecialg -T")
                else badargs = TRUE;
            }
        }
        else
        {
            ++argnum;
            if (argnum == 1) outfilename = arg;
            else             badargs = TRUE;
        }
    }

    if ((gswitch!=0) + (sswitch!=0) + (zswitch!=0) > 1)
        gt_abort(">E genspecialg: -gsz are incompatible\n");

    if ((gswitch!=0) + (sswitch!=0) + (dswitch!=0) > 1)
        gt_abort(">E genspecialg: -gsd are incompatible\n");

    nreq = (pswitch!=0) + (cswitch!=0) + (eswitch!=0) + (kswitch!=0)
           + (bswitch!=0) + (Qswitch!=0) + (Pswitch!=0) + (fswitch!=0)
           + (Cswitch!= 0) + (Gswitch!=0) + (Jswitch!=0)
           + (Tswitch!= 0);
    if (nreq > 1)
        gt_abort(">E genspecialg: must have exactly one of -bcfkpCGJPQT\n");
    else if (nreq < 1) 
	badargs = TRUE;

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (gswitch)      codetype = GRAPH6;
    else if (zswitch) codetype = DIGRAPH6;
    else              codetype = SPARSE6;

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

    SG_INIT(sg);

    if (pswitch)
        makepath(size,zswitch,&sg);
    else if (cswitch)
        makecycle(size,zswitch,&sg);
    else if (kswitch)
        makecomplete(size,zswitch,&sg);
    else if (eswitch)
        makeempty(size,zswitch,&sg);
    else if (Qswitch)
        makehypercube(size,zswitch,&sg);
    else if (bswitch)
    {
	if (nbargs != 2) gt_abort(">E genspecialg: -b needs two arguments\n");
        makecompletebipartite(bargs[0],bargs[1],zswitch,&sg);
    }
    else if (Jswitch)
    {
	if (nJargs != 2) gt_abort(">E genspecialg: -J needs two arguments\n");
        makeJohnson(Jargs[0],Jargs[1],zswitch,&sg);
    }
    else if (Pswitch)
    {
	if (nPargs != 2) gt_abort(">E genspecialg: -P needs two arguments\n");
        makegenpetersen(Pargs[0],Pargs[1],zswitch,&sg);
    }
    else if (Cswitch)
        makecirculant(args[0],args+1,nCargs-1,zswitch,&sg);
    else if (Gswitch)
    {
	if (nGargs < 2)
            gt_abort(">E genspecialg: -G needs at least two arguments\n");
        makegrid(args,nGargs,zswitch,&sg);
    }
    else if (Tswitch)
    {
	if (nTargs < 1)
            gt_abort(">E genspecialg: -T needs at least one argument\n");
        maketheta(args,nTargs,zswitch,&sg);
    }
    else if (fswitch)
	makeflowersnark(size,zswitch,&sg);

    sortlists_sg(&sg);
    if (dswitch)                   writedread(outfile,&sg,zswitch);
    else if (codetype == GRAPH6)   writeg6_sg(outfile,&sg);
    else if (codetype == DIGRAPH6) writed6_sg(outfile,&sg);
    else                           writes6_sg(outfile,&sg);

    if (!quiet)
        fprintf(stderr,">Z %d vertices %lu edges\n",sg.nv,
                       (unsigned long)(zswitch ? sg.nde : sg.nde/2));

    exit(0);
}
