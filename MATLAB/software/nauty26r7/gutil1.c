/* gutil1.c: Some graph utilities. */

#include "gtools.h"
#include "gutils.h"

/**************************************************************************/

void
degstats(graph *g, int m, int n, unsigned long *edges, int *mindeg,
     int *mincount, int *maxdeg, int *maxcount, boolean *eulerian)
/* Compute degree-related graph properties.
   *edges = number of edges
   *mindeg, *mincount = minimum degree and how many there are
   *maxdeg, *maxcount = maximum degree and how many there are
   *eulerian = whether the graph has only even degrees
*/
{
    setword *pg;
    int i,j,d,dor;
    int mind,mindc,maxd,maxdc;
    unsigned long ned;

    mind = n;
    mindc = 0;
    maxd = 0;
    maxdc = 0;
    ned = 0;
    dor = 0;

    pg = (setword*)g;
    for (i = 0; i < n; ++i)
    {
        d = 0;
        for (j = 0; j < m; ++j, ++pg)
            if (*pg) d += POPCOUNT(*pg);

        if (d == mind)
            ++mindc;
        else if (d < mind)
        {
            mind = d;
            mindc = 1;
        }

        if (d == maxd)
            ++maxdc;
        else if (d > maxd)
        {
            maxd = d;
            maxdc = 1;
        }

        dor |= d;
        ned += d;
    }

    *mindeg = mind;
    *mincount = mindc;
    *maxdeg = maxd;
    *maxcount = maxdc;
    *edges = ned / 2;
    *eulerian = (dor & 1) == 0;
}

/**************************************************************************/

void
degstats2(graph *g, boolean digraph, int m, int n,
     unsigned long *edges, int *loops,
     int *minindeg, int *minincount, int *maxindeg, int *maxincount,
     int *minoutdeg, int *minoutcount, int *maxoutdeg, int *maxoutcount,
     boolean *eulerian)
/* Compute degree-related graph properties.
   *edges = number of edges (including loops), directed edges for digraphs
   *loops = number of loops
   *minindeg, *minincount = minimum in-degree and how many there are
   *maxindeg, *maxincount = maximum in-degree and how many there are
   *minoutdeg,*minoutcount,*maxoutdeg,*maxoutcount = similarly for out-degree
   *eulerian = whether the undirected graph has only even degrees,
               or the directed graph has indegree=outdegree at each vertex.
   A loop contributes 2 to the degrees of an undirected graph and
     1 to each degree for a directed graph.
*/
{
    setword *pg;
    int i,j,d,dor;
    int mind,mindc,maxd,maxdc;
    unsigned long ned;
    int nloops;
#if MAXN
    int indeg[MAXN];
    int outdeg[MAXN];
#else
    DYNALLSTAT(int,indeg,indeg_sz);
    DYNALLSTAT(int,outdeg,outdeg_sz);
#endif

#if !MAXN
    if (digraph)
    {
        DYNALLOC1(int,indeg,indeg_sz,n,"degstats2");
        DYNALLOC1(int,outdeg,outdeg_sz,n,"degstats2");
    }
#endif

    if (!digraph)
    {
        mind = n+2;
        mindc = 0;
        maxd = 0;
        maxdc = 0;
        ned = 0;
        dor = 0;
        nloops = 0;

        pg = (setword*)g;
        for (i = 0; i < n; ++i)
        {
            d = 0;
	    if (ISELEMENT(pg,i))
	    {
		++d;
		++nloops;
	    }
            for (j = 0; j < m; ++j, ++pg)
                if (*pg) d += POPCOUNT(*pg);

            if (d == mind)
                ++mindc;
            else if (d < mind)
            {
                mind = d;
                mindc = 1;
            }
    
            if (d == maxd)
                ++maxdc;
            else if (d > maxd)
            {
                maxd = d;
                maxdc = 1;
            }
    
            dor |= d;
            ned += d;
        }

        *minindeg = *minoutdeg = mind;
        *minincount = *minoutcount = mindc;
        *maxindeg = *maxoutdeg = maxd;
        *maxincount = *maxoutcount = maxdc;
        *edges = ned / 2;
        *eulerian = (dor & 1) == 0;
        *loops = nloops;
    }
    else
    {
	for (i = 0; i < n; ++i) indeg[i] = outdeg[i] = 0;

        nloops = 0;
	ned = 0;
	for (i = 0, pg = (setword*)g; i < n; ++i, pg += m)
	{
	    if (ISELEMENT(pg,i)) ++nloops;
	    for (j = -1; (j = nextelement(pg,m,j)) >= 0;)
	    {
		++outdeg[i];
		++indeg[j];
	    }
	    ned += outdeg[i];
	}
	*edges = ned;
	*loops = nloops;

        mind = maxd = indeg[0];
        mindc = maxdc = 1;

	for (i = 1; i < n; ++i)
	{
	    d = indeg[i];

            if (d == mind)
                ++mindc;
            else if (d < mind)
            {
                mind = d;
                mindc = 1;
            }
    
            if (d == maxd)
                ++maxdc;
            else if (d > maxd)
            {
                maxd = d;
                maxdc = 1;
	    }
	}
	*minindeg = mind;
        *minincount = mindc;
	*maxindeg = maxd;
        *maxincount = maxdc;

        mind = maxd = outdeg[0];
        mindc = maxdc = 1;

	for (i = 1; i < n; ++i)
	{
	    d = outdeg[i];

            if (d == mind)
                ++mindc;
            else if (d < mind)
            {
                mind = d;
                mindc = 1;
            }
    
            if (d == maxd)
                ++maxdc;
            else if (d > maxd)
            {
                maxd = d;
                maxdc = 1;
	    }
	}
	*minoutdeg = mind;
        *minoutcount = mindc;
	*maxoutdeg = maxd;
        *maxoutcount = maxdc;

	for (i = 0; i < n; ++i)
	    if (indeg[i] != outdeg[i]) break;
	*eulerian = (i == n);
    }
}

/*********************************************************************/

boolean
isconnected1(graph *g, int n)
/* test if g is connected (m=1) */
{
    setword seen,expanded,toexpand;
    int i;

    seen = bit[0];
    expanded = 0;

    while ((toexpand = (seen & ~expanded)) != 0)
    {
        i = FIRSTBITNZ(toexpand);
        expanded |= bit[i];
        seen |= g[i];
    }

    return  POPCOUNT(seen) == n;
}

/**************************************************************************/

boolean
isconnected(graph *g, int m, int n)
/* Test if g is connected */
{
    int i,head,tail,w;
    set *gw;
#if MAXN
    int queue[MAXN],visited[MAXN];
#else
    DYNALLSTAT(int,queue,queue_sz);
    DYNALLSTAT(int,visited,visited_sz);
#endif

    if (m == 1) return isconnected1(g,n);

#if !MAXN
    DYNALLOC1(int,queue,queue_sz,n,"isconnected");
    DYNALLOC1(int,visited,visited_sz,n,"isconnected");
#endif

    for (i = 0; i < n; ++i) visited[i] = 0;

    queue[0] = 0;
    visited[0] = 1;

    head = 0;
    tail = 1;
    while (head < tail)
    {
        w = queue[head++];
        gw = GRAPHROW(g,w,m);
        for (i = -1; (i = nextelement(gw,m,i)) >= 0;)
        {
            if (!visited[i])
            {
                visited[i] = 1;
                queue[tail++] = i;
            }
        }
    }

    return tail == n;
}

/**************************************************************************/
 
boolean
issubconnected(graph *g, set *sub, int m, int n)
/* Test if the subset of g induced by sub is connected. Empty is connected. */
{
    int i,head,tail,w,subsize;
    set *gw;
#if MAXN
    int queue[MAXN],visited[MAXN];
    setword subw[MAXM];
#else
    DYNALLSTAT(int,queue,queue_sz);
    DYNALLSTAT(int,visited,visited_sz);
    DYNALLSTAT(set,subw,subw_sz);
 
    DYNALLOC1(int,queue,queue_sz,n,"issubconnected");
    DYNALLOC1(int,visited,visited_sz,n,"issubconnected");
    DYNALLOC1(set,subw,subw_sz,m,"issubconnected");
#endif
 
    subsize = 0;
    for (i = 0; i < m; ++i) subsize += (sub[i] ? POPCOUNT(sub[i]) : 0);
 
    if (subsize <= 1) return TRUE;
 
    for (i = 0; i < n; ++i) visited[i] = 0;
 
    i = nextelement(sub,m,-1);
    queue[0] = i;
    visited[i] = 1;
 
    head = 0;
    tail = 1;
    while (head < tail)
    {
        w = queue[head++];
        gw = GRAPHROW(g,w,m);
        for (i = 0; i < m; ++i) subw[i] = gw[i] & sub[i];
 
        for (i = -1; (i = nextelement(subw,m,i)) >= 0;)
        {
            if (!visited[i])
            {
                visited[i] = 1;
                queue[tail++] = i;
            }
        }
    }
 
    return tail == subsize;
}
 
/**********************************************************************/
 
boolean
isbiconnected1(graph *g, int n)
/* Test if g is biconnected; version for m=1. */
{
    int sp,v,w;
    setword sw;
    int numvis;
    setword visited;
    int num[WORDSIZE],lp[WORDSIZE],stack[WORDSIZE];
 
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
 
boolean
isbiconnected(graph *g, int m, int n)
/* test if g is biconnected */
{
    int sp,v,vc;
    int numvis;
    set *gv;
#if MAXN
    int num[MAXN],lp[MAXN],stack[MAXN];
#else
    DYNALLSTAT(int,num,num_sz);
    DYNALLSTAT(int,lp,lp_sz);
    DYNALLSTAT(int,stack,stack_sz);
#endif

    if (n <= 2) return FALSE;
    if (m == 1) return isbiconnected1(g,n);

#if !MAXN
    DYNALLOC1(int,num,num_sz,n,"isbiconnected");
    DYNALLOC1(int,lp,lp_sz,n,"isbiconnected");
    DYNALLOC1(int,stack,stack_sz,n,"isbiconnected");
#endif
 
    num[0] = 0;
    for (v = 1; v < n; ++v) num[v] = -1;
    lp[0] = 0;
    numvis = 1;
    sp = 0;
    v = 0;
    vc = -1;
    gv = (set*)g;
 
    for (;;)
    {
        vc = nextelement(gv,m,vc);
        if (vc < 0)
        {
            if (sp <= 1)  return numvis == n;
            vc = v;
            v = stack[--sp];
            gv = GRAPHROW(g,v,m);
            if (lp[vc] >= num[v])  return FALSE;
            if (lp[vc] < lp[v])    lp[v] = lp[vc];
        }
        else if (num[vc] < 0)
        {
            stack[++sp] = vc;
            v = vc;
            gv = GRAPHROW(g,v,m);
            vc = -1;
            lp[v] = num[v] = numvis++;
        }
        else if (vc != v)
        {
            if (num[vc] < lp[v])  lp[v] = num[vc];
        }
    }
}

/**************************************************************************/

boolean
twocolouring(graph *g, int *colour, int m, int n)
/* If g is bipartite, set colour[*] to 0 or 1 to indicate an example
   of 2-colouring and return TRUE.  Otherwise return FALSE.
   Colour 0 is assigned to the first vertex of each component.  */
{
    int i,head,tail,v,w,need;
    set *gw;
    setword xg;
#if MAXN
    int queue[MAXN];
#else
    DYNALLSTAT(int,queue,queue_sz);
#endif

#if !MAXN
    DYNALLOC1(int,queue,queue_sz,n,"twocolouring");
#endif

    for (i = 0; i < n; ++i) colour[i] = -1;

    if (m == 1)
    {
        for (v = 0; v < n; ++v)
            if (colour[v] < 0)
            {
                queue[0] = v;
                colour[v] = 0;
    
                head = 0;
                tail = 1;
                while (head < tail) 
                {
                    w = queue[head++];
                    need = 1 - colour[w];
                    xg = g[w];
		    while (xg)
                    {
			TAKEBIT(i,xg);
                        if (colour[i] < 0)
                        {
                            colour[i] = need;
                            queue[tail++] = i;
                        }
                        else if (colour[i] != need)
                            return FALSE;
                    }
                }
            }
    }
    else
    {
        for (v = 0; v < n; ++v)
            if (colour[v] < 0)
            {
                queue[0] = v;
                colour[v] = 0;
    
                head = 0;
                tail = 1;
                while (head < tail) 
                {
                    w = queue[head++];
                    need = 1 - colour[w];
                    gw = GRAPHROW(g,w,m);
                    for (i = -1; (i = nextelement(gw,m,i)) >= 0;)
                    {
                        if (colour[i] < 0)
                        {
                            colour[i] = need;
                            queue[tail++] = i;
                        }
                        else if (colour[i] != need)
                            return FALSE;
                    }
                }
            }
     }

    return TRUE;
}

/**************************************************************************/

boolean
isbipartite(graph *g, int m, int n)
/* Test if g is bipartite */
{
#if MAXN
    int colour[MAXN];
#else
    DYNALLSTAT(int,colour,colour_sz);

    DYNALLOC1(int,colour,colour_sz,n,"isbipartite");
#endif

    return twocolouring(g,colour,m,n);
}

/**************************************************************************/

int
bipartiteside(graph *g, int m, int n)
/* If g is not bipartite, return 0.
   Otherwise return the size of the smallest of the two parts of
    some 2-coluring.  Note that this is not isomorphism-invariant
    if g is disconnected. */
{
    boolean isbip;
    int i,sz;
#if MAXN
    int colour[MAXN];
#else
    DYNALLSTAT(int,colour,colour_sz);
#endif

#if !MAXN
    DYNALLOC1(int,colour,colour_sz,n,"isbipartite");
#endif

    isbip = twocolouring(g,colour,m,n);
    if (!isbip) return 0;

    sz = 0;
    for (i = 0; i < n; ++i) sz += colour[i];

    if (sz+sz <= n) return sz;
    else            return n - sz;
}

/**************************************************************************/

int
girth(graph *g, int m, int n)
/* Find the girth of graph g.  0 means acyclic. */
{
    int i,head,tail,v,w;
    int best,c,dw1;
    set *gw;
#if MAXN
    int dist[MAXN],queue[MAXN];
#else   
    DYNALLSTAT(int,queue,queue_sz);
    DYNALLSTAT(int,dist,dist_sz);
    
    DYNALLOC1(int,queue,queue_sz,n,"girth");
    DYNALLOC1(int,dist,dist_sz,n,"girth");
#endif  
    
    best = n+3;

    for (v = 0; v < n; ++v)
    {
        for (i = 0; i < n; ++i) dist[i] = -1;

        queue[0] = v;
        dist[v] = 0;

        head = 0;
        tail = 1;
        while (head < tail)
        {
            w = queue[head++];
            gw = GRAPHROW(g,w,m);
            dw1 = dist[w] + 1;
            for (i = -1; (i = nextelement(gw,m,i)) >= 0;)
            {
                if (dist[i] < 0)
                {
                    dist[i] = dw1;
                    queue[tail++] = i;
                }
                else if (dist[i] >= dist[w])
                {
                    c = dw1 + dist[i];
                    if (c < best) best = c;
                    if ((c & 1) != 0 || c > best) break;
                }
            }
            if (i >= 0) break;
        }
        if (best == 3) return 3;
    }

    return (best > n ? 0 : best);
}

/**************************************************************************/

void
find_dist(graph *g, int m, int n, int v, int *dist)
/* Put in dist[0..n-1] the distance of each vertex from v.
   Vertices in a different component are given the distance n. */
{
    int i,head,tail,w;
    set *gw;
#if MAXN
    int queue[MAXN];
#else   
    DYNALLSTAT(int,queue,queue_sz);
#endif  
    
#if !MAXN
    DYNALLOC1(int,queue,queue_sz,n,"isconnected");
#endif  
    
    for (i = 0; i < n; ++i) dist[i] = n;

    queue[0] = v;
    dist[v] = 0;

    head = 0;
    tail = 1;
    while (tail < n && head < tail)
    {
        w = queue[head++];
        gw = GRAPHROW(g,w,m);
        for (i = -1; (i = nextelement(gw,m,i)) >= 0;)
        {
            if (dist[i] == n)
            {
                dist[i] = dist[w] + 1;
                queue[tail++] = i;
            }
        }
    }
}

/**************************************************************************/

void
find_dist2(graph *g, int m, int n, int v, int w, int *dist)
/* Put in dist[0..n-1] the distance of each vertex from {v,w}.
   Vertices in a different component are given the distance n. */
{
    int i,head,tail,x;
    set *gx;
#if MAXN
    int queue[MAXN];
#else   
    DYNALLSTAT(int,queue,queue_sz);
#endif  
    
#if !MAXN
    DYNALLOC1(int,queue,queue_sz,n,"isconnected");
#endif  
    
    for (i = 0; i < n; ++i) dist[i] = n;

    queue[0] = v;
    queue[1] = w;
    dist[v] = dist[w] = 0;

    head = 0;
    tail = 2;
    while (tail < n && head < tail)
    {
        x = queue[head++];
        gx = GRAPHROW(g,x,m);
        for (i = -1; (i = nextelement(gx,m,i)) >= 0;)
        {
            if (dist[i] == n)
            {
                dist[i] = dist[x] + 1;
                queue[tail++] = i;
            }
        }
    }
}

/**************************************************************************/

void
diamstats(graph *g, int m, int n, int *radius, int *diameter)
/* Find the radius and diameter.  Both -1 if g is disconnected.
   We use an O(mn) algorithm, which is pretty disgraceful. */
{
    int v,i,head,tail,w;
    int ecc,diam,rad;
    set *gw;
#if MAXN
    int queue[MAXN],dist[MAXN];
#else
    DYNALLSTAT(int,queue,queue_sz);
    DYNALLSTAT(int,dist,dist_sz);
#endif

    /* if (m == 1) {diamstats1(g,n,radius,diameter); return; } */

#if !MAXN
    DYNALLOC1(int,queue,queue_sz,n,"isconnected");
    DYNALLOC1(int,dist,dist_sz,n,"isconnected");
#endif

    diam = -1;
    rad = n;

    for (v = 0; v < n; ++v)
    {
        for (i = 0; i < n; ++i) dist[i] = -1;

        queue[0] = v;
        dist[v] = 0;

        head = 0;
        tail = 1;
        while (tail < n && head < tail)
        {
            w = queue[head++];
            gw = GRAPHROW(g,w,m);
            for (i = -1; (i = nextelement(gw,m,i)) >= 0;)
            {
                if (dist[i] < 0)
                {
                    dist[i] = dist[w] + 1;
                    queue[tail++] = i;
                }
            }
        }

        if (tail < n)
        {
            *diameter = *radius = -1;
            return;
        }

        ecc = dist[queue[n-1]];

        if (ecc > diam) diam = ecc;
        if (ecc < rad)  rad  = ecc;
    }

    *diameter = diam;
    *radius = rad;
}

/**************************************************************************/

static long
maxclnode1(graph *g, setword cliq, setword cov, int maxv)
/* Internal search node.  cov has all the vertices outside cliq that
 * cover all of cliq.  maxv is the last vertex of cliq.
 */
{
    long ans;
    int i;
    setword w;

    if (cov == 0) return 1;

    ans = 0;
    w = cov & BITMASK(maxv);
    while (w)
    {
        TAKEBIT(i,w);
	ans += maxclnode1(g,cliq|bit[i],cov&g[i]&~bit[i],i);
    }
    return ans;
}

long
maxcliques(graph *g, int m, int n)
/* Find the number of maximal cliques */
{
    int i;
    long ans;

    if (m == 1)
    {
        ans = 0;
	for (i = 0; i < n; ++i)
	    ans += maxclnode1(g,bit[i],g[i],i);
    }
    else
    {
	fprintf(stderr,">E maxcliques() is only implemented for m=1\n");
	exit(1);
    }

    return ans;
}
