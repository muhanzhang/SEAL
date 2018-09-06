/* cubhamg.c : pick those inputs that are nonhamiltonian and
                have max degree <= 3.

 Usage:
cubhamg [-#] [-v|-V] [-n#-#|-y#-#|-i|-I|-o|-x|-e|-E] [-b|-t] [infile [outfile]]

        infile is the name of the input file in graph6/sparse6 format
        outfile is the name of the output file in the same format

	stdin and stdout are the defaults for infile and outfile

	The output file will have a header >>graph6<< or >>sparse6<<
        if and only if the input file does.

        Optional switches:

        -#  A parameter useful for tuning (default 100)
	-v  Report nonhamiltonian graphs and noncubic graphs
	-V  .. in addition give a cycle for the hamiltonian ones
	-n#-#  If the two numbers are v and i, then the i-th edge
	    out of vertex v is required to be not in the cycle.
	    It must be that i=1..3 and v=0..n-1.
	-y#-#  If the two numbers are v and i, then the i-th edge
	    out of vertex v is required to be in the cycle.
	    It must be that i=1..3 and v=0..n-1.
            You can use any number of -n/-y switches to force
            edges.  Out of range first arguments are ignored.
            If -y and -n give same edge, -y wins.
        -i  Test + property: for each edge e, there is a hamiltonian
            cycle using e.
	-I  Test ++ property: for each pair of edges e,e', there is
            a hamiltonian cycle which uses both e and e'.
        -o  Test - property: for each edge e, there is a hamiltonian 
            cycle avoiding e.
        -x  Test +- property: for each pair of edges e,e', there is
            a hamiltonian cycle which uses e but avoids e'.
        -e  Test 3/4 property: for each edge e, at least 3 of the 4
            paths of length 3 passing through e lie on hamiltonian cycles.
        -E  Test 3/4+ property: for each edge e failing the 3/4 property,
            all three ways of joining e to the rest of the graph are
            hamiltonian avoiding e.
        -T# Specify a timeout, being a limit on how many search tree
            nodes are made.  If the timeout occurs, the graph is 
            written to the output as if it is nonhamiltonian.
        -R# Specify the number of repeat attempts for each stage.
        -F  Analyze covering paths from 2 or 4 vertices of degree 2.

	-b  Require biconnectivity
        -t  Require triconnectivity  (note: quadratic algorithm)

        -y, -n, -#, -R and -T are ignored for -i, -I, -x, -o, -e, -E, -F

	B. D. McKay, Nov 1995 + Aug 1996 + Feb 2002 + Jul 2008 + Nov 2015

**************************************************************************/

#ifndef MAXN
#define MAXN 152  /* 2 more than largest graph size! */
#endif

#if MAXN==0
#error MAXN cannot be zero for cubhamg
#endif

#include "gtools.h" 
#include "naurng.h"

/**************************************************************************/

#define RANPERM TRUE

/* cubham.h */

#define MAXNE ((3 * MAXN) / 2)
#define FALSE 0
#define TRUE 1
#define YES 1
#define DUNNO 0
#define NO (-1)
#define HABORT 5

#define BADLIM 100  /* limit on number of bad things to report per graph */
#define MAXA  100   /* max number of -y or -n switches */

typedef int cubgraph[MAXN][4];
typedef int vertvec[MAXN];
typedef int edgevec[MAXNE+1];
static long nodecount,maxnodes,totalnodes;
static long timeout;
static long repeats;
static int verbose;

#define NO_LIMIT 0x7FFFFFFFL

static long standard[]
  = {30,40,50,60,100,200,300,400,500,1000,2000,3000,5000,
     10000,20000,30000,100000,300000,1000000};
#define NUMMAXNODES (sizeof(standard)/sizeof(long))
static nauty_counter numtries[NUMMAXNODES+1];

/* cubham.c */

/***************************************************************************

        *****   Instructions for using cubham.  *****

    cubham can find hamiltonian cycles in graphs of maximum 
    degree at most 3, with specifed edges either required or forbidden.  
    To use it:

    (1) Check that MAXN (defined in cubham.h) is at least equal to the
    number of vertices.

    (2) #include cubham.h.

    (3) The cubic graph must be stored in an object of type cubgraph.  Say
    you call it g.  The vertices are numbered 0,1,..., and the neighbours of
    vertex i are g[i][0], g[i][1] and g[i][2].  Entries of the form g[i][3]
    are unused.  If the degree is less than 3, -1 is used as padding.

    (4) Call   cubinit(g,eno,v1,v2,nv,&ne), where
        g     = the cubic graph                         (input)
        eno   = another object of type cubgraph         (output)
        v1,v2 = objects of type edgevec                 (output)
        nv    = the number of vertices                  (input)
	ne    = the number of edges                     (output)

    This numbers the edges 0,1,...ne-1, and defines
        eno[i][j]      = the number of the edge {i, g[i][j]}
        {v1[j], v2[j]} = the vertices of the j-th edge.

    (5) Call  cubham(g,eno,initclass,v1,v2,cycle,outclass,nv,ne), where
        g,eno,v1,v2,nv = as above                       (input)
        initclass = initial edge classification         (input)
        outclass = final edge classifiction             (output)
        cycle = the hamiltonian cycle found, if any     (output)

    The value returned by cubham is either YES (hamiltonian cycle found)
        or NO (there isn't any).

    The initial edge classification is specified by you in the edgevec
    initclass[].  For each edge j, (0 <= j < 3*nv/2), set
    initclass[j] = NO, YES or DUNNO, if edge j must not be in the cycle,
    must be in the cycle, or don't care, respectively.  Passing NULL as
    the initclass parameter is equivalent to setting each edge to DUNNO.
    All initial classifications are legal, even those clearly impossible.

    The final edge classification is set by cubham in the edgevec
    outclass[], if a hamiltonian cycle is found.  Each entry should be
    either NO or YES.  No final classification is provided if NULL is passed
    as the outclass parameter.

    The hamiltonian cycle itself, if any, is returned as cycle[0], cycle[1],
    ..., cycle[nv-1], cycle[0].  If the cycle is not needed, you can pass
    NULL for this parameter.

    Step (4) need not be repeated if the same graph is processed again with
    a different initial classification.

****************************************************************************/

#define POP(x) (onstack[x = *(--stackptr)] = 0)
#define PUSH(x) if(onstack[x]!=stacklev){onstack[x]=stacklev;*(stackptr++)=x;}
#define RESETSTACK {stacklev++; stackptr = stack;}

typedef struct
{
    edgevec class;
    vertvec din,dout,farend;
} nodedata;

static nodedata hcnodat;
static cubgraph eno;
static vertvec onstack,stack;       /* stack contains vertex numbers */
static int *stackptr,stacklev;      /* stackptr points above top */
static int classstack[4*MAXNE];      /* stack of classifications */
   /* x >= 0        : edge number
     (x < 0 above y : farend[-x-1] = y  */
static int *classstackptr;       /* points above top of classstack */
static int classout(cubgraph,nodedata*,int,int,int);
static int classin(cubgraph,cubgraph,nodedata*,int,int,int,int*,int);

#define MAXES 0

#if MAXES
static int maxlevel,maxclassstack;
#endif

static void
dummy(void)
{
}

static void
check_it(int index, cubgraph g, cubgraph eno, edgevec v1, edgevec v2,
         int *din, int *dout, int *class, int *farend, int nin, int nv,
         int stable)
/* Check some things */
{
    int xdin[MAXN],xdout[MAXN],xnin,v,i,j,k,l,*gv,xin,xout,has1;

    for (i = 0; i < nv; ++i) xdin[i] = xdout[i] = 0;
    xnin = has1 = 0;

    for (v = 0; v < nv; ++v)
    {
	gv = g[v];
	xin = xout = 0;
	for (i = 0; i < 3; ++i)
	    if (gv[i] >= 0)
	    {
		j = eno[v][i];
		if (class[j] == NO) ++xout;
		if (class[j] == YES) ++xin;
	    }
	if (xout != dout[v] || xin != din[v])
	{
	    fprintf(stderr,">E%d degrees of %d: din,dout=%d,%d really %d,%d\n",
			index,v,din[v],dout[v],xin,xout);
	    dummy();
	}
	if (xin == 1) ++has1;
	xnin += xin;
    }

    xnin /= 2;
    if (xnin != nin)
    {
        fprintf(stderr,">E%d nin=%d actually %d\n",index,nin,xnin);
        dummy();
    }
    if (nin != 0 && !has1)
    {
        fprintf(stderr,">E%d nin=%d has no in=1\n",index,nin);
        dummy();
    }

    for (i = 0; i < nv; ++i)
	if (din[i] == 0)
	{
            if (farend[i] != i)
    	    {
        	fprintf(stderr,">E%d farend[isolate %d]=%d\n",
			index,i,farend[i]);
        	dummy();
	    }
	}
	else if (din[i] == 1)
	{
	    k = -1;
	    j = i;
	    do
	    {
	        for (l = 0; l < 3; ++l)
		    if (g[j][l] >= 0 && g[j][l] != k && class[eno[j][l]] == YES)
			break;
		k = j;
		if (l < 3) j = g[j][l];
	    } while (l < 3);
            if (farend[i] != j)
    	    {
        	fprintf(stderr,">E%d farend[%d]=%d really %d\n",
			index,i,farend[i],j);
        	dummy();
	    }
        }

    if (stable)
	for (i = 0; i < nv; ++i)
	    if ((dout[i] == 1 && din[i] != 2) || (din[i] == 2 && dout[i] != 1)
		|| dout[i] > 1 || din[i] > 2)
	    {
		fprintf(stderr,">E%d din[%d]=%d dout[%d]=%d\n",
		        index,i,din[i],i,dout[i]);
		dummy();
	    }
}

static void
cubinit(cubgraph g, cubgraph eno, edgevec v1, edgevec v2, int nv, int ne)
/* initialise edge numbers, etc. */
{
        int *gpx,*gpy,*enop,x,y,i,j,n,en;

        n = nv;
        en = 0;
        for (x = 0; x < n; ++x)
        {
            gpx = g[x];
            enop = eno[x];
            for (i = 0; i < 3; ++i)
                if ((y = gpx[i]) < 0)
		    enop[i] = ne;
		else if (y > x)
                {
                    v1[en] = x;
                    v2[en] = y;
                    enop[i] = en++;
                }
                else
                {
                    gpy = g[y];
                    for (j = 0; gpy[j] != x; j++)
                        {}
                    enop[i] = eno[y][j];
                }
        }

        if (en != ne)
            fprintf(stderr,"%% cubinit got en=%d when ne=%d\n",en,ne);
}

static int
propagate(cubgraph g, cubgraph eno, nodedata *ndptr, int *nin, int nv)
/* propagate classifications: */
/*   ans = YES, NO or DUNNO */
{
        int v,w,i,status;
        nodedata *np;
        int *gp,*enop,*class,*din,*dout;

        status = DUNNO;
        np = ndptr;
        class = np->class;
        din = np->din;
        dout = np->dout;

        while (status == DUNNO && stackptr > stack)
        {
            POP(v);
            gp = g[v];
            enop = eno[v];
            if (dout[v] == 0)
            {
                if (din[v] == 2)
                {
                    if (class[enop[0]] == DUNNO)      i = 0;
                    else if (class[enop[1]] == DUNNO) i = 1;
                    else                              i = 2;
                    w = gp[i];
                    status = classout(g,np,v,w,enop[i]);
                    PUSH(w);
                }
                else if (din[v] == 3)
                    status = NO;
            }
            else if (dout[v] == 1)
            {
                for (i = 0; i < 3; ++i)
                if (class[enop[i]] == DUNNO)
                {
                    w = gp[i];
                    if ((status = classin(g,eno,np,v,w,enop[i],nin,nv))
                                                                 != DUNNO)
                        break;
                    else
                        PUSH(w);
                }
            }
            else
                status = NO;
        }

        if (status != NO && *nin == nv) return YES;
        else			        return status;
}

static int
classout(cubgraph g, nodedata *nodat, int v, int w, int en) 
/* classify edge en = vw out */
{
        nodedata *np;

        np = nodat;
        ++np->dout[v];
        ++np->dout[w];
        np->class[en] = NO;
	*classstackptr++ = en;
#if MAXES
	if (classstackptr-classstack > maxclassstack)
	    maxclassstack = classstackptr-classstack;
#endif

        return DUNNO;
}

static int
classin(cubgraph g, cubgraph eno, nodedata *nodat,
                                    int v, int w, int en, int *nin, int nv)
/* classify edge en = vw in */
{
        nodedata *np;
        int *farend,*gp,fv,fw,i;

        np = nodat;
        ++np->din[v];
        ++np->din[w];
        np->class[en] = YES;
	*classstackptr++ = en; 
#if MAXES
	if (classstackptr-classstack > maxclassstack)
	    maxclassstack = classstackptr-classstack;
#endif

        ++*nin;
        if (*nin == nv)
	{
	    return DUNNO;
	}

        farend = np->farend;
        fv = farend[v];
        fw = farend[w];
	*classstackptr++ = farend[fv];
	*classstackptr++ = -fv-1;
	*classstackptr++ = farend[fw];
	*classstackptr++ = -fw-1;
#if MAXES
	if (classstackptr-classstack > maxclassstack)
	    maxclassstack = classstackptr-classstack;
#endif
        farend[fv] = fw;
        farend[fw] = fv;

        gp = g[fv];
        if      (gp[0] == fw) i = 0;
        else if (gp[1] == fw) i = 1;
        else if (gp[2] == fw) i = 2;
        else                  return DUNNO;

        i = eno[fv][i];
        if (np->class[i] == DUNNO)
        {
            PUSH(fv);
            PUSH(fw);
            if (*nin == nv - 1)
                return classin(g,eno,np,fv,fw,i,nin,nv);
            else
                return classout(g,np,fv,fw,i);
        }

        return DUNNO;
}

static int
hamnode(cubgraph g, cubgraph eno, edgevec v1, edgevec v2,
	nodedata *nodat, int level, int nin, int nv)
/* main node for recursion */
{
        int i,p,q,status;
        int v,w,en,*gv,*enov;
	int *csptr;

#if MAXES
	if (level > maxlevel) maxlevel = level;
#endif

        if (++nodecount > maxnodes && maxnodes != NO_LIMIT) return HABORT;
        status = propagate(g,eno,nodat,&nin,nv);

        if (status != DUNNO) return status;

        for (v = nv; --v >= 0;)
            if (nodat->din[v] == 1) break;

        if (v < 0) v = 0;

        gv = g[v];
        enov = eno[v];

        for (i = 0; i < 3; ++i)
        {
            en = enov[i];
            if (nodat->class[en] == DUNNO)
            {
                w = gv[i];
                csptr = classstackptr;
                status = classout(g,nodat,v,w,en);
                if (status == YES) break;
	 	if (status == NO)
		{
		    while (classstackptr > csptr)
		    {
			p = *--classstackptr;
			if (p >= 0)
			{
			    if (nodat->class[p] == YES)
			    {
				--nodat->din[v1[p]];
				--nodat->din[v2[p]];
			    }
			    else
			    {
				--nodat->dout[v1[p]];
				--nodat->dout[v2[p]];
			    }
			    nodat->class[p] = DUNNO;
			}
			else
			{
			    q = *--classstackptr;
			    nodat->farend[-p-1] = q;
			}
		    }
                    continue;
		}
                RESETSTACK;
                PUSH(v);
                PUSH(w);
                status = hamnode(g,eno,v1,v2,nodat,level+1,nin,nv);
                if (status == YES) break;
		while (classstackptr > csptr)
		{
		    p = *--classstackptr;
		    if (p >= 0)
		    {
		        if (nodat->class[p] == YES)
		        {
		    	    --nodat->din[v1[p]];
		    	    --nodat->din[v2[p]];
		        }
		        else
		        {
		    	    --nodat->dout[v1[p]];
		    	    --nodat->dout[v2[p]];
		        }
		        nodat->class[p] = DUNNO;
		    }
		    else
		    {
			q = *--classstackptr;
			nodat->farend[-p-1] = q;
		    }
		}
                if (status == HABORT) return HABORT;
            }
        }

        if (status == DUNNO)
            fprintf(stderr,"hamnode returning DUNNO, this can't happen\n");
        return status;
}

static int
cubham(cubgraph g, cubgraph eno, edgevec initclass, edgevec v1, edgevec v2,
       vertvec cycle, edgevec outclass, int nv, int ne) 
/* external interface */
{
        int i,j,status,nin,v,w;

        for (i = ne; --i >= 0;)
            hcnodat.class[i] = DUNNO;
	if (3*nv > 2*ne) hcnodat.class[ne] = NO;

        for (i = nv; --i >= 0;)
        {
            hcnodat.din[i] = hcnodat.dout[i] = 0;
            hcnodat.farend[i] = i;
            onstack[i] = 0;
        }
        nin = 0;
        stacklev = 0;
        RESETSTACK;

	for (i = nv; --i >= 0;)
	{
	    if (g[i][1] < 0) return NO;
	    if (g[i][2] < 0)
	    {
		hcnodat.dout[i] = 1;
		PUSH(i);
	    }
	}

        status = DUNNO;
	classstackptr = classstack;
        if (initclass)
            for (i = 0; i < ne; ++i)
                if (initclass[i] != DUNNO)
                {
                    v = v1[i];
                    w = v2[i];
                    if (initclass[i] == NO)
                    {
                        if (hcnodat.class[i] == YES)
                            status = NO;
                        else if (hcnodat.class[i] == DUNNO)
                        {
                            if (hcnodat.dout[v] == 0)
                            {
                                status = classout(g,&hcnodat,v,w,i);
                                PUSH(v);
                                PUSH(w);
                            }
                            else
                                status = NO;
			}
                    }
                    else if (initclass[i] == YES)
                    {
                        if (hcnodat.class[i] == NO)
                            status = NO;
                        else if (hcnodat.class[i] == DUNNO)
			{
                            if (hcnodat.din[v] < 2)
                            {
                                status = classin(g,eno,&hcnodat,v,w,i,&nin,nv);
                                PUSH(v);
                                PUSH(w);
                            }
                            else
                                status = NO;
			}
                    }

                    if (status != DUNNO) break;
                }

        if (status == DUNNO)
            status = hamnode(g,eno,v1,v2,&hcnodat,0,nin,nv);

        if (status == YES && cycle)
        {
            w = -1;
            v = 0;
            cycle[0] = 0;
            for (i = 1; i < nv; ++i)
            {
                for (j = 0; g[v][j] == w || hcnodat.class[eno[v][j]] != YES; ++j)
                    {}
                w = v;
                v = g[v][j];
                cycle[i] = v;
            }
        }
        if (status == YES && outclass)
            for (i = 0; i < ne; ++i)
                outclass[i] = hcnodat.class[i];

        return status;
}

/********************************************************************/

static int
isham(cubgraph cub,
      int n, int ne, int weight,
      int *vv, int *vi, int nvv,
      int *yy, int *yi, int nyy, int *cyc)
/* test if hamiltonian; optionally return a cycle 
   Forbid the vi[i]-th nbr of vv[i], for i=0..nvv-1 
   Force the yi[i]-th nbr of yy[i], for i=0..nyy-1 
     WARNING: vi[i]/yi[i] is numbered starting at 1  */
{
        int i,j,k;
        int nmax,ch;
        cubgraph cubcopy;
        edgevec v1,v2,initclass,outclass;
        int perm[MAXN],pinv[MAXN];
	double tmp;

#if !RANPERM

        maxnodes = NO_LIMIT;
        nodecount = 0;
        cubinit(cub,eno,v1,v2,n,ne);

        for (i = 0; i < ne; ++i)
            initclass[i] = DUNNO;

	for (i = 0; i < nvv; ++i)
	    if (vv[i] < n) initclass[eno[vv[i]][vi[i]-1]] = NO;
	for (i = 0; i < nyy; ++i)
	    if (yy[i] < n) initclass[eno[yy[i]][yi[i]-1]] = YES;

        ch = cubham(cub,eno,initclass,v1,v2,cyc,outclass,n,ne);

        totalnodes += nodecount;
        ++numtries[0];
       
#else
        ch = HABORT;
	maxnodes = -1;
        for (nmax = 0; ch == HABORT && maxnodes != timeout; ++nmax)
        {
            if (nmax/repeats < NUMMAXNODES)
            {
		tmp = (double)standard[nmax/repeats] * (double)weight 
                                             * (double)n/ 10000.0;
		if (tmp >= (double)NO_LIMIT) maxnodes = NO_LIMIT;
		else                        maxnodes = tmp;
		if (timeout > 0 && timeout < maxnodes) maxnodes = timeout;
            }
            else if (timeout > 0)
		maxnodes = timeout;
	    else
                maxnodes = NO_LIMIT;

	    if (nmax != 0)
	    {
                for (i = n; --i > 0;)
                {
                    k = KRAN(i+1);
                    j = perm[i];
                    perm[i] = perm[k];
                    perm[k] = j;
                }
	    }
	    else
	    {
		for (i = 0; i < n; ++i)
		    perm[i] = i;
	    }

            for (i = 0; i < n; ++i)
            {
                j = perm[i];
                cubcopy[j][0] = cub[i][0] < 0 ? -1 : perm[cub[i][0]];
                cubcopy[j][1] = cub[i][1] < 0 ? -1 : perm[cub[i][1]];
                cubcopy[j][2] = cub[i][2] < 0 ? -1 : perm[cub[i][2]];
            }

            cubinit(cubcopy,eno,v1,v2,n,ne);
            nodecount = 0;

            for (i = 0; i < ne; ++i)
                initclass[i] = DUNNO;

            for (i = 0; i < nvv; ++i)
                if (vv[i] < n) initclass[eno[perm[vv[i]]][vi[i]-1]] = NO;
            for (i = 0; i < nyy; ++i)
                if (yy[i] < n) initclass[eno[perm[yy[i]]][yi[i]-1]] = YES;

            ch = cubham(cubcopy,eno,initclass,v1,v2,cyc,outclass,n,ne);
            totalnodes += nodecount;
            ++numtries[nmax/repeats];
        }

	if (cyc != NULL && ch == YES)
	{
	    for (i = 0; i < n; ++i)
		pinv[perm[i]] = i;
	    for (i = 0; i < n; ++i)
		cyc[i] = pinv[cyc[i]];
	}
#endif

        return ch;
}

/**************************************************************************/

static int
optadd(cubgraph cub, int v1, int v2)
/* v1 and v2 must have degree 2 and be distinct.
   Add edge v1-v2 if not already present; return index of edge in cub[v1]. */
{
    if (cub[v1][0] == v2) return 0;
    if (cub[v1][1] == v2) return 1;
    cub[v1][2] = v2;
    cub[v2][2] = v1;
    return 2;
}

/**************************************************************************/

static void
dofragment(nauty_counter id, cubgraph cub, int n, int ne, int weight)
/* Test for coverage by one or two paths between vertices of degree 2 */
{
    int i,i1,i2,i3,i4;
    int v1,v2,v3,v4,j1,j3;
    int deg2[MAXN],ndeg2;
    int yy[3],yi[3],newne;
    int cyc[MAXN];
    int status;

    ndeg2 = 0;
    for (i = 0; i < n; ++i)
    {
        if (cub[i][0] < 0 || cub[i][1] < 0)
            gt_abort(">E -F forbids degree 0,1\n");

	if (cub[i][2] < 0) deg2[ndeg2++] = i;
    }

    printf("Input " COUNTER_FMT ":",id);
    for (i = 0; i < ndeg2; ++i) printf(" %d",deg2[i]);
    printf("\n");

    printf(" Pairs: ");
    for (i1 = 0; i1 < ndeg2; ++i1)
    for (i2 = i1+1; i2 < ndeg2; ++i2)
    {
	v1 = deg2[i1]; v2 = deg2[i2];
	j1 = optadd(cub,v1,v2);
	yy[0] = v1; yi[0] = j1+1;
	newne = ne + (j1==2);
	status = isham(cub,n,newne,weight,NULL,NULL,0,yy,yi,1,cyc);
        if (status == HABORT)
	    printf(" T%d-%d",v1,v2);
        if (status == NO)
	    printf(" N%d-%d",v1,v2);
        else 
        {
	    printf(" Y%d-%d",v1,v2);
            if (verbose > 1)
	    {
		printf("[");
	        for (i = 0; i < n; ++i) printf(" %d",cyc[i]);
	        printf("]\n");
	    }
        }
	cub[v1][2] = cub[v2][2] = -1;
    }
    printf("\n");

    printf(" Quartets: ");
    for (i1 = 0; i1 < ndeg2; ++i1)
    for (i2 = i1+1; i2 < ndeg2; ++i2)
    for (i3 = i1+1; i3 < ndeg2; ++i3)
    for (i4 = i3+1; i4 < ndeg2; ++i4)
    {
        if (i3 == i2 || i4 == i2) continue;

	v1 = deg2[i1]; v2 = deg2[i2];
	j1 = optadd(cub,v1,v2);
	v3 = deg2[i3]; v4 = deg2[i4];
	j3 = optadd(cub,v3,v4);
	yy[0] = v1; yi[0] = j1+1;
        yy[1] = v3; yi[1] = j3+1;
	newne = ne + (j1==2) + (j3==2);
	status = isham(cub,n,newne,weight,NULL,NULL,0,yy,yi,2,cyc);
        if (status == HABORT)
	    printf(" T%d-%d,%d-%d",v1,v2,v3,v4);
        if (status == NO)
	    printf(" N%d-%d,%d-%d",v1,v2,v3,v4);
        else
        {
	    printf(" Y%d-%d,%d-%d",v1,v2,v3,v4);
	    if (verbose > 1)
	    {
		printf("[");
	        for (i = 0; i < n; ++i) printf(" %d",cyc[i]);
	        printf("]\n");
	    }
        }
	cub[v1][2] = cub[v2][2] = -1;
	cub[v3][2] = cub[v4][2] = -1;
    }
    printf("\n");
}

/********************************************************************/

static int
hasinout(cubgraph cub,
         int n, int ne, int *x0, int *x1, int *y0, int *y1, int limit)
/* test if cub has in-out (+-) property */
{
	edgevec v1,v2,initclass,outclass;
	set *d0,*di,*dii;
	int i,ii,j,jj;
	int me,nbad;
        DYNALLSTAT(graph,done,done_sz);
	
	me = (ne + WORDSIZE - 1) / WORDSIZE;

        DYNALLOC2(graph,done,done_sz,ne,me,"hasinout");

	d0 = (set*)done;
	EMPTYSET(d0,me);
	for (j = 0; j < ne; ++j)
	    ADDELEMENT(d0,j);

	for (i = 1, di = d0 + me; i < ne; ++i, di += me)
	{
	    for (j = 0; j < me; ++j)
		di[j] = d0[j];
	    DELELEMENT(di,i);
	}
	DELELEMENT(d0,0);

	cubinit(cub,eno,v1,v2,n,ne);
	for (i = 0; i < ne; ++i)
	    initclass[i] = DUNNO;

	maxnodes = NO_LIMIT;
	nbad = 0;

	for (i = 0, di = (set*)done; i < ne; ++i, di += me)
	    for (j = -1; (j = nextelement(di,me,j)) >= 0;)
	    {
		initclass[i] = NO;
		initclass[j] = YES;
		++numtries[0];
		if (cubham(cub,eno,initclass,v1,v2,NULL,
							outclass,n,ne) == NO)
		{
		    x0[nbad] = v1[i]; x1[nbad] = v2[i];
		    y0[nbad] = v1[j]; y1[nbad] = v2[j];
		    ++nbad;
		    if (nbad >= limit) return nbad;
		}
		else
		{
                    for (ii = i, dii = di; ii < ne; ++ii, dii += me)
                    if (outclass[ii] == NO)
                        for (jj = 0; jj < ne; ++jj)
                            if (outclass[jj] == YES)
                                DELELEMENT(dii,jj);
		}

		initclass[i] = DUNNO;
		initclass[j] = DUNNO;
	    }

	return nbad;
}

/**************************************************************************/

static int
hasinin(cubgraph cub,
        int n, int ne, int *x0, int *x1, int *y0, int *y1, int limit)
/* test if cub has in-in (++) property */
{
        edgevec v1,v2,initclass,outclass;
	set *d0,*di,*dii;
	int i,ii,j,jj;
	int me,nbad;
        DYNALLSTAT(graph,done,done_sz);
	
	me = (ne + WORDSIZE - 1) / WORDSIZE;

        DYNALLOC2(graph,done,done_sz,ne,me,"hasinin");

	d0 = (set*)done;
	EMPTYSET(d0,me);
	for (j = 0; j < ne; ++j)
	    ADDELEMENT(d0,j);

	for (i = 1, di = d0 + me; i < ne; ++i, di += me)
	{
	    for (j = 0; j < me; ++j)
		di[j] = d0[j];
	    DELELEMENT(di,i);
	}
	DELELEMENT(d0,0);

	cubinit(cub,eno,v1,v2,n,ne);
	for (i = 0; i < ne; ++i)
	    initclass[i] = DUNNO;

	maxnodes = NO_LIMIT;
	nbad = 0;

	for (i = 0, di = (set*)done; i < ne; ++i, di += me)
	    for (j = i; (j = nextelement(di,me,j)) >= 0;)
	    {
		initclass[i] = YES;
		initclass[j] = YES;
		++numtries[0];
		if (cubham(cub,eno,initclass,v1,v2,NULL,outclass,n,ne) == NO)
		{
                    x0[nbad] = v1[i]; x1[nbad] = v2[i];
                    y0[nbad] = v1[j]; y1[nbad] = v2[j];
                    ++nbad;
                    if (nbad >= limit) return nbad;
		}
		else
	 	{
                    for (ii = i, dii = di; ii < ne; ++ii, dii += me)
                    if (outclass[ii] == YES)
                        for (jj = ii; jj < ne; ++jj)
                            if (outclass[jj] == YES)
                                DELELEMENT(dii,jj);
		}

		initclass[i] = DUNNO;
		initclass[j] = DUNNO;
	    }
	return nbad;
}

/**************************************************************************/

static int
hasin(cubgraph cub, int n, int ne, int *x0, int *x1, int limit)
/* test if cub has "in" property */
{
	edgevec v1,v2,initclass,outclass;
	boolean done[MAXNE];
	int i,ii;
	int nbad;
	
	cubinit(cub,eno,v1,v2,n,ne);
	for (i = 0; i < ne; ++i)
	{
	    initclass[i] = DUNNO;
	    done[i] = FALSE;
	}

	maxnodes = NO_LIMIT;
	nbad = 0;

	for (i = 0; i < ne; ++i)
	    if (!done[i])
	    {
		initclass[i] = YES;
		++numtries[0];
		if (cubham(cub,eno,initclass,v1,v2,NULL,outclass,n,ne) == NO)
		{
		    x0[nbad] = v1[i]; x1[nbad] = v2[i];
		    ++nbad;
		    if (nbad >= limit) return nbad;
		}
		else
		{
                    for (ii = i; ii < ne; ++ii)
                    if (outclass[ii] == YES) done[ii] = TRUE;
		}

		initclass[i] = DUNNO;
	    }
	return nbad;
}

/**************************************************************************/

static boolean
eplus(cubgraph acub, int n, int ne, int x, int y, int *pwhy)
{
        cubgraph cub;
        edgevec v1,v2,initclass;
	int i,a,b,c,d,xy,why;

	if (3*n != 2*ne)
	{
	    fprintf(stderr,
		"cubhamg: eplus() not implemented for noncubic graphs\n");
	    exit(1);
	}

        for (i = 0; i < n; ++i)
        {
            cub[i][0] = acub[i][0];
            cub[i][1] = acub[i][1];
            cub[i][2] = acub[i][2];
        }

	if      (cub[x][0] == y) 
	{
	    xy = 0;
	    a = cub[x][1];
	    b = cub[x][2];
	}
	else if (cub[x][1] == y)
        {
            xy = 1;
            a = cub[x][0];
            b = cub[x][2];
        }
	else
        {
            xy = 2;
            a = cub[x][0];
            b = cub[x][1];
        }

	if      (cub[y][0] == x)
        {
            c = cub[y][1];
            d = cub[y][2];
        }
        else if (cub[y][1] == x)
        {
            c = cub[y][0];
            d = cub[y][2];
        }
        else
        {
            c = cub[y][0];
            d = cub[y][1];
        }

	for (i = 0; i < ne; ++i)
	    initclass[i] = DUNNO;

	why = 0;
	cubinit(cub,eno,v1,v2,n,ne);
	initclass[eno[x][xy]] = NO;
	if (cubham(cub,eno,initclass,v1,v2,NULL,NULL,n,ne) == YES)
	    why |= 1;
	initclass[eno[x][xy]] = DUNNO;

#define CHANGE(z,p,q) if (cub[z][0]==p) cub[z][0]=q;\
     else if (cub[z][1]==p) cub[z][1]=q; else cub[z][2]=q;

	CHANGE(x,b,c);
	CHANGE(y,c,b);
	CHANGE(b,x,y);
	CHANGE(c,y,x);

        cubinit(cub,eno,v1,v2,n,ne);
        initclass[eno[x][xy]] = NO;
        if (cubham(cub,eno,initclass,v1,v2,NULL,NULL,n,ne) == YES)
            why |= 2;
        initclass[eno[x][xy]] = DUNNO;

	CHANGE(x,c,d);
	CHANGE(y,d,c);
	CHANGE(c,x,y);
	CHANGE(d,y,x);

        cubinit(cub,eno,v1,v2,n,ne);
        initclass[eno[x][xy]] = NO;
        if (cubham(cub,eno,initclass,v1,v2,NULL,NULL,n,ne) == YES)
            why |= 4;

	*pwhy = why;
	return why == 7;
}

/**************************************************************************/

static int
hase34(cubgraph cub,
       int n, int ne, int *x0, int *x1, int *why, boolean plus, int limit)
/* test if cub has "e34" property */
{
	edgevec v1,v2,initclass,outclass;
	int ea[4*MAXNE],eb[4*MAXNE],ec[4*MAXNE];
	boolean done[4*MAXNE];
	int count[MAXNE];
	int i,ii;
	int vi,nbad,pluswhy;
	static int pop[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
	
	if (3*n != 2*ne)
	{
	    fprintf(stderr,
		"cubhamg: hase34() not implemented for noncubic graphs\n");
	    exit(1);
	}
	cubinit(cub,eno,v1,v2,n,ne);
	for (i = 0; i < ne; ++i)
	{
	    initclass[i] = DUNNO;
	    count[i] = 0;
	}

	for (i = ii = 0; i < ne; ++i)
	{
	    eb[ii] = eb[ii+1] = eb[ii+2] = eb[ii+3] = i;
	    done[ii] = done[ii+1] = done[ii+2] = done[ii+3] = FALSE;
	    
	    vi = v1[i];
	    if (eno[vi][0] != i) ea[ii++] = eno[vi][0];
	    if (eno[vi][1] != i) ea[ii++] = eno[vi][1];
	    if (eno[vi][2] != i) ea[ii++] = eno[vi][2];
	    ea[ii] = ea[ii-2];
	    ea[ii+1] = ea[ii-1];
	
	    vi = v2[i];
            if (eno[vi][0] != i) ec[ii++] = eno[vi][0];
            if (eno[vi][1] != i) ec[ii++] = eno[vi][1];
            if (eno[vi][2] != i) ec[ii++] = eno[vi][2];
            ec[ii-4] = ec[ii-3] = ec[ii-2]; 
            ec[ii-2] = ec[ii-1];
	}

	maxnodes = NO_LIMIT;
	nbad = 0;

	for (i = 0; i < 4*ne; ++i)
	    if (!done[i])
	    {
		initclass[ea[i]] = YES;
		initclass[eb[i]] = YES;
		initclass[ec[i]] = YES;
		++numtries[0];
		if (cubham(cub,eno,initclass,v1,v2,NULL,
						outclass,n,ne) == YES)
		    for (ii = i; ii < 4*ne; ++ii)
			if (!done[ii] && outclass[ea[ii]] == YES &&
			    outclass[eb[ii]] == YES && outclass[ec[ii]] == YES)
			{
			    done[ii] = TRUE;
			    count[eb[ii]] |= 1 << (ii & 3);
			}
		++numtries[0];
		initclass[ea[i]] = DUNNO;
                initclass[eb[i]] = DUNNO;
                initclass[ec[i]] = DUNNO;
	    }

	pluswhy = 0;
	for (i = 0; i < ne; ++i)
	    if (pop[count[i]] < 3 
                && (!plus || !eplus(cub,n,ne,v1[i],v2[i],&pluswhy)))
	    {
                x0[nbad] = v1[i]; x1[nbad] = v2[i];
		why[nbad] = (pluswhy << 4) | count[i];
                ++nbad;
                if (nbad >= limit) return nbad;
            }

	return nbad;
}

/**************************************************************************/

#if 0
static int
oldhasout(cubgraph cub, int n, int ne, int *x0, int *x1, int limit)
/* test if cub has "out" property */
{
        edgevec v1,v2,initclass,outclass;
        boolean done[MAXNE];
        int i,ii;
        int nbad;
        
        cubinit(cub,eno,v1,v2,n,ne);
        for (i = 0; i < ne; ++i)
        {
            initclass[i] = DUNNO;
            done[i] = FALSE;
        }
 
        maxnodes = NO_LIMIT;
	nbad = 0;
 
        for (i = 0; i < ne; ++i)
            if (!done[i])
            {
                initclass[i] = NO;
                ++numtries[0];
                if (cubham(cub,eno,initclass,v1,v2,NULL,outclass,n,ne) == NO)
                {
                    x0[nbad] = v1[i]; x1[nbad] = v2[i];
                    ++nbad;
                    if (nbad >= limit) return nbad;
                }
		else
		{
		    for (ii = i; ii < ne; ++ii)
                    if (outclass[ii] == NO) done[ii] = TRUE;
		}
 
                initclass[i] = DUNNO;
            }
        return nbad;
}
#endif

/**************************************************************************/

static int
hasout(cubgraph cub, int n, int ne, int weight, int *x0, int *x1, int limit)
/* test if cub has "out" property 
   Returns are -2 for timeout, -1 for nonhamiltonian, otherwise
   number of edges not present in any cycle */
{
	cubgraph done;
	int cyc[MAXN];
        int vv[2],vi[2];
        int nbad,ch;
	int i,j,ii,jj,x,y,z;
 
	nbad = 0;
 
	for (i = 0; i < n; ++i) done[i][0] = done[i][1] = done[i][2] = 0;

	ch = isham(cub,n,ne,weight,vv,vi,0,NULL,NULL,0,cyc);
	if (ch == HABORT) return -2;
	if (ch == NO) return -1;

	for (i = 0; i < n; ++i)
	{
	    x = cyc[i];
            y = cyc[i==n-1?0:i+1];
	    z = cyc[i==0?n-1:i-1];
	    for (j = 0; j < 3; ++j)
		if (cub[x][j] >= 0 &&
		    cub[x][j] != y && cub[x][j] != z) break;
	    if (j < 3) done[x][j] = 1;
	}

	for (ii = 0; ii < n; ++ii)
	for (jj = 0; jj < 3; ++jj)
	{
	    if (cub[ii][jj] > ii && !done[ii][jj])
	    {
		vv[0] = ii;
		vi[0] = jj+1;
		ch = isham(cub,n,ne,weight,vv,vi,1,NULL,NULL,0,cyc);
                if (ch == HABORT) return -2;
		if (ch == NO)
		{
		    x0[nbad] = ii;
		    x1[nbad] = cub[ii][jj];
		    ++nbad;
		}
		else
		{
		    for (i = 0; i < n; ++i)
		    {
	    	        x = cyc[i];
            	        y = cyc[i==n-1?0:i+1];
	    	        z = cyc[i==0?n-1:i-1];
	    	        for (j = 0; j < 3; ++j)
			    if (cub[x][j] >= 0 &&
		    	    cub[x][j] != y && cub[x][j] != z) break;
	    	        if (j < 3) done[x][j] = 1;
		    }

		}
	    }
	}

        return nbad;
}

/**************************************************************************/

static boolean
biconnected_cub(cubgraph cub, int n)
/* test whether cub is biconnected */
{
        int i,sp,v,w,x;
        set visited[MAXM];
        int numvis,num[MAXN],lp[MAXN],stack[MAXN];
	int m,*gv;

        if (n <= 2) return FALSE;

	m = (n + WORDSIZE - 1) / WORDSIZE;
	EMPTYSET(visited,m);
	ADDELEMENT(visited,0);

        stack[0] = 0;
        num[0] = 0;
        lp[0] = 0;
        numvis = 1;
        sp = 0;
        v = 0;
	gv = (int*)cub[v];

        for (;;)
        {
	    for (i = 0; i < 3; ++i)
		if (gv[i] >= 0 && !ISELEMENT(visited,gv[i])) break;

            if (i < 3)
            {
                w = v;
                v = gv[i];  /* visit next child */
                stack[++sp] = v;
		gv = (int*)cub[v];
                ADDELEMENT(visited,v);
                lp[v] = num[v] = numvis++;

                for (i = 0; i < 3; ++i)
		{
		    x = gv[i];
		    if (x >= 0 && x != w && ISELEMENT(visited,x))
			if (num[x] < lp[v])  lp[v] = num[x];
		}
            }
            else
            {
                w = v;                  /* back up to parent */
                if (sp <= 1)          return numvis == n;
                v = stack[--sp];
		gv = (int*)cub[v];
                if (lp[w] >= num[v])  return FALSE;
                if (lp[w] < lp[v])    lp[v] = lp[w];
            }
        }
}

/**************************************************************************/

static boolean
biconnected_cub_v(cubgraph cub, int vv, int n)
/* test whether cub-vv is biconnected */
{
        int i,sp,v,w,x,start;
        set visited[MAXM];
        int numvis,num[MAXN],lp[MAXN],stack[MAXN];
	int m,*gv;

        if (n <= 3) return FALSE;
	start = (vv == 0 ? 1 : 0);

	m = (n + WORDSIZE - 1) / WORDSIZE;
	EMPTYSET(visited,m);
	ADDELEMENT(visited,start);
	ADDELEMENT(visited,vv);

        stack[0] = start;
        num[start] = 0;
        lp[start] = 0;
        numvis = 1;
        sp = 0;
        v = start;
	gv = (int*)cub[v];

        for (;;)
        {
	    for (i = 0; i < 3; ++i)
		if (gv[i] >= 0 && !ISELEMENT(visited,gv[i])) break;

            if (i < 3)
            {
                w = v;
                v = gv[i];  /* visit next child */
                stack[++sp] = v;
		gv = (int*)cub[v];
                ADDELEMENT(visited,v);
                lp[v] = num[v] = numvis++;

                for (i = 0; i < 3; ++i)
		{
		    x = gv[i];
		    if (x >= 0 && x != w && x != vv && ISELEMENT(visited,x))
			if (num[x] < lp[v])  lp[v] = num[x];
		}
            }
            else
            {
                w = v;                  /* back up to parent */
                if (sp <= 1)          return numvis == n-1;
                v = stack[--sp];
		gv = (int*)cub[v];
                if (lp[w] >= num[v])  return FALSE;
                if (lp[w] < lp[v])    lp[v] = lp[w];
            }
        }
}

/**************************************************************************/

static boolean
biconnected_v(graph *g, int vv, int m, int n)
/* test whether g-vv is biconnected */
/* version for arbitrary sizes */
{
        int i,sp,v,w;
	setword ww;
        set sw[MAXM],visited[MAXM];
        int numvis,num[MAXN],lp[MAXN],stack[MAXN];
	int start;
	set *gv;

        if (n <= 3) return FALSE;

	start = vv == 0 ? 1 : 0;
	EMPTYSET(visited,m);
	ADDELEMENT(visited,start);
	ADDELEMENT(visited,vv);

        stack[0] = start;
        num[start] = 0;
        lp[start] = 0;
        numvis = 1;
        sp = 0;
        v = start;
	gv = (set*)g + m*v;

        for (;;)
        {
	    for (i = 0; i < m; ++i)
		if ((ww = gv[i] & ~visited[i])) break; /* = */

            if (i < m)
            {
                w = v;
                v = TIMESWORDSIZE(i) + FIRSTBIT(ww);   /* visit next child */
                stack[++sp] = v;
		gv = (set*)g + m*v;
                ADDELEMENT(visited,v);
                lp[v] = num[v] = numvis++;
		for (i = 0; i < m; ++i)
		    sw[i] = gv[i] & visited[i];
		DELELEMENT(sw,w);
		DELELEMENT(sw,vv);

                for (i = 0; i < m; ++i)
		{
		    ww = sw[i];
                    while (ww)
                    {
			TAKEBIT(w,ww);
			w += TIMESWORDSIZE(i);
                        if (num[w] < lp[v])  lp[v] = num[w];
                    }
		}
            }
            else
            {
                w = v;                  /* back up to parent */
                if (sp <= 1)     return numvis == n-1;
                v = stack[--sp];
		gv = (set*)g + m*v;
                if (lp[w] >= num[v])  return FALSE;
                if (lp[w] < lp[v])    lp[v] = lp[w];
            }
        }
}

/**************************************************************************/

static boolean
triconnected_cub(cubgraph cub, int n)
/* test whether triconnected; awfully inefficient */
{
	int vv;

	if (n <= 3) return FALSE;

	for (vv = n-1; --vv >= 0;)
	    if (!biconnected_cub_v(cub,vv,n)) return FALSE;

	return TRUE;
}

/**************************************************************************/

static boolean
gtocub(graph *g, int m, int n, cubgraph cub, int *pne)
/* Convert nauty-format graph into cubgraph.  Returns FALSE if there
 * are any vertices of degree 4 or more.
 */
{
	int i,j,nde;
	set *gi;

	nde = 0;
	for (i = 0, gi = (set*)g; i < n; ++i, gi += m)
	{
	    cub[i][0] = cub[i][1] = cub[i][2] = -1;
	    j = nextelement(gi,m,-1);
	    if (j < 0) continue;
	    cub[i][0] = j;
	    ++nde;
	    j = nextelement(gi,m,j);
            if (j < 0) continue;
            cub[i][1] = j;
	    ++nde;
	    j = nextelement(gi,m,j);
            if (j < 0) continue;
            cub[i][2] = j;
	    ++nde;
	    if (nextelement(gi,m,j) >= 0) return FALSE;
	}

	*pne = nde / 2;
	return TRUE;
}

/**************************************************************************/

static boolean
sgtocub(sparsegraph *sg, cubgraph cub, int *pne)
/* Convert sparse-format graph into cubgraph.  Returns FALSE if there
 * are any vertices of degree 4 or more.
 */
{
	int *d,*e;
	size_t *v;
	int n,i,j,vi;
	unsigned long int nde;

	nde = 0;
	n = sg->nv;
	SG_VDE(sg,v,d,e);
	for (i = 0; i < n; ++i)
	{
	    if (d[i] >= 4) return FALSE;
	    vi = v[i];

	    cub[i][0] = cub[i][1] = cub[i][2] = -1;
	    for (j = 0; j < d[i]; ++j)
	       cub[i][j] = e[vi+j];
	    nde += d[i];
	}

	*pne = nde / 2;
	return TRUE;
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
        char *infilename,*outfilename;
        FILE *infile,*outfile,*msgfile;
        boolean badargs,biconn,triconn,fragment;
	boolean in,out,inin,inout,e34,e34plus,testing;
        nauty_counter numread,noncub,nonham,nonconn,numto;
        int ch,weight,m,n,ne,i,namesgot,nl;
	int nbad,limit,x0[BADLIM],x1[BADLIM],y0[BADLIM],y1[BADLIM];
	sparsegraph sg;
	int vv[MAXA],vi[MAXA],nvv,cyc[MAXN];
        int yy[MAXA],yi[MAXA],nyy;
	double t0,t1;
	cubgraph cub;
        char *arg;
	int codetype;

        infilename = outfilename = NULL;
        badargs = FALSE;
	e34plus = e34 = in = out = inin = inout = FALSE;
	fragment = biconn = triconn = testing = FALSE;
	verbose = 0;
	weight = 100;
	nvv = nyy = 0;
	timeout = 0;
	repeats = 1;

     /* parse argument list */

        namesgot = 0;
        for (i = 1; i < argc && !badargs; ++i)
        {
            arg = argv[i];
            if (arg[0] == '-' && arg[1] != '\0')
            {
                if (arg[1] == 'v')
                    verbose = (verbose == 0 ? 1 : verbose);
		else if (arg[1] == 'V')
		    verbose = 2;
                else if (arg[1] == 'i')
                    in = TRUE;
                else if (arg[1] == 'I')
                    inin = TRUE;
                else if (arg[1] == 'o')
                    out = TRUE;
		else if (arg[1] == 'x')
		    inout = TRUE;
		else if (arg[1] == 'e')
		    e34 = TRUE;
		else if (arg[1] == 'E')
		    e34 = e34plus = TRUE;
                else if (arg[1] == 'b')
                    biconn = TRUE;
                else if (arg[1] == 't')
                    triconn = TRUE;
                else if (arg[1] == 'Q')
                    testing = TRUE;
                else if (arg[1] == 'F')
                    fragment = TRUE;
		else if (arg[1] == 'n')
		{
		    if (nvv == MAXA)
		    {
			fprintf(stderr,">E cubhamg : too many -n switches\n");
			exit(1);
		    }
		    if (sscanf(arg+2,"%d-%d",&vv[nvv],&vi[nvv]) != 2)
			badargs = TRUE;
		    else if (vi[nvv] < 1 || vi[nvv] > 3)
		    {
			fprintf(stderr,
                               ">E cubhamg : second arg of -n must be 1..3\n");
			exit(1);
		    }
		    ++nvv;
		}
		else if (arg[1] == 'y')
		{
		    if (nyy == MAXA)
		    {
			fprintf(stderr,">E cubhamg : too many -y switches\n");
			exit(1);
		    }
		    if (sscanf(arg+2,"%d-%d",&yy[nyy],&yi[nyy]) != 2)
			badargs = TRUE;
		    else if (yi[nyy] < 1 || yi[nyy] > 3)
		    {
			fprintf(stderr,
                               ">E cubhamg : second arg of -y must be 1..3\n");
			exit(1);
		    }
		    ++nyy;
		}
		else if (arg[1] >= '0' && arg[1] <= '9')
		    sscanf(arg+1,"%d",&weight);
		else if (arg[1] == 'T')
		{
		    if (sscanf(arg+2,"%ld",&timeout) != 1 || timeout < 0)
		    {
			fprintf(stderr,">E cubhamg : bad -T value\n");
			exit(1);
		    }
		    else if (timeout > NO_LIMIT-1)
			timeout = NO_LIMIT - 1;
		}
		else if (arg[1] == 'R')
		{
		    if (sscanf(arg+2,"%ld",&repeats) != 1 || repeats < 1)
		    {
			fprintf(stderr,">E cubhamg : bad -R value\n");
			exit(1);
		    }
		}
		else
		    badargs = TRUE;
            }
            else
            {
                if (namesgot == 0)
                {
                    namesgot = 1;
                    infilename = arg;
                }
                else if (namesgot == 1)
                {
                    namesgot = 2;
                    outfilename = arg;
                }
                else
                    badargs = TRUE;
            }
        }

	if (badargs)
	{
	    fprintf(stderr,
         ">E Usage: cubhamg [-#] [-v | -V] [-n#-#] [-y#-#] [infile [outfile]]\n");
	    exit(1);
	}

     /* open input file */

        if (infilename && infilename[0] == '-') infilename = NULL;
        infile = opengraphfile(infilename,&codetype,FALSE,1);
        if (!infile) exit(1);
        if (!infilename) infilename = "stdin";
 
        if (infilename == NULL) infilename = "stdin";

        NODIGRAPHSYET(codetype);

     /* open output file */

        if (outfilename == NULL || outfilename[0] == '-')
        {
            outfile = stdout;
            outfilename = "stdout";
            msgfile = stderr;
        }
        else
        {
            msgfile = stdout;
            if ((outfile = fopen(outfilename,"w")) == NULL)
            {
                fprintf(stderr,
                    ">E cubhamg: can't open %s for writing\n",outfilename);
                ABORT(">E cubhamg");
            }
        }

	if (triconn) biconn = FALSE;
	if (e34) in = out = inin = inout = FALSE;
	if (inout) in = out = inin = FALSE;
	if (inin) in = out = FALSE;
	if (out) in = FALSE;
	if (testing) out = TRUE;

	if (codetype&HAS_HEADER)
        {
            if (codetype&SPARSE6) writeline(outfile,SPARSE6_HEADER);
            else                  writeline(outfile,GRAPH6_HEADER);
        }

//if (codetype&SPARSE6) outproc = writes6;
//else                  outproc = writeg6;

	numread = noncub = nonham = nonconn = 0;
	numto = 0;
	INITRANBYTIME;
	t0 = CPUTIME;

	limit = verbose ? BADLIM : 1;

	SG_INIT(sg);
	while (read_sg(infile,&sg))
	{
	    ++numread;
	    n = sg.nv;

	    if (n >= MAXN-1)
	    {
		fprintf(stderr,
                    ">E cubhamg: input " COUNTER_FMT " too big\n",numread);
		exit(1);
	    }
#if MAXM==1
	    if (3*n >= 2*WORDSIZE)
	    {
		fprintf(stderr,
                    ">E cubhamg: must compile with MAXM>1 for ne>WORDSIZE\n");
		exit(1);
	    }
#endif

	    if (!sgtocub(&sg,cub,&ne))
	    {
		if (verbose)
		    fprintf(msgfile,"Input " COUNTER_FMT
                                    " has maxdeg>3.\n",numread);
		++noncub;
	    }
	    else if (biconn && !biconnected_cub(cub,n))
		++nonconn;
	    else if (triconn && !triconnected_cub(cub,n))
		++nonconn;
            else if (e34)
            {
                if ((nbad = hase34(cub,n,ne,x0,x1,y0,e34plus,limit)) > 0)
                {
                    if (verbose)
		    {
                        fprintf(msgfile,"Input " COUNTER_FMT
                                        " fails property e34%s:",
                                               numread,e34plus ? "+" : "");
			if (e34plus)
			    for (i = 0; i < nbad; ++i)
			        fprintf(msgfile," %d-%d[%02x]",
						x0[i],x1[i],y0[i]);
			else
			    for (i = 0; i < nbad; ++i)
			        fprintf(msgfile," %d-%d",x0[i],x1[i]);
			fprintf(msgfile,"\n");
		    }
                    ++nonham;
                    writelast(outfile);
                }
            }
	    else if (inout)
	    {
		if ((nbad = hasinout(cub,n,ne,x0,x1,y0,y1,limit)) > 0)
		{
		    if (verbose)
                    { 
                        fprintf(msgfile,"Input " COUNTER_FMT
                                        " fails property -+:",numread);
                        for (i = 0; i < nbad; ++i) 
                            fprintf(msgfile," %d-%d/%d-%d",
				    x0[i],x1[i],y0[i],y1[i]); 
                        fprintf(msgfile,"\n");
                    }
		    ++nonham;
                    writelast(outfile);
		}
	    }
            else if (inin)
            {
                if ((nbad = hasinin(cub,n,ne,x0,x1,y0,y1,limit)) > 0)
                {
                    if (verbose)
		    {  
                        fprintf(msgfile,"Input " COUNTER_FMT
                                        " fails property ++:",numread); 
                        for (i = 0; i < nbad; ++i) 
                            fprintf(msgfile," %d-%d/%d-%d",
                                    x0[i],x1[i],y0[i],y1[i]);
                        fprintf(msgfile,"\n"); 
                    }
                    ++nonham;
                    writelast(outfile);
                }
            }
            else if (out)
            { 
                if ((nbad = hasout(cub,n,ne,weight,x0,x1,limit)) != 0)
                { 
                    if (verbose) 
		    {
			if (nbad == -2)
                            fprintf(msgfile,"Input " COUNTER_FMT
                                        " timed out\n",numread);
			else if (nbad == -1)
                            fprintf(msgfile,"Input " COUNTER_FMT
                                        " is nonhamiltonian\n",numread);
			else
			{
                            fprintf(msgfile,"Input " COUNTER_FMT
                                            " fails property -:",numread);
                            for (i = 0; i < nbad; ++i)
                                fprintf(msgfile," %d-%d",x0[i],x1[i]);
                            fprintf(msgfile,"\n");
		 	}
                    }
		    if (nbad != -1)
		    {
                        ++nonham; 
                        writelast(outfile);
		    }
                } 
            }
            else if (in)
            {
                if ((nbad = hasin(cub,n,ne,x0,x1,limit)) > 0)
                {
                    if (verbose)
                    {    
                        fprintf(msgfile,"Input " COUNTER_FMT
                                        " fails property +:",numread); 
                        for (i = 0; i < nbad; ++i) 
                            fprintf(msgfile," %d-%d",x0[i],x1[i]); 
                        fprintf(msgfile,"\n"); 
                    }
                    ++nonham;
                    writelast(outfile);
                }
            }
	    else if (fragment)
            {
		dofragment(numread,cub,n,ne,weight);
            }
	    else if ((ch = isham(cub,n,ne,weight,vv,vi,nvv,yy,yi,nyy,cyc)) == NO)
	    {
	   	if (verbose)
                    fprintf(msgfile,"Input " COUNTER_FMT
                                    " is not hamiltonian.\n",numread);
                ++nonham;
		writelast(outfile);
	    }
	    else if (ch == HABORT)
	    {
                if (verbose)
                    fprintf(msgfile,"Input " COUNTER_FMT " timed out.\n",numread);
                ++numto;
                 writelast(outfile);
            }
	    else if (verbose >= 2)
	    {
		fprintf(msgfile,"Cycle in input " COUNTER_FMT ":\n",numread);
		if (n <= 100) nl = 26;
		else          nl = 19;
		for (i = 0; i < n; ++i)
	 	{
		    if (i > 0 && i % nl == 0) fprintf(msgfile,"\n ");
		    fprintf(msgfile," %d",cyc[i]);
		}
		fprintf(msgfile,"\n");
	    }
	}
	t1 = CPUTIME;

	fprintf(msgfile,">C " COUNTER_FMT
                        " graphs read from %s\n",numread,infilename);
	if (noncub > 0)
	    fprintf(msgfile,">C " COUNTER_FMT " graphs with maxdeg > 3\n",noncub);
	if (nonconn > 0)
	    fprintf(msgfile,">C " COUNTER_FMT " graphs not %sconnected\n",
			    nonconn,biconn ? "bi" : "tri");
	if (numto > 0)
	    fprintf(msgfile,">C " COUNTER_FMT " graphs timed out\n",numto);
	fprintf(msgfile,
 	    ">Z " COUNTER_FMT " %s graphs written to %s; %3.2f sec\n",
              nonham+numto, e34 ? (e34plus ? "non-e34+" : "non-e34") :
	            inout ? "non-inout" : 
                        out ? "non-out" :
		            in ? "non-in" : "nonhamiltonian",
                outfilename,t1-t0);
	if (verbose)
	{
	    fprintf(msgfile,"Tries:"); 
	    for (i = 0; i <= NUMMAXNODES && numtries[i] > 0; ++i)
		fprintf(msgfile," " COUNTER_FMT,numtries[i]);
	    fprintf(msgfile,"\n");
	}
#if MAXES
	fprintf(msgfile,"Maximum level = %d; Maximum classstack = %d\n",
		maxlevel,maxclassstack);
#endif

	return 0;
}
