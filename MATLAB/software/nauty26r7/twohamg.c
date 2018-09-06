/* twohamg.c  */

#define USAGE "twohamg [-sgvq] [-L#] [infile [outfile]]"

#define HELPTEXT \
" Partition quartic graphs into two hamiltonian cycles.\n\
  Output those which cannot be partitioned.\n\
\n\
    -s  force output to sparse6 format\n\
    -g  force output to graph6 format\n\
        If neither -s or -g are given, the output format is\n\
        determined by the header or, if there is none, by the\n\
        format of the first input graph. Also see -S. \n\
\n\
    The output file will have a header if and only if the input file does.\n\
\n\
    -p  Read a cubic graph and use its prism. Vertex i of the input becomes\n\
           vertices 2*i,2*i+1 in the prism.\n\
    -x  Test for decompositions using each 2-path\n\
    -X  As -x but only output if two 2-paths are missed at some vertex\n\
    -y  Test for decompositions using each non-triangular 3-path\n\
    -t# With -x and -X, consider only paths with center #\n\
        With -y, consider only paths starting at #\n\
    -Y  With -p, only consider paths whose central edge is vertical\n\
    -v  Give a partition for those graphs who have one and a message\n\
        for those which don't. With -x, list exceptional 2-paths.\n\
    -L# Limit to 1000*# iterations; write with message if timeout.\n\
        Graphs that time out are written to the output.\n\
\n\
    -q  suppress auxiliary information\n"

#define DEBUG 0
#define FISHTAIL 1   /* Whether to use fish-tail rule */

/*************************************************************************/

#include "gtools.h"

#define WHITE 0  /* unknown colour */
#define BLUE 1
#define RED 2

#define NO 0
#define YES 1
#define TIMEOUT 2
#define NOTHING 3

typedef struct addrval_struct
{
    int *addr;
    int value;
} addrval;

static long seed = 314159265;

/* valstack is a stack to save and restore degrees, farends, and colours */
DYNALLSTAT(addrval,valstack,valstack_sz);
static addrval *valstacktop;
#define SETVAL(ad,val) { valstacktop->addr = &(ad); \
     valstacktop->value = ad; ++valstacktop; ad = val; }
#define UNSETVAL { --valstacktop; *(valstacktop->addr) = valstacktop->value; }

typedef struct p4
{
    int e1,e2,e3;
    int v1,v2,v3,v4;
    boolean ok;
} p4;

DYNALLSTAT(int,bluefarend,bluefarend_sz);   /* Parallel to sg.v */
DYNALLSTAT(int,redfarend,redfarend_sz);   /* Parallel to sg.v */
DYNALLSTAT(int,eno,eno_sz);         /* Edge numbers, parallel to sg.e */
DYNALLSTAT(int,cross,cross_sz);         /* Parallel to sg.eno */
DYNALLSTAT(int,colour,colour_sz);   /* indexed by edge number */
DYNALLSTAT(int,v1,v1_sz); 
DYNALLSTAT(int,v2,v2_sz); 
DYNALLSTAT(int,bluedeg,bluedeg_sz); /* Parallel to sg.v */
DYNALLSTAT(int,reddeg,reddeg_sz);   /* Parallel to sg.v */
DYNALLSTAT(int,beste,beste_sz);  /* List of possible best edges */
DYNALLSTAT(p4,p4list,p4list_sz);  /* List of non-triangular p4s */

/* vstack is a stack of interesting vertices; onstack says whether a
   vertex is on vstack so that we don't put it there twice */
DYNALLSTAT(int,vstack,vstack_sz);
DYNALLSTAT(boolean,onstack,onstack_sz);
static int *top;
#define PUSH(v) if (!onstack[v]) { *(top++) = (v); onstack[v] = TRUE; }
#define POP(v) { v = *(--top); onstack[v] = FALSE; }
#define STACKISEMPTY (top == vstack)

static nauty_counter nodes,limit,totallimit;
static nauty_counter locallimit[]
  = {5,7,10,20,30,40,50,60,70,100,200,400,600,1000,2000,4000,8000,
     20000,50000,100000,200000,5000000,0};
#define NUMLIMITS (sizeof(locallimit)/sizeof(locallimit[0]))
static nauty_counter nphases[100];  /* Must be greater than NUMLIMITS */
static nauty_counter ntimeouts;

/**************************************************************************/

static void
makeprism_sg(sparsegraph *sg, sparsegraph *sh)
/* Make the cartesian product sh = sg x K2.
   sh must be initialized. */
{
    size_t *v,*vv,vi,vvi;
    int *d,*dd,*e,*ee;
    int i,j,n;

    SG_VDE(sg,v,d,e);
    n = sg->nv;

    SG_ALLOC(*sh,2*sg->nv,2*(n+sg->nde),"makeprism_sg");
    sh->nv = 2*n;
    sh->nde = 2*(n+sg->nde);
    SG_VDE(sh,vv,dd,ee);
 
    vvi = 0;
    for (i = 0; i < n; ++i)
    {
	vi = v[i];
	dd[2*i] = dd[2*i+1] = d[i] + 1;
	vv[2*i] = vvi;
	for (j = 0; j < d[i]; ++j) ee[vvi++] = 2*e[vi++];
	ee[vvi++] = 2*i+1;
	vi = v[i];
	vv[2*i+1] = vvi;
	for (j = 0; j < d[i]; ++j) ee[vvi++] = 2*e[vi++]+1;
	ee[vvi++] = 2*i;
    }
}

/**************************************************************************/

static void
dumpdata(int id, int nblue, int nred, int n)
{
    int i,j;

    printf("%d: nblue=%d nred=%d -------------------------\n",id,nblue,nred);

#if DEBUG > 1
    for (i = 0; i < n; ++i)
    {
        printf("%2d: ",i);
	for (j = 0; j < 4; ++j)
	    printf(" %2d%c",eno[4*i+j],"wbr"[colour[eno[4*i+j]]]);
        printf("   b=%d r=%d bfe=%d rfe=%d\n",
	    bluedeg[i],reddeg[i],bluefarend[i],redfarend[i]);
    }
#endif
}

/**************************************************************************/

static void
initialise_g(int n, int *e)
/* Allocate all and initialise eno[],v1[],v2[]
   e is a vector like sg.e */
{
    int ne,i,j,k,l;
 
    DYNALLOC1(addrval,valstack,valstack_sz,10*n,"malloc");
    DYNALLOC1(int,bluefarend,bluefarend_sz,n,"malloc");
    DYNALLOC1(int,redfarend,redfarend_sz,n,"malloc");
    DYNALLOC1(int,eno,eno_sz,4*n,"malloc"); 
    DYNALLOC1(int,v1,v1_sz,2*n,"malloc");
    DYNALLOC1(int,v2,v2_sz,2*n,"malloc");
    DYNALLOC1(int,colour,colour_sz,2*n,"malloc");
    DYNALLOC1(int,bluedeg,bluedeg_sz,n,"malloc");
    DYNALLOC1(int,reddeg,reddeg_sz,n,"malloc");
    DYNALLOC1(int,vstack,vstack_sz,n,"malloc");
    DYNALLOC1(boolean,onstack,onstack_sz,n,"malloc");
    DYNALLOC1(int,beste,beste_sz,2*n,"malloc");

  /* Randomize e; seems to be no purpose for this any more. */

    for (i = 0; i < n; ++i)
    {
	j = KRAN(4);
	k = e[4*i+j]; e[4*i+j] = e[4*i+3]; e[4*i+3] = k; 
	j = KRAN(3);
	k = e[4*i+j]; e[4*i+j] = e[4*i+2]; e[4*i+2] = k;
	j = KRAN(2);
	k = e[4*i+j]; e[4*i+j] = e[4*i+1]; e[4*i+1] = k;
    }

    ne = 0;
    for (i = 0; i < n; ++i)
    {
	for (j = 0; j < 4; ++j)
	{
	    k = e[4*i+j];
	    if (k > i)
	    {
		v1[ne] = i;
		v2[ne] = k;
		eno[4*i+j] = ne++;
            }
	    else    /* Note: this code assumes a simple graph */
	    {
		for (l = 0; l < 4; ++l)
		    if (e[4*k+l] == i) break;
		eno[4*i+j] = eno[4*k+l];
	    }
	}
    }
    if (ne != 2*n) gt_abort(">E ne is incorrect");

#if DEBUG
    { int ii;
        printf("===== n=%d === ne=%d ===================\n",n,ne);
    
        for (ii = 0; ii < n; ++ii)
	    printf("%2d: %2d %2d %2d %2d    %2d %2d %2d %2d\n",
	       ii,e[4*ii],e[4*ii+1],e[4*ii+2],e[4*ii+3],
                  eno[4*ii],eno[4*ii+1],eno[4*ii+2],eno[4*ii+3]);
    }
#endif
}

/**************************************************************************/

static void
initialise_colouring(int n)
/* Make all edges white */
{
    int i;

    for (i = 0; i < 2*n; ++i) colour[i] = WHITE;

    for (i = 0; i < n; ++i)
    {
        bluefarend[i] = i;
        redfarend[i] = i;
        bluedeg[i] = 0;
        reddeg[i] = 0;
        onstack[i] = FALSE;
    }    

    top = vstack;
    valstacktop = valstack;
}

/**************************************************************************/

static boolean
makeblue(int edge, boolean lastok)
/* Colour WHITE edge BLUE, return success. 
   lastok indicates if it is ok to add the final blue edge */
{
    int w1,w2,f1,f2;

    w1 = v1[edge];
    w2 = v2[edge];
    if (bluedeg[w1] == 2 || bluedeg[w2] == 2) return FALSE;

    f1 = bluefarend[w1];
    f2 = bluefarend[w2];
    if (f1 == w2 && !lastok) return FALSE;

    SETVAL(colour[edge],BLUE);
    SETVAL(bluedeg[w1],bluedeg[w1]+1);
    SETVAL(bluedeg[w2],bluedeg[w2]+1);
    SETVAL(bluefarend[f1],f2);
    SETVAL(bluefarend[f2],f1);
    PUSH(w1);
    PUSH(w2);
    PUSH(f1);
    // PUSH(f2);   not sure this one is needed

    return TRUE;
}

/**************************************************************************/

static boolean
makered(int edge, boolean lastok)
/* Colour WHITE edge RED, return success. 
   lastok indicates if it is ok to add the final red edge */
{
    int w1,w2,f1,f2;

    w1 = v1[edge];
    w2 = v2[edge];
    if (reddeg[w1] == 2 || reddeg[w2] == 2) return FALSE;

    f1 = redfarend[w1];
    f2 = redfarend[w2];
    if (f1 == w2 && !lastok) return FALSE;

    SETVAL(colour[edge],RED);
    SETVAL(reddeg[w1],reddeg[w1]+1);
    SETVAL(reddeg[w2],reddeg[w2]+1);
    SETVAL(redfarend[f1],f2);
    SETVAL(redfarend[f2],f1);
    PUSH(w1);
    PUSH(w2);
    PUSH(f1);
    // PUSH(f2);   not sure this one is needed

    return TRUE;
}

/**************************************************************************/

static boolean
propagate(int n, int *e, int *nblue, int *nred)
/* Look at active vertices and propagate colourings */
{
    int i,j,v,w;

    while (!STACKISEMPTY)
    {
	POP(v);

        if (reddeg[v] == 2 && bluedeg[v] < 2)
        {
	    for (i = 0; i < 4; ++i)
	    {
		j = eno[4*v+i];
		if (colour[j] == WHITE)
		{
		    if (!makeblue(j,*nblue==n-1)) return FALSE;
		    ++*nblue;
		}
	    }
        }
	else if (bluedeg[v] == 2 && reddeg[v] < 2)
        {
	    for (i = 0; i < 4; ++i)
	    {
		j = eno[4*v+i];
		if (colour[j] == WHITE)
		{
		    if (!makered(j,*nred==n-1)) return FALSE;
		    ++*nred;
		}
	    }
        }

	if (bluedeg[v] == 1)
	{
	    w = bluefarend[v];
	    for (i = 0; i < 4; ++i)
	        if (e[4*v+i] == w) break;
            if (i < 4)
	    {
		j = eno[4*v+i];
		if (colour[j] == WHITE)
		{
		    if (*nblue == n-1)
		    {
			if (!makeblue(j,TRUE)) return FALSE;
			++*nblue;
		    }
		    else
		    {
		        if (!makered(j,*nred==n-1)) return FALSE;
		        ++*nred;
		    }
		}
	    }
	}

	if (reddeg[v] == 1)
	{
	    w = redfarend[v];
	    for (i = 0; i < 4; ++i)
	        if (e[4*v+i] == w) break;
            if (i < 4)
	    {
		j = eno[4*v+i];
		if (colour[j] == WHITE)
		{
                    if (*nred == n-1)
                    {
                        if (!makered(j,TRUE)) return FALSE;
                        ++*nred;
                    }    
                    else
                    {
                        if (!makeblue(j,*nblue==n-1)) return FALSE;
                        ++*nblue;
                    }
		}
	    }
        }
    }

    return TRUE;
}

/**************************************************************************/

static int
fishtail(int n, int *nblue, int *nred)
/* Try to colour edges by the fishtail method
   Return NO (contradiction), YES (did some good), or NOTHING.  */
{
    int i,j;
    int e1,w1,e2,w2,e3,w3;
    int edge,col;
    int ans;

    ans = NOTHING;

    for (i = 0; i < n; ++i)
    {
	edge = -1;

	if (reddeg[i] == 1 && bluedeg[i] == 0 && *nblue < n-2)
	{
	    j = 4*i;
	    if (colour[eno[j]] != WHITE) ++j;
	    e1 = eno[j];
	    w1 = (v1[e1] == i ? v2[e1] : v1[e1]);
	    ++j;
	    if (colour[eno[j]] != WHITE) ++j;
	    e2 = eno[j];
	    w2 = (v1[e2] == i ? v2[e2] : v1[e2]);
	    ++j;
	    if (colour[eno[j]] != WHITE) ++j;
	    e3 = eno[j];
	    w3 = (v1[e3] == i ? v2[e3] : v1[e3]);

	    if (bluedeg[w1] == 1 && bluefarend[w1] == w2)
	    {
		edge = e3;
                col = BLUE;
	    }
	    else if (bluedeg[w1] == 1 && bluefarend[w1] == w3)
	    {
		edge = e2;
                col = BLUE;
	    }
	    else if (bluedeg[w2] == 1 && bluefarend[w2] == w3)
	    {
		edge = e1;
                col = BLUE;
	    }
	}
	else if (reddeg[i] == 0 && bluedeg[i] == 1 && *nred < n-2)
	{
	    j = 4*i;
	    if (colour[eno[j]] != WHITE) ++j;
	    e1 = eno[j];
	    w1 = (v1[e1] == i ? v2[e1] : v1[e1]);
	    ++j;
	    if (colour[eno[j]] != WHITE) ++j;
	    e2 = eno[j];
	    w2 = (v1[e2] == i ? v2[e2] : v1[e2]);
	    ++j;
	    if (colour[eno[j]] != WHITE) ++j;
	    e3 = eno[j];
	    w3 = (v1[e3] == i ? v2[e3] : v1[e3]);

	    if (reddeg[w1] == 1 && redfarend[w1] == w2)
	    {
		edge = e3;
                col = RED;
	    }
	    else if (reddeg[w1] == 1 && redfarend[w1] == w3)
	    {
		edge = e2;
                col = RED;
	    }
	    else if (reddeg[w2] == 1 && redfarend[w2] == w3)
	    {
		edge = e1;
                col = RED;
	    }
	}
	if (edge >= 0)
	{
	    if (col == BLUE)
	    {
		if (!makeblue(edge,*nblue == n-1)) return NO;
		++*nblue;
	    }
	    else
	    {
		if (!makered(edge,*nred == n-1)) return NO;
		++*nred;
	    }
	    ans = YES;
#if DEBUG
            printf("Fishtail: %d=(%d-%d) -> %c\n",
                               edge,v1[edge],v2[edge],"wbr"[col]);
#endif
        }
    }

    return ans;
}

/**************************************************************************/

static int
searchnode(int level, int n, int *e, int nblue, int nred)
{
    boolean ok;
    int i,status,nbest;
    addrval *valptr;
    int best,score,bestscore;
    int fe,fc;
    long ran;
     
    ok = propagate(n,e,&nblue,&nred);

#if DEBUG
    if (ok) dumpdata(level,nblue,nred,n);
    else { printf("FAIL\n"); return NO; }
#else
    if (!ok) return NO;
#endif

    if (nblue == n && nred == n) return YES;

#if FISHTAIL
    status = fishtail(n,&nblue,&nred);
    if (status == NO) return NO;
    if (status == YES) return searchnode(level+1,n,e,nblue,nred);
#endif

    valptr = valstacktop;

    bestscore = -1;
    nbest = 0;
    for (i = 0; i < 2*n; ++i)
	if (colour[i] == WHITE)
        {
	    score = (bluedeg[v1[i]] == 1) + (bluedeg[v2[i]] == 1) +
                    (reddeg[v1[i]] == 1) + (reddeg[v2[i]] == 1);
	    if (score > bestscore)
	    {
		bestscore = score;
		beste[0] = i;
		nbest = 1;
	    }
	    else if (score == bestscore)
	        beste[nbest++] = i;
        }

    if (bestscore == 0 && nred + nblue > 0)
	return FALSE;   /* Disconnected */

    ran = KRAN(2*nbest);
    best = beste[ran/2];

    if ((ran&1))
    {
        if (makeblue(best,nblue==n-1))
        {
#if DEBUG
            printf("   setting %d(%d-%d) blue\n",best,v1[best],v2[best]);
#endif
	    status = searchnode(level+1,n,e,nblue+1,nred);
	    if (status != NO) return status;
	    while (valstacktop > valptr) UNSETVAL;
        }

        if (++nodes == limit) return TIMEOUT;

        if (makered(best,nred==n-1))
        {
#if DEBUG
            printf("   setting %d(%d-%d) red\n",best,v1[best],v2[best]);
#endif
	    status = searchnode(level+1,n,e,nblue,nred+1);
	    if (status != NO) return status;
	    while (valstacktop > valptr) UNSETVAL;
        }
    }
    else
    {
        if (makered(best,nred==n-1))
        {
#if DEBUG
            printf("   setting %d(%d-%d) red\n",best,v1[best],v2[best]);
#endif
	    status = searchnode(level+1,n,e,nblue,nred+1);
	    if (status != NO) return status;
	    while (valstacktop > valptr) UNSETVAL;
        }

        if (++nodes == limit) return TIMEOUT;

        if (makeblue(best,nblue==n-1))
        {
#if DEBUG
            printf("   setting %d(%d-%d) blue\n",best,v1[best],v2[best]);
#endif
	    status = searchnode(level+1,n,e,nblue+1,nred);
	    if (status != NO) return status;
	    while (valstacktop > valptr) UNSETVAL;
        }
    }

    return NO;
}

/**************************************************************************/

static int
dispatchsearch(int n, int *e, int nblue, int nred)
/* Multiple tries at solution */
{
    int i,status;
    addrval *valptr;
    boolean ok;
    nauty_counter remaininglimit;

    ok = propagate(n,e,&nblue,&nred);

#if DEBUG
    if (ok) dumpdata(0,nblue,nred,n);
    else { printf("FAIL\n"); ++nphases[0]; return NO; }
#else
    if (!ok) { ++nphases[0]; return NO; }
#endif

    valptr = valstacktop;

    remaininglimit = totallimit;

    for (i = 0; i < NUMLIMITS; ++i)
    {
	if (totallimit > 0)
	{
	    limit = locallimit[i];
            if (limit == 0 || limit > remaininglimit) limit = remaininglimit;
	    remaininglimit -= limit;
            if (limit == 0) limit = 1;
	}
	else
	    limit = locallimit[i];

	nodes = 0;

	status = searchnode(1,n,e,nblue,nred);
	if (status != TIMEOUT)
        {
	    ++nphases[i+1];
	    return status;
	}
	while (valstacktop > valptr) UNSETVAL;
    }
    
    ++ntimeouts;
    return TIMEOUT;
}

/**************************************************************************/

static int
isdecomposable(sparsegraph sg, int *ham1, int *ham2)
/* Test if sg is decomposable as two hamiltonian cycles */
{
    int nblue,nred,n,e0,status;
    int i0,i,j,lastv,laste;

    n = sg.nv;
    initialise_g(n,sg.e);
    initialise_colouring(n);

    e0 = KRAN(2*n);

    nblue = nred = 0;
    if (!makeblue(e0,nblue==n-1))   /* make edge e0 blue */
	return NO;
    ++nblue;

#if DEBUG
    printf("   setting %d(%d-%d) blue\n",e0,v1[0],v2[0]);
#endif

    // status = searchnode(1,n,sg.e,nblue,nred);
    status = dispatchsearch(n,sg.e,nblue,nred);

    if (status != YES) return status;

    if (!ham1 || !ham2) return YES;

    for (i0 = 0; colour[i0] != BLUE; ++i0) {} 
    ham1[0] = v1[i0];
    laste = i0;
    lastv = ham1[1] = v2[i0];
    for (i = 2; i < n; ++i)
    {
	for (j = 0; j < 4; ++j)
	    if (eno[4*lastv+j] != laste && colour[eno[4*lastv+j]] == BLUE)
		break;
	laste = eno[4*lastv+j];
	lastv = ham1[i] = (v1[laste] == lastv ? v2[laste] : v1[laste]);
    }

    for (i0 = 0; colour[i0] != RED; ++i0) {}
    ham2[0] = v1[i0];
    laste = i0;
    lastv = ham2[1] = v2[i0];
    for (i = 2; i < n; ++i)
    {
	for (j = 0; j < 4; ++j)
	    if (eno[4*lastv+j] != laste && colour[eno[4*lastv+j]] == RED)
		break;
	laste = eno[4*lastv+j];
	lastv = ham2[i] = (v1[laste] == lastv ? v2[laste] : v1[laste]);
    }

    return YES;
}

/**************************************************************************/

static int
iscrossdecomposable(sparsegraph sg, int vertex)
/* Test if sg is decomposable as two hamiltonian cycles
   in all ways through each vertex (or the given vertex if it
   is >= 0.  Details left in cross[].
   Return -1: not decomposable at all
           0: misses some 2-paths including two at some vertices
           1: misses some 2-paths but at most one per vertex
           2: fully decomposable
           3: timeout
*/
{
    int i,j,k,l,n,c;
    int ans,status;
    int imin,imax;

    n = sg.nv;
    DYNALLOC1(int,cross,cross_sz,4*n,"malloc");

    if (vertex >= n) gt_abort(">E vertex given to -t is too large");

    for (i = 0; i < 4*n; ++i) cross[i] = 0;

    status = isdecomposable(sg,FALSE,FALSE);
    if (status == NO) return -1;
    else if (status == TIMEOUT) return 3;

    for (i = 0; i < n; ++i)
    {
	c = colour[eno[4*i]];
	for (j = 1; j < 4; ++j)
	    if (colour[eno[4*i+j]] == c) cross[4*i+j] = 1;
    }

    ans = 2;

    if (vertex >= 0) { imin = imax = vertex; }
    else             { imin = 0; imax = n-1; }

    for (i = imin; i <= imax; ++i)
    for (j = 1; j < 4; ++j)
        if (!cross[4*i+j])
	{
	    initialise_colouring(n);
	    if (makeblue(eno[4*i],FALSE) && makeblue(eno[4*i+j],n==2))
	    {
		status = dispatchsearch(n,sg.e,2,0);
		if (status == TIMEOUT) return 3;
		if (status == YES)
		{
    		    for (k = 0; k < n; ++k)
    		    {
		        c = colour[eno[4*k]];
		        for (l = 1; l < 4; ++l)
	                if (colour[eno[4*k+l]] == c) cross[4*k+l] = 1;
                    }
		}
		else
		    ans = 1;
	    }
	    else
	        ans = 1;
	}

    if (ans == 1)
        for (i = 0; i < n; ++i)
	{
	    if (cross[4*i+1] + cross[4*i+2] + cross[4*i+3] <= 1)
		ans = 0;
	}

    return ans;
}

/**************************************************************************/

static int
p4decomposition(sparsegraph sg, int vertex, boolean vertical)
/* Test which non-triangular P4s extend to a decomposition.
   Return -2: timeout
          -1: not decomposable at all
          >=0: number of P4s missed
   If 0 <= vertex < n, only P4s starting at that vertex are considered.
   If vertical (only for prisms), the central edge must be vertical.
   The paths missed can be found in p4list[].
*/
{
    int i,j,k,l,n,c;
    int ans,status;
    int imin,imax,nump4;
    int v1,v2,v3,v4,e1,e2,e3,j1,j2,j3;

    n = sg.nv;
    if (vertex >= n) gt_abort(">E vertex given to -t is too large");

    status = isdecomposable(sg,FALSE,FALSE);
    if (status == NO) return -1;
    else if (status == TIMEOUT) return -2;

    if (vertex >= 0)
    {
        imin = imax = vertex; 
        DYNALLOC1(p4,p4list,p4list_sz,36,"malloc");
    }
    else
    {
	imin = 0; imax = n-1;
        DYNALLOC1(p4,p4list,p4list_sz,18*n,"malloc");
    }

    nump4 = 0;
    for (v1 = imin; v1 <= imax; ++v1)
    for (j1 = 0; j1 < 4; ++j1)
    {
	v2 = sg.e[4*v1+j1];
	e1 = eno[4*v1+j1];
	for (j2 = 0; j2 < 4; ++j2)
	{
	    v3 = sg.e[4*v2+j2];
	    if (v3 == v1) continue;
	    if (vertical && (v2|1) != (v3|1)) continue;
	    e2 = eno[4*v2+j2];
	    for (j3 = 0; j3 < 4; ++j3)
	    {
		v4 = sg.e[4*v3+j3];
		if (v4 == v1 || v4 == v2) continue;
		e3 = eno[4*v3+j3];
		if (vertex >= 0 || v1 < v4)
		{
		    p4list[nump4].v1 = v1;
		    p4list[nump4].v2 = v2;
		    p4list[nump4].v3 = v3;
		    p4list[nump4].v4 = v4;
		    p4list[nump4].e1 = e1;
		    p4list[nump4].e2 = e2;
		    p4list[nump4].e3 = e3;
		    p4list[nump4].ok = FALSE;
		    ++nump4;
		}
	    }
	}
    }

    for (i = 0; i < nump4; ++i)
        if (colour[p4list[i].e1] == colour[p4list[i].e2] &&
	    colour[p4list[i].e1] == colour[p4list[i].e3])
	    p4list[i].ok = TRUE;

    for (i = 0; i < nump4; ++i)
        if (!p4list[i].ok)
	{
	    initialise_colouring(n);
	    if (makeblue(p4list[i].e1,FALSE) &&
	        makeblue(p4list[i].e2,n==2) &&
		makeblue(p4list[i].e3,n==3))
	    {
		status = dispatchsearch(n,sg.e,3,0);
		if (status == TIMEOUT) return -2;
		if (status == YES)
		{
    		    for (k = 0; k < nump4; ++k)
    		    {
			if (p4list[k].ok) continue;
			if (colour[p4list[k].e1] == colour[p4list[k].e2]
                         && colour[p4list[k].e1] == colour[p4list[k].e3])
                           p4list[k].ok = TRUE;
                    }
		}
	    }
	}

    ans = 0;
    for (i = 0; i < nump4; ++i) if (!p4list[i].ok) ++ans;

    return ans;
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    sparsegraph sg,sh,*this;
    int n,codetype;
    int argnum,i,j,outcode;
    char *arg,sw;
    boolean badargs;
    boolean sswitch,gswitch,qswitch,vswitch,xswitch,Xswitch;
    boolean pswitch,Lswitch,tswitch,yswitch,Yswitch;
    long Lvalue;
    double t;
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    nauty_counter nin,nout;
    int status,xstatus,ystatus,tvalue,imin,imax;
    DYNALLSTAT(int,ham1,ham1_sz);
    DYNALLSTAT(int,ham2,ham2_sz);

    HELP; PUTVERSION;

    INITSEED;
    ran_init(seed);

    sswitch = gswitch = yswitch = qswitch = vswitch = FALSE;
    tswitch = xswitch = Xswitch = Lswitch = pswitch = FALSE;
    Yswitch = FALSE;
    infilename = outfilename = NULL;

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
                     SWBOOLEAN('s',sswitch)
                else SWBOOLEAN('g',gswitch)
                else SWBOOLEAN('q',qswitch)
                else SWBOOLEAN('x',xswitch)
                else SWBOOLEAN('y',yswitch)
                else SWBOOLEAN('Y',Yswitch)
                else SWBOOLEAN('X',Xswitch)
                else SWBOOLEAN('v',vswitch)
		else SWBOOLEAN('p',pswitch)
		else SWLONG('L',Lswitch,Lvalue,"-L")
		else SWINT('t',tswitch,tvalue,"-t")
                else badargs = TRUE;
            }
        }
        else
        {
            ++argnum;
            if      (argnum == 1) infilename = arg;
            else if (argnum == 2) outfilename = arg;
            else                  badargs = TRUE;
        }
    }

    if (sswitch && gswitch) 
        gt_abort(">E twohamg: -s and -g are incompatible");
    if ((xswitch != 0) + (Xswitch != 0)
               + (yswitch != 0) + (Yswitch != 0) > 1)
        gt_abort(">E twohamg: -x, -X, -y, -Y are incompatible");
    if (Yswitch && !pswitch)
        gt_abort(">E twohang: -Y is not allowed without -p");

    if (!Lswitch) totallimit = 0;
    else          totallimit = Lvalue * (nauty_counter)1000;

    if (!tswitch) tvalue = -1;

    if (badargs || argnum > 2)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!qswitch)
    {
        fprintf(stderr,">A twohamg");
        if (sswitch || gswitch || vswitch || xswitch || Xswitch 
                || tswitch || pswitch || Lswitch || yswitch || Yswitch)
            fprintf(stderr," -");
        if (sswitch) fprintf(stderr,"s");
        if (gswitch) fprintf(stderr,"g");
        if (xswitch) fprintf(stderr,"x");
        if (Xswitch) fprintf(stderr,"X");
        if (yswitch) fprintf(stderr,"y");
        if (Yswitch) fprintf(stderr,"Y");
        if (pswitch) fprintf(stderr,"p");
        if (vswitch) fprintf(stderr,"v");
	if (Lswitch) fprintf(stderr,"L%ld",Lvalue);
	if (tswitch) fprintf(stderr,"t%d",tvalue);
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,FALSE,1);
    if (!infile) exit(1);
    if (!infilename) infilename = "stdin";

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

    NODIGRAPHSYET(codetype);

    if (sswitch || (!gswitch && (codetype&SPARSE6)))
        outcode = SPARSE6;
    else 
        outcode = GRAPH6;

    if (codetype&HAS_HEADER)
    {
        if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
        else                    writeline(outfile,GRAPH6_HEADER);
    }

    nin = nout = 0;
    t = CPUTIME;
    SG_INIT(sg);
    SG_INIT(sh);

    while (read_sg(infile,&sg) != NULL)
    {
        ++nin;
        if (pswitch)
	{
            n = sg.nv;
            for (i = 0; i < n; ++i)
                if (sg.d[i] != 3) break;
            if (i != n)
            {
                fprintf(stderr,
                        ">W input " COUNTER_FMT " is not cubic\n",nin);
                continue;
            }
	    makeprism_sg(&sg,&sh);
	    n = sh.nv;
	    this = &sh;
	}
	else
	{
            n = sg.nv;
            for (i = 0; i < n; ++i)
                if (sg.d[i] != 4) break;
            if (i != n)
            {
                fprintf(stderr,
                        ">W input " COUNTER_FMT " is not quartic\n",nin);
                continue;
            }
            for (i = 0; i < n; ++i)
                if (sg.v[i] != 4*i) break;
            if (i != n)
                gt_abort("readg6_sg() result is not in standard form");
	    this = &sg;
	}

	if (xswitch || Xswitch)
	{
	    xstatus = iscrossdecomposable(*this,tvalue);
	    if ((xswitch && xstatus < 2) || (Xswitch && xstatus < 1))
	    {
		if (outcode == SPARSE6) writes6_sg(outfile,&sg);   
                else if (outcode == GRAPH6) writeg6_sg(outfile,&sg);
                ++nout;

		if (vswitch)
		{
		    if (xstatus < 0)
		    {
                        fprintf(stderr,">H Graph " COUNTER_FMT
                                     " is indecomposable\n",nin);
		    }
		    else
		    {
		        fprintf(stderr,">X" COUNTER_FMT ": ",nin);
			if (tswitch) { imin = imax = tvalue; }
			else         { imin = 0; imax = n-1; }
		        for (i = imin; i <= imax; ++i)
		        for (j = 1; j < 4; ++j)
			    if (!cross[4*i+j])
			    {
			        fprintf(stderr," %d-%d-%d",
			            (v1[eno[4*i]] == i
                                       ? v2[eno[4*i]] : v1[eno[4*i]]),
			            i,
			            (v1[eno[4*i+j]] == i
                                       ? v2[eno[4*i+j]] : v1[eno[4*i+j]]));
		            }
		        fprintf(stderr,"\n");
		    }
		}
	    }
	    else if (xstatus == 3)
	    {
		fprintf(stderr,">H Graph " COUNTER_FMT
                                     " timed out\n",nin);
		if (outcode == SPARSE6) writes6_sg(outfile,&sg);
                else if (outcode == GRAPH6) writeg6_sg(outfile,&sg);
                ++nout;
	    }
	}
	else if (yswitch || Yswitch)
        {
	    ystatus = p4decomposition(*this,tvalue,Yswitch);
	    if (ystatus != 0)
	    {
		if (outcode == SPARSE6) writes6_sg(outfile,&sg);   
                else if (outcode == GRAPH6) writeg6_sg(outfile,&sg);
                ++nout;
	    }

	    if (ystatus == -2)
	    {
		fprintf(stderr,">H Graph " COUNTER_FMT
                                     " timed out\n",nin);
	    }
	    else if (vswitch)
	    {
		if (ystatus == -1)
		{
                    fprintf(stderr,">H Graph " COUNTER_FMT
                                     " is indecomposable\n",nin);
		}
		else if (ystatus > 0)
		{
		    fprintf(stderr,">X" COUNTER_FMT ": ",nin);
		    for (i = j = 0; j < ystatus; ++i)
		        if (!p4list[i].ok)
			{
			    ++j;
			    fprintf(stderr," %d-%d-%d-%d",
				p4list[i].v1,p4list[i].v2,
				p4list[i].v3,p4list[i].v4);
		        }
		    fprintf(stderr,"\n");
		}
	    }
        }
	else
	{
            if (vswitch)
            {
	        DYNALLOC1(int,ham1,ham1_sz,n,"malloc");
	        DYNALLOC1(int,ham2,ham2_sz,n,"malloc");
                status = isdecomposable(*this,ham1,ham2);
            }
	    else
	        status = isdecomposable(*this,NULL,NULL);
    
            if (status != YES)
            {
                if (outcode == SPARSE6) writes6_sg(outfile,&sg);
                else if (outcode == GRAPH6) writeg6_sg(outfile,&sg);

		if (status == TIMEOUT)
		    fprintf(stderr,">H Graph " COUNTER_FMT
                                     " timed out\n",nin);
                else if (vswitch)
                    fprintf(stderr,">H Graph " COUNTER_FMT
                                     " is indecomposable\n",nin);
	        ++nout;
            }
            else if (vswitch)
            {
#if DEBUG
	        dumpdata(0,0,0,n);
#endif
                fprintf(stderr,">H" COUNTER_FMT ": ",nin);
	        for (i = 0; i < n; ++i) fprintf(stderr," %d",ham1[i]);
	        fprintf(stderr,"\n ");
	        for (i = 0; i < n; ++i) fprintf(stderr," %d",ham2[i]);
	        fprintf(stderr,"\n");
            }
	}
    }
    t = CPUTIME - t;

    if (!qswitch)
    {
	fprintf(stderr,">T to=" COUNTER_FMT " phases=",ntimeouts);
	for (i = 0; i <= NUMLIMITS; ++i)
	    fprintf(stderr," " COUNTER_FMT,nphases[i]);
	fprintf(stderr,"\n");

        fprintf(stderr,
                ">Z " COUNTER_FMT " graphs read from %s, "
                COUNTER_FMT " written to %s; %3.2f sec.\n",
                nin,infilename,nout,outfilename,t);
    }

    exit(0);
}
