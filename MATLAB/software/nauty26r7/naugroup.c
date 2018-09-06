/* naugroup.c

Procedures for handling groups found by nauty.
*/

#include "naugroup.h"

static permrec *freelist = NULL;
static int freelist_n = 0;

static grouprec *group = NULL;
static int group_depth = 0;
DYNALLSTAT(cosetrec,coset,coset_sz);
static permrec *gens;
DYNALLSTAT(set,workset,workset_sz);
DYNALLSTAT(int,allp,allp_sz);
DYNALLSTAT(int,id,id_sz);

/**************************************************************************/

permrec
*newpermrec(int n)
/* Get a permrec of order n.  This procedure and the next one are
designed to be efficient if lots of group ops are done with the
same value of n. */
{
    permrec *p;

    if (freelist_n != n)
    {
	while (freelist != NULL)
	{
	    p = freelist;
	    freelist = freelist->ptr;
	    free(p);
	}
	freelist_n = n;
    }

    if (freelist != NULL)
    {
	p = freelist;
	freelist = freelist->ptr;
	return p;
    }

    p = (permrec*) malloc(sizeof(permrec)+(freelist_n-2)*sizeof(int)); 

    if (p == NULL)
    {
	fprintf(ERRFILE,">E malloc failed in newpermrec()\n");
	exit(1);
    }

    return p;
}

/**************************************************************************/

void
freepermrec(permrec *p, int n)
/* Free a permrec of given size. */
{
    permrec *q;

    if (p == NULL) return;

    if (freelist_n != n)
    {
	while (freelist)
	{
	    q = freelist;
	    freelist = freelist->ptr;
	    free(q);
	}
	freelist_n = n;
    }

    p->ptr = freelist;
    freelist = p;
}

/**************************************************************************/

grouprec *
groupptr(boolean cutloose)
/* Give the address of the group structure, cutting it loose
   if requested. */
{
    grouprec *p;

    p = group;

    if (cutloose)
    {
	group = NULL;
	group_depth = 0;
	coset = NULL;
	coset_sz = 0;
    }

    return p;
}

/**************************************************************************/

void
freegroup(grouprec *grp)
/* Free (or pretend to free) group structure. */
{
    int i,j;
    cosetrec *p;
    permrec *q,*qq;

    for (i = 0; i < grp->depth; ++i)
    {
	p = grp->levelinfo[i].replist;
	if (p != NULL)
	    for (j = grp->levelinfo[i].orbitsize; --j >= 0; )
	    {
		freepermrec(p[j].rep,grp->n);
		p[j].rep = NULL;
	    }
    }

    if (grp->depth > 0)
    {
        p = grp->levelinfo[0].replist;
        if (p != NULL && p != coset)
	{
	    free(p);
	    grp->levelinfo[0].replist = NULL;
	}

        q = grp->levelinfo[0].gens;
        while (q != NULL)
        {
	    qq = q;
	    q = q->ptr;
	    freepermrec(qq,grp->n);
        }
	grp->levelinfo[0].gens = NULL;
    }
}

/**************************************************************************/

void
groupautomproc(int count, int *perm, int *orbits,
                                        int numorbits, int stabvertex, int n)
{
    permrec *p;
    int i;

    p = newpermrec(n);
    for (i = 0; i < n; ++i) p->p[i] = perm[i];
    p->ptr = gens;
    gens = p;
}

/**************************************************************************/

void
grouplevelproc(int *lab, int *ptn, int level, int *orbits, statsblk *stats,
               int tv, int index, int tcellsize, int numcells, int cc, int n)
{
    int depth;
    size_t sz;

    if (numcells == n)   /* first call */
    {
	depth = level - 1;

	if (group) freegroup(group);

	if (depth > group_depth || !group)
	{
	    if (depth <= 1) sz = sizeof(grouprec);
	    else            sz = sizeof(grouprec) + (depth-1)*sizeof(levelrec);
	    if (group) group = (grouprec*)realloc((void*)group,sz);
	    else       group = (grouprec*)malloc(sz);
	    if (group == NULL)
	    {
		fprintf(ERRFILE,">E malloc failed in grouplevelproc\n");
		exit(1);
	    }
	    group_depth = depth;
	}

	group->n = n;
	group->depth = depth;
	gens = NULL;
	return;
    }

    group->levelinfo[level-1].fixedpt = tv;
    group->levelinfo[level-1].orbitsize = index;
    group->levelinfo[level-1].gens = gens;
    group->levelinfo[level-1].replist = NULL;

    if (level == 1) group->numorbits = stats->numorbits;
}

/**************************************************************************/

void
makecosetreps(grouprec *grp)
/* Make all coset representatives for this group */
{
    int i,j,k,n,depth;
    int l,index;
    int *p,*q;
    permrec *gen,*g;
    cosetrec *cr;
    int head,tail;
    DYNALLSTAT(int,queue,queue_sz);
    DYNALLSTAT(int,lab,lab_sz);

    n = grp->n;
    depth = grp->depth;

    DYNALLOC1(int,queue,queue_sz,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");

    j = 0;
    for (i = 0; i < depth; ++i)
	j += grp->levelinfo[i].orbitsize;

    if (j > 0) DYNALLOC1(cosetrec,coset,coset_sz,j,"malloc");

    cr = coset;
    for (i = 0; i < depth; ++i)
    {
	grp->levelinfo[i].replist = cr;
	cr += grp->levelinfo[i].orbitsize;
    }

    for (i = 0; i < depth; ++i)
    {
	cr = grp->levelinfo[i].replist;
	gen = grp->levelinfo[i].gens;
	for (j = 0; j < n; ++j) lab[j] = -1;
	queue[0] = grp->levelinfo[i].fixedpt;
	lab[queue[0]] = 0;
	cr[0].image = queue[0];
	cr[0].rep = NULL;
	head = 0;
	tail = 1;
	index = 0;
	while (head < tail)
	{
	    j = queue[head++];
	    p = (cr[lab[j]].rep ? cr[lab[j]].rep->p : NULL);
	    for (g = gen; g != NULL; g = g->ptr)
	    {
		k = g->p[j];
		if (lab[k] < 0)
		{
		    ++index;
		    lab[k] = index;
		    queue[tail++] = k;
		    cr[index].image = k;
		    cr[index].rep = newpermrec(n);
		    q = cr[index].rep->p;
		    if (p == NULL)
			for (l = 0; l < n; ++l) q[l] = g->p[l];
		    else
			for (l = 0; l < n; ++l) q[l] = g->p[p[l]];
		}
	    }
	}
    }
}

/**************************************************************************/

int
permcycles(int *p, int n, int *len, boolean sort)
/* Puts in len[0..] the cycle lengths of p.  If sort, sort them. 
   Return the number of cycles. */
{
    int m,i,j,k,h,nc,leni;

    m = (n + WORDSIZE - 1) / WORDSIZE;
    DYNALLOC1(set,workset,workset_sz,m,"malloc");

    EMPTYSET(workset,m);

    nc = 0;
    for (i = 0; i < n; ++i)
        if (!ISELEMENT(workset,i))
	{
	    k = 1;
	    for (j = p[i]; j != i; j = p[j]) 
	    {
		ADDELEMENT(workset,j);
		++k;
	    }
	    len[nc++] = k;
	}

    if (sort && nc > 1)
    {
	j = nc / 3;
        h = 1;
        do
            h = 3 * h + 1;
        while (h < j);

	do
        {
            for (i = h; i < nc; ++i)
            {
                leni = len[i];
                for (j = i; len[j-h] > leni; )
                {
                    len[j] = len[j-h];
                    if ((j -= h) < h) break;
                }
                len[j] = leni;
            }
            h /= 3;
        }
        while (h > 0);
    }

    return nc;
}

/**************************************************************************/

static void
groupelts(levelrec *lr, int n, int level, void (*action)(int*,int),
          int *before, int *after, int *id)
/* Recursive routine used by allgroup. */
{
    int i,j,orbsize;
    int *p,*cr;
    cosetrec *coset;
    
    coset = lr[level].replist;
    orbsize = lr[level].orbitsize;

    for (j = 0; j < orbsize; ++j)
    {
	cr = (coset[j].rep == NULL ? NULL : coset[j].rep->p);
	if (before == NULL)
	    p = cr;
	else if (cr == NULL)
	    p = before;
	else
	{
	    p = after;
	    for (i = 0; i < n; ++i) p[i] = cr[before[i]];
	}

	if (level == 0) 
	    (*action)((p == NULL ? id : p),n);
	else
	    groupelts(lr,n,level-1,action,p,after+n,id);
    }
}

/**************************************************************************/

void
allgroup(grouprec *grp, void (*action)(int*,int))
/* Call action(p,n) for every element of the group, including the identity. 
   The identity is always the first call. */
{
    int i,depth,n;

    depth = grp->depth;
    n = grp->n;

    DYNALLOC1(int,id,id_sz,n,"malloc");
    for (i = 0; i < n; ++i) id[i] = i;

    if (depth == 0)
    {
	(*action)(id,n);
	return;
    }

    DYNALLOC1(int,allp,allp_sz,n*depth,"malloc");

    groupelts(grp->levelinfo,n,depth-1,action,NULL,allp,id);
}

/**************************************************************************/

static void
groupelts2(levelrec *lr, int n, int level,
    void (*action)(int*,int,int*), int *before,
    int *after, int *id, int *abort)
/* Recursive routine used by allgroup2. */
{
    int i,j,orbsize;
    int *p,*cr;
    cosetrec *coset;
    
    coset = lr[level].replist;
    orbsize = lr[level].orbitsize;

    for (j = 0; j < orbsize; ++j)
    {
	cr = (coset[j].rep == NULL ? NULL : coset[j].rep->p);
	if (before == NULL)
	    p = cr;
	else if (cr == NULL)
	    p = before;
	else
	{
	    p = after;
	    for (i = 0; i < n; ++i) p[i] = cr[before[i]];
	}

	if (level == 0) 
	    (*action)((p == NULL ? id : p),n,abort);
	else
	    groupelts2(lr,n,level-1,action,p,after+n,id,abort);
	if (*abort) return;
    }
}

/**************************************************************************/

int
allgroup2(grouprec *grp, void (*action)(int*,int,int*))
/* Call action(p,n,&abort) for every element of the group, including
   the identity.  The identity is always the first call.
   If action() stores a non-zero value in abort, group generation is
   aborted and the abort value is returned by this procedure.  If no
   non-zero value is ever returned in abort by action(), this
   procedure returns 0. */
{
    int i,depth,n,abort;

    depth = grp->depth;
    n = grp->n;

    DYNALLOC1(int,id,id_sz,n,"malloc");
    for (i = 0; i < n; ++i) id[i] = i;

    abort = 0;
    if (depth == 0)
    {
	(*action)(id,n,&abort);
	return abort;
    }

    DYNALLOC1(int,allp,allp_sz,n*depth,"malloc");

    groupelts2(grp->levelinfo,n,depth-1,action,NULL,allp,id,&abort);

    return abort;
}

/**************************************************************************/

static void
groupelts3(levelrec *lr, int n, int level,
    void (*action)(int*,int,int*,void*), int *before,
    int *after, int *id, int *abort, void *userptr)
/* Recursive routine used by allgroup3. */
{
    int i,j,orbsize;
    int *p,*cr;
    cosetrec *coset;
    
    coset = lr[level].replist;
    orbsize = lr[level].orbitsize;

    for (j = 0; j < orbsize; ++j)
    {
	cr = (coset[j].rep == NULL ? NULL : coset[j].rep->p);
	if (before == NULL)
	    p = cr;
	else if (cr == NULL)
	    p = before;
	else
	{
	    p = after;
	    for (i = 0; i < n; ++i) p[i] = cr[before[i]];
	}

	if (level == 0) 
	    (*action)((p == NULL ? id : p),n,abort,userptr);
	else
	    groupelts3(lr,n,level-1,action,p,after+n,id,abort,userptr);
	if (*abort) return;
    }
}

/**************************************************************************/

int
allgroup3(grouprec *grp, void (*action)(int*,int,int*,void*), void *userptr)
/* Call action(p,n,&abort,userptr) for every element of the group,
   including the identity.  The identity is always the first call.
   If action() stores a non-zero value in abort, group generation is
   aborted and the abort value is returned by this procedure.  If no
   non-zero value is ever returned in abort by action(), this
   procedure returns 0. The pointer userptr is not interpretted and
   is passed to action() to use as it likes. */
{
    int i,depth,n,abort;

    depth = grp->depth;
    n = grp->n;

    DYNALLOC1(int,id,id_sz,n,"malloc");
    for (i = 0; i < n; ++i) id[i] = i;

    abort = 0;
    if (depth == 0)
    {
	(*action)(id,n,&abort,userptr);
	return abort;
    }

    DYNALLOC1(int,allp,allp_sz,n*depth,"malloc");

    groupelts3(grp->levelinfo,n,depth-1,action,NULL,allp,id,&abort,userptr);

    return abort;
}

