/* schreier.c - procedures for manipulating a permutation group using
 * the random schreier algorithm.  There is a separate file schreier.txt
 * which describes the usage.
 *
 * Written for nauty and traces, Brendan McKay 2010-2013.
 */

#include "schreier.h" 

TLS_ATTR long long multcount = 0;
TLS_ATTR long long filtercount = 0;

static permnode id_permnode; 
  /* represents identity, no actual content, doesn't need TLS_ATTR */
#define ID_PERMNODE (&id_permnode)

#if !MAXN
DYNALLSTAT(int,workperm,workperm_sz);
DYNALLSTAT(int,workperm2,workperm2_sz);
DYNALLSTAT(int,workpermA,workpermA_sz);
DYNALLSTAT(int,workpermB,workpermB_sz);
DYNALLSTAT(set,workset,workset_sz);
DYNALLSTAT(set,workset2,workset2_sz);
#else
static TLS_ATTR int workperm[MAXN];
static TLS_ATTR int workperm2[MAXN];
static TLS_ATTR int workpermA[MAXN];
static TLS_ATTR int workpermB[MAXN];
static TLS_ATTR set workset[MAXM];
static TLS_ATTR set workset2[MAXM];
#endif

static TLS_ATTR schreier *schreier_freelist = NULL;
	/* Freelist of scheier structures connected by next field.
         * vec, pwr and orbits fields are assumed allocated. */
static TLS_ATTR permnode *permnode_freelist = NULL;
	/* Freelist of permnode structures connected by next field.
         * p[] is assumed extended. */

static TLS_ATTR int schreierfails = SCHREIERFAILS;

#define TMP

static boolean filterschreier(schreier*,int*,permnode**,boolean,int,int);
#define PNCODE(x) ((int)(((size_t)(x)>>3)&0xFFFUL))

/* #define TESTP(id,p,n) testispermutation(id,p,n) */
#define TESTP(id,p,n)

/************************************************************************/

static void
testispermutation(int id, int *p, int n)
/* For debugging purposes, crash with a message if p[0..n-1] is
   not a permutation. */
{
    int i,m;
    DYNALLSTAT(set,seen,seen_sz);

    for (i = 0; i < n; ++i)
        if (p[i] < 0 || p[i] > n) break;

    if (i < n)
    {
	fprintf(stderr,">E Bad permutation (id=%d): n=%d p[%d]=%d\n",
		id,n,i,p[i]);
	exit(1);
    }

    m = SETWORDSNEEDED(n);
    DYNALLOC1(set,seen,seen_sz,m,"malloc seen");
    EMPTYSET(seen,m);

    for (i = 0; i < n; ++i)
    {
	if (ISELEMENT(seen,p[i]))
	{
	    fprintf(stderr,
		">E Bad permutation (id=%d): n=%d p[%d]=%d is a repeat\n",
		id,n,i,p[i]);
	    exit(1);
	}
        ADDELEMENT(seen,p[i]);
    }
}
    
/************************************************************************/

int
schreier_fails(int nfails)
/* Set the number of consecutive failures for filtering;
 * A value of <= 0 defaults to SCHREIERFAILS.
 * The function value is the previous setting. */
{
    int prev;

    prev = schreierfails;

    if (nfails <= 0) schreierfails = SCHREIERFAILS;
    else             schreierfails = nfails;

    return prev;
}

/************************************************************************/

static void
clearfreelists(void)
/* Clear the schreier and permnode freelists */
{
    schreier *sh,*nextsh;
    permnode *p,*nextp;

    nextsh = schreier_freelist;
    while (nextsh)
    {
	sh = nextsh;
	nextsh = sh->next;
	free(sh->vec);
	free(sh->pwr);
	free(sh->orbits);
	free(sh);
    }
    schreier_freelist = NULL;

    nextp = permnode_freelist;
    while (nextp)
    {
	p = nextp;
	nextp = p->next;
	free(p);
    }
    permnode_freelist = NULL;
}

/************************************************************************/

static permnode
*newpermnode(int n)
/* Allocate a new permode structure, with initialized next/prev fields */
{
    permnode *p;

    while (permnode_freelist)
    {
	p = permnode_freelist;
	permnode_freelist = p->next;
	if (p->nalloc >= n && p->nalloc <= n+100)
	{
	    p->next = p->prev = NULL;
	    p->mark = 0;
            return p;
	}
	else
	    free(p);
    }

    p = (permnode*) malloc(sizeof(permnode)+(n-2)*sizeof(int)); 

    if (p == NULL)
    {
        fprintf(ERRFILE,">E malloc failed in newpermnode()\n");
        exit(1);
    }

    p->next = p->prev = NULL;
    p->nalloc = n;

    return p;
}

/************************************************************************/

static schreier
*newschreier(int n)
/* Allocate a new schreier structure, with initialised next field */
{
    schreier *sh;

    while (schreier_freelist)
    {
	sh = schreier_freelist;
	schreier_freelist = sh->next;
	if (sh->nalloc >= n && sh->nalloc <= n+100)
	{
	    sh->next = NULL;
	    return sh;
	}
	else
	{
	    free(sh->vec);
	    free(sh->pwr);
	    free(sh->orbits);
	    free(sh);
	}
    }

    sh = (schreier*) malloc(sizeof(schreier));

    if (sh == NULL)
    {
        fprintf(ERRFILE,">E malloc failed in newschreier()\n");
        exit(1);
    }

    sh->vec = (permnode**) malloc(sizeof(permnode*)*n);
    sh->pwr = (int*) malloc(sizeof(int)*n);
    sh->orbits = (int*) malloc(sizeof(int)*n);

    if (sh->vec == NULL || sh->pwr == NULL || sh->orbits == NULL)
    {
        fprintf(ERRFILE,">E malloc failed in newschreier()\n");
        exit(1);
    }

    sh->next = NULL;
    sh->nalloc = n;

    return sh;
}

/************************************************************************/

void
freeschreier(schreier **gp, permnode **gens)
/* Free schreier structure and permutation ring.  Assume this is everything. */
/* Use NULL for arguments which don't need freeing. */
{
    schreier *sh,*nextsh;
    permnode *p,*nextp;

    if (gp && *gp)
    {
        nextsh = *gp;
        while (nextsh)
        {
	    sh = nextsh;
	    nextsh = sh->next;
	    sh->next = schreier_freelist;
	    schreier_freelist = sh;
        }
        *gp = NULL;
    }

    if (gens && *gens)
    {
        p = *gens;
        do
        {
            nextp = p->next;
	    p->next = permnode_freelist;
	    permnode_freelist = p;
	    p = nextp;
        } while (p != *gens);
        *gens = NULL;
    }
}

/************************************************************************/

permnode*
findpermutation(permnode *pn, int *p, int n)
/* Return a pointer to permutation p in the circular list,
 * or NULL if it isn't present. */
{
    permnode *rn;
    int i;

    if (!pn) return NULL;

    rn = pn;
    do
    {
	for (i = 0; i < n; ++i) 
	    if (rn->p[i] != p[i]) break;
	if (i == n) return rn;
	rn = rn->next;
    } while (rn != pn);

    return NULL;
}

/************************************************************************/

void
addpermutation(permnode **ring, int *p, int n)
/* Add new permutation to circular list, marked.
 * and return pointer to it in *ring. */
{
    permnode *pn,*rn;

    pn = newpermnode(n);
    rn = *ring;

    memcpy(pn->p,p,n*sizeof(int));
	
    if (!rn)
	pn->next = pn->prev = pn;
    else
    {
	pn->next = rn->next;
	pn->prev = rn;
	rn->next = pn->next->prev = pn;
    }

    pn->refcount = 0;
    pn->mark = 1;
    *ring = pn;
}

/************************************************************************/

static void
addpermutationunmarked(permnode **ring, int *p, int n)
/* Add new permutation to circular list, not marked.
 * and return pointer to it in *ring. */
{
    TESTP(3,p,n);
    addpermutation(ring,p,n);
    (*ring)->mark = 0;
}

/************************************************************************/

boolean
addgenerator(schreier **gp, permnode **ring, int *p, int n)
/* Add new permutation to group, unless it is discovered to be
 * already in the group.  It is is possible to be in the group
 * and yet this fact is not discovered.
 * Return TRUE if the generator (or an equivalent) is added or the
 * group knowledge with the current partial base is improved. */
{
    TESTP(2,p,n);
    return filterschreier(*gp,p,ring,FALSE,-1,n);
}

/************************************************************************/

boolean
condaddgenerator(schreier **gp, permnode **ring, int *p, int n)
/* Add new permutation to group, unless it is discovered to be
 * already in the group.  It is is possible to be in the group
 * and yet this fact is not discovered, but this version will 
 * always notice if this permutation precisely is present.
 * Return TRUE if the generator (or an equivalent) is added or the
 * group knowledge with the current partial base is improved. */
{
    TESTP(4,p,n);
    if (findpermutation(*ring,p,n)) 
	return FALSE;
    else
        return filterschreier(*gp,p,ring,FALSE,-1,n);
}

/************************************************************************/

static void
delpermnode(permnode **ring)
/* Delete permnode at head of circular list, making the next node head. */
{
    permnode *newring;

    if (!*ring) return;

    if ((*ring)->next == *ring)
	newring = NULL;
    else
    {
	newring = (*ring)->next;
	newring->prev = (*ring)->prev;
        (*ring)->prev->next = newring;
    }

    (*ring)->next = permnode_freelist;
    permnode_freelist = *ring;

    *ring = newring;
}

/************************************************************************/

void
deleteunmarked(permnode **ring)
/* Delete all permutations in the ring that are not marked */
{
    permnode *pn,*firstmarked;

    pn = *ring;
    firstmarked = NULL;

    while (pn != NULL && pn != firstmarked)
    {
	if (pn->mark)
	{
	    if (!firstmarked) firstmarked = pn;
	    pn = pn->next;
	}
	else
	    delpermnode(&pn);
    }

    *ring = pn;
}

/************************************************************************/

static void
clearvector(permnode **vec, permnode **ring, int n)
/* clear vec[0..n-1], freeing permnodes that have no other references
 * and are not marked */
{
    int i;

    for (i = 0; i < n; ++i)
        if (vec[i])
        {
            if (vec[i] != ID_PERMNODE)
	    {
		--(vec[i]->refcount);
		if (vec[i]->refcount == 0 && !vec[i]->mark)
		{
		    *ring = vec[i];
		    delpermnode(ring);
		}
	    }
            vec[i] = NULL;
        }
}

/************************************************************************/

static void
initschreier(schreier *sh, int n)
/* Initialise schreier structure to trivial orbits and empty vector */
{
    int i;

    sh->fixed = -1;
    for (i = 0; i < n; ++i)
    {
	sh->vec[i] = NULL; 
	sh->orbits[i] = i;
    }
}

/************************************************************************/

void
newgroup(schreier **sh, permnode **ring, int n)
/* Make the trivial group, allow for ring to be set elsewhere */
{
    *sh = newschreier(n);
    initschreier(*sh,n);
    if (ring) *ring = NULL;
}

/************************************************************************/

static void
applyperm(int *wp, int *p, int k, int n)
/* Apply the permutation p, k times to each element of wp */
{
    int i,j,cyclen,kk,m;

    TESTP(1,p,n);

    if (k <= 5)
    {
        if (k == 0)
            return;
        else if (k == 1)
            for (i = 0; i < n; ++i) wp[i] = p[wp[i]];
        else if (k == 2)
            for (i = 0; i < n; ++i) wp[i] = p[p[wp[i]]];
        else if (k == 3)
            for (i = 0; i < n; ++i) wp[i] = p[p[p[wp[i]]]];
        else if (k == 4)
            for (i = 0; i < n; ++i) wp[i] = p[p[p[p[wp[i]]]]];
        else if (k == 5)
            for (i = 0; i < n; ++i) wp[i] = p[p[p[p[p[wp[i]]]]]];
    }
    else if (k <= 19)
    {
#if !MAXN
        DYNALLOC1(int,workpermA,workpermA_sz,n,"applyperm");
#endif
        for (i = 0; i < n; ++i) workpermA[i] = p[p[p[i]]];
        for (; k >= 6; k -= 6)
            for (i = 0; i < n; ++i) wp[i] = workpermA[workpermA[wp[i]]];
        if (k == 1)
            for (i = 0; i < n; ++i) wp[i] = p[wp[i]];
        else if (k == 2)
            for (i = 0; i < n; ++i) wp[i] = p[p[wp[i]]];
        else if (k == 3)
            for (i = 0; i < n; ++i) wp[i] = workpermA[wp[i]];
        else if (k == 4)
            for (i = 0; i < n; ++i) wp[i] = p[workpermA[wp[i]]];
        else if (k == 5)
            for (i = 0; i < n; ++i) wp[i] = p[p[workpermA[wp[i]]]];
    }
    else
    {
        m = SETWORDSNEEDED(n);
#if !MAXN
        DYNALLOC1(int,workpermA,workpermA_sz,n,"applyperm");
        DYNALLOC1(int,workpermB,workpermB_sz,n,"applyperm");
        DYNALLOC1(set,workset2,workset2_sz,m,"applyperm");
#endif

        EMPTYSET(workset2,m);

      /* We will construct p^k in workpermB one cycle at a time. */

        for (i = 0; i < n; ++i)
        {
            if (ISELEMENT(workset2,i)) continue;
            if (p[i] == i)
                workpermB[i] = i;
            else
            {
                cyclen = 1;
                workpermA[0] = i;
                for (j = p[i]; j != i; j = p[j])
                {
                    workpermA[cyclen++] = j;
                    ADDELEMENT(workset2,j);
                }
                kk = k % cyclen;
                for (j = 0; j < cyclen; ++j)
                {
                    workpermB[workpermA[j]] = workpermA[kk];
                    if (++kk == cyclen) kk = 0;
                }
            }
        }
        for (i = 0; i < n; ++i) wp[i] = workpermB[wp[i]];
    }
}

/************************************************************************/

static boolean
filterschreier(schreier *gp, int *p, permnode **ring,
               boolean ingroup, int maxlevel, int n)
/* Filter permutation p up to level maxlevel of gp.
 * Use ingroup=TRUE if p is known to be in the group, otherwise
 * at least one equivalent generator is added unless it is proved
 * (nondeterministically) that it is in the group already. 
 * maxlevel < 0 means no limit, maxlevel=0 means top level only, etc.
 * Return TRUE iff some change is made. */
{
    int i,j,j1,j2,lev;
    int ipwr;
    schreier *sh;
    int *orbits,*pwr;
    permnode **vec,*curr;
    boolean changed,lchanged,ident;
#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n,"filterschreier");
#endif

++filtercount;

    memcpy(workperm,p,n*sizeof(int));

    if (*ring && p == (*ring)->p)
    {
	ingroup = TRUE;
	curr = *ring;
    }
    else
	curr = NULL;

 /* curr is the location of workperm in ring, if anywhere */

    sh = gp;
    changed = FALSE;
    if (maxlevel < 0) maxlevel = n+1;
 
    for (lev = 0; lev <= maxlevel; ++lev)
    {
	for (i = 0; i < n; ++i) if (workperm[i] != i) break;
	ident = (i == n);
	if (ident) break;

	lchanged = FALSE;
	orbits = sh->orbits;
	vec = sh->vec;
        pwr = sh->pwr;
	for (i = 0; i < n; ++i)
	{
	    j1 = orbits[i];
            while (orbits[j1] != j1) j1 = orbits[j1];
	    j2 = orbits[workperm[i]];
            while (orbits[j2] != j2) j2 = orbits[j2];

	    if (j1 != j2)
	    {
		lchanged = TRUE;
		if (j1 < j2) orbits[j2] = j1;
		else         orbits[j1] = j2;
	    }
	}
	if (lchanged)
	    for (i = 0; i < n; ++i) orbits[i] = orbits[orbits[i]];

	if (lchanged) changed = TRUE;
	
	if (sh->fixed >= 0)
	{
	    for (i = 0; i < n; ++i)
	        if (vec[i] && !vec[workperm[i]])
	        {
		    changed = TRUE;
		    ipwr = 0;
		    for (j = workperm[i]; !vec[j] ; j = workperm[j]) ++ipwr;

		    for (j = workperm[i]; !vec[j] ; j = workperm[j])
		    {
			if (!curr)
			{
			    if (!ingroup) addpermutation(ring,workperm,n);
			    else  addpermutationunmarked(ring,workperm,n);
			    ingroup = TRUE;
			    curr = *ring;
			}
			vec[j] = curr;
			pwr[j] = ipwr--;
			++curr->refcount;
		    }
	        }

	    j = workperm[sh->fixed];

	    while (j != sh->fixed)
	    {
		applyperm(workperm,vec[j]->p,pwr[j],n);
		++multcount;
		curr = NULL;
		j = workperm[sh->fixed];
	    }
	    sh = sh->next;
	}
	else
	    break;
    }

    if (!ident && !ingroup)
    {
	changed = TRUE;
	addpermutation(ring,p,n);
    }

    return changed;
}

/************************************************************************/

boolean
expandschreier(schreier *gp, permnode **ring, int n)
/* filter random elements until schreierfails failures.
 * Return true if it ever expanded. */
{
    int i,j,nfails,wordlen,skips;
    boolean changed;
    permnode *pn;
#if !MAXN
    DYNALLOC1(int,workperm2,workperm2_sz,n,"expandschreier");
#endif

    pn = *ring;
    if (pn == NULL) return FALSE;

    nfails = 0;
    changed = FALSE;

    for (skips = KRAN(17); --skips >= 0; ) pn = pn->next;
 
    memcpy(workperm2,pn->p,n*sizeof(int));

    while (nfails < schreierfails)
    {
	wordlen = 1 + KRAN(3);
	for (j = 0; j < wordlen; ++j)
	{
	    for (skips = KRAN(17); --skips >= 0; ) pn = pn->next;
	    for (i = 0; i < n; ++i) workperm2[i] = pn->p[workperm2[i]];
	}
	if (filterschreier(gp,workperm2,ring,TRUE,-1,n))
	{
	    changed = TRUE;
	    nfails = 0;
	}
	else
	    ++nfails;
    }

    return changed;
}

/************************************************************************/

int*
getorbits(int *fix, int nfix, schreier *gp, permnode **ring, int n)
/* Get a pointer to the orbits for this partial base. The pointer
 * remains valid until pruneset(), getorbits(), getorbitsmin()
 * or grouporder() is called with an incompatible base (neither a
 * prefix nor an extension). The contents of the array pointed to
 * MUST NOT BE MODIFIED by the calling program.
 */
{
    int k;
    schreier *sh,*sha;

    sh = gp;
    for (k = 0; k < nfix; ++k)
    {
	if (sh->fixed != fix[k]) break;
	sh = sh->next;
    }

    if (k == nfix) return sh->orbits;
 
    sh->fixed = fix[k];
    clearvector(sh->vec,ring,n);
    sh->vec[fix[k]] = ID_PERMNODE;  

    for (sha = sh->next; sha ; sha = sha->next) clearvector(sha->vec,ring,n);
    
    for (++k; k <= nfix; ++k)
    {
	if (!sh->next) sh->next = newschreier(n);
	sh = sh->next;
	initschreier(sh,n);
	if (k < nfix)
	{
	    sh->fixed = fix[k];
	    sh->vec[fix[k]] = ID_PERMNODE;
	}
	else
	    sh->fixed = -1;
    }

    if (*ring) expandschreier(gp,ring,n);
    return sh->orbits;
}

/************************************************************************/

int
getorbitsmin(int *fix, int nfix, schreier *gp, permnode **ring,
         int **orbits, int *cell, int ncell, int n, boolean changed)
/* If the basis elements fix[0..nfix-1] are minimal in their orbits,
 * as far as we know, return value nfix and set *orbits to point
 * to orbits fixing fix[0..nfix-1]. If fix[i] is seen to be not
 * minimal for some i <= nfix-1, return i and set *orbits to point
 * to orbits fixing fix[0..i-1]. If the partial base is already
 * known, or fix[0..nfix-1] can already be seen to be non-minimal,
 * do this work without more filtering. This shortcut is turned
 * off if changed==TRUE. Otherwise, filter until schreierfails
 * failures.
 * The pointer returned remains valid until pruneset(), getorbits(),
 * getorbitsmin() or grouporder() is called with an incompatible base
 * (neither a prefix nor an extension). The contents of the array
 * pointed to MUST NOT BE MODIFIED by the calling program.
 * If cell != NULL, return early if possible when cell[0..ncell-1]
 * are all in the same orbit fixing fix[0..nfix-1]. Otherwise
 * cell,ncell play no part in the computation.
 */
{
    schreier *sh,*sha;
    int *fixorbs;
    int i,j,k,icell,nfails,wordlen,skips;
    permnode *pn;
#if !MAXN
    DYNALLOC1(int,workperm2,workperm2_sz,n,"expandschreier");
#endif

    sh = gp;
    k = 0;
    if (!changed)
        for (k = 0; k < nfix; ++k)
        {
	    if (sh->orbits[fix[k]] != fix[k])
	    {
	        *orbits = sh->orbits;
	        return k;
	    }
	    if (sh->fixed != fix[k]) break;
	    sh = sh->next;
        }

    if (k == nfix)
    {
	*orbits = sh->orbits;
	return nfix;
    }
 
    sh->fixed = fix[k];
    clearvector(sh->vec,ring,n);
    sh->vec[fix[k]] = ID_PERMNODE;  

    for (sha = sh->next; sha ; sha = sha->next) clearvector(sha->vec,ring,n);
    
    for (++k; k <= nfix; ++k)
    {
	if (!sh->next) sh->next = newschreier(n);
	sh = sh->next;
	initschreier(sh,n);
	if (k < nfix)
	{
	    sh->fixed = fix[k];
	    sh->vec[fix[k]] = ID_PERMNODE;
	}
	else
	    sh->fixed = -1;
    }
    *orbits = fixorbs = sh->orbits;

    if (cell)
    {
	for (icell = 1; icell < ncell; ++icell)
	    if (fixorbs[cell[icell]] != fixorbs[cell[0]]) break;

	if (icell >= ncell) return nfix;
    }

    if (*ring)
    {
        pn = *ring;

        nfails = 0;

        for (skips = KRAN(17); --skips >= 0; ) pn = pn->next;

        memcpy(workperm2,pn->p,n*sizeof(int));

        while (nfails < schreierfails)
        {
            wordlen = 1 + KRAN(3);
            for (j = 0; j < wordlen; ++j)
            {
                for (skips = KRAN(17); --skips >= 0; ) pn = pn->next;
                for (i = 0; i < n; ++i) workperm2[i] = pn->p[workperm2[i]];
            }
            if (filterschreier(gp,workperm2,ring,TRUE,-1,n))
            {
                nfails = 0;
                sh = gp;
                for (k = 0; k < nfix; ++k)
                {
	            if (sh->orbits[fix[k]] != fix[k])
	            {
	                *orbits = sh->orbits;
	                return k;
	            }
	            sh = sh->next;
		}
		if (cell)
		{
		    for ( ; icell < ncell; ++icell)
	    	 	if (fixorbs[cell[icell]] != fixorbs[cell[0]]) break;

		    if (icell >= ncell) return nfix;
		}
	    }
            else
                ++nfails;
        }
    }

    return nfix;
}

/************************************************************************/

void
pruneset(set *fixset, schreier *gp, permnode **ring, set *x, int m, int n)
/* Remove from x any point not minimal for the orbits for this base.
 * If the base is already known, just provide the orbits without
 * more filtering. Otherwise, filter until schreierfails failures.
 */ 
{
    int i,k;
    schreier *sh,*sha;
    int *orbits;

#if !MAXN
    DYNALLOC1(set,workset,workset_sz,m,"pruneset");
#endif
    for (i = 0; i < m; ++i) workset[i] = fixset[i];
	
    sh = gp;
    while (sh->fixed >= 0 && ISELEMENT(workset,sh->fixed))
    {
	DELELEMENT(workset,sh->fixed);
	sh = sh->next;
    }

    k = nextelement(workset,m,-1);
    if (k < 0)
	orbits = sh->orbits;
    else
    {
        sh->fixed = k;
        clearvector(sh->vec,ring,n);
        sh->vec[k] = ID_PERMNODE;  

        for (sha = sh->next; sha ; sha = sha->next)
	    clearvector(sha->vec,ring,n);
    
	while ((k = nextelement(workset,m,k)) >= 0)
        {
	    if (!sh->next) sh->next = newschreier(n);
	    sh = sh->next;
	    initschreier(sh,n);
	    sh->fixed = k;
	    sh->vec[k] = ID_PERMNODE;
        }
	if (!sh->next) sh->next = newschreier(n);
	sh = sh->next;
        initschreier(sh,n);
	sh->fixed = -1;

        if (*ring) expandschreier(gp,ring,n);
        orbits = sh->orbits;
    }

    for (k = -1; (k = nextelement(x,m,k)) >= 0; )
	if (orbits[k] != k) DELELEMENT(x,k);
}

/************************************************************************/

int
schreier_gens(permnode *ring)
/* Returns the number of generators in the ring */
{
    int j;
    permnode *pn;

    if (!ring) j = 0;
    else for (j = 1, pn = ring->next; pn != ring; pn = pn->next) ++j;
 
    return j;
}

/************************************************************************/

void
dumpschreier(FILE *f, schreier *gp, permnode *ring, int n)
/* Dump the whole schreier structure to file f. */
{
    schreier *sh;
    permnode *pn;
    int i,j,jj,k;


    fprintf(f,"Schreier structure n=%d; ",n);

    jj = -1;
    for (j = 0, sh = gp; sh; sh = sh->next)
    {
	++j;
	if (sh->fixed < 0 && jj < 0) jj = j;
    }
    fprintf(f," levels=%d (%d used); ",j,jj);

    if (!ring) j = 0;
    else for (j = 1, pn = ring->next; pn != ring; pn = pn->next) ++j;
    fprintf(f,"gens=%d; ",j);

    for (j = 0, sh = schreier_freelist; sh; sh = sh->next) ++j;
    for (k = 0, pn = permnode_freelist; pn; pn = pn->next) ++k;
    fprintf(f,"freelists: %d,%d\n",j,k);

    if (ring)
    {
	fprintf(f,"Generators:\n");
        pn = ring;
	do
	{
	    fprintf(f,"  %03x ref=%lu mk=%d alloc=%d p=",PNCODE(pn),
			pn->refcount,pn->mark,pn->nalloc);
	    for (i = 0; i < n; ++i) fprintf(f," %d",pn->p[i]);
	    fprintf(f,"\n");
	    pn = pn->next;
	} while (pn != ring);
    }

    if (gp)
    {
	fprintf(f,"Levels:\n");
	for (sh = gp; sh; sh = sh->next)
	{
	    fprintf(f,"fixed=%2d alloc=%d vec=",sh->fixed,sh->nalloc);
	    for (i = 0; i < n; ++i)
	    {
		if (sh->vec[i] == ID_PERMNODE) fprintf(f," %d=e",i);
		else if (sh->vec[i])
		{
		    k = sh->pwr[i];
		    j = (sh->vec[i])->p[i];
		    fprintf(f," %03x",PNCODE(sh->vec[i]));
		    if (k == 1)
			fprintf(f,"(%d,%d)",i,j);
		    else
		    {
			fprintf(f,"^%d",k);
			while (--k >= 1) j = (sh->vec[i])->p[j];
			fprintf(f,"(%d,%d)",i,j);
		    }
		}
	    }
	    fprintf(f,"\n  Orb=");
	    j = 0;
	    for (i = 0; i < n; ++i)
	    {
		fprintf(f," %d",sh->orbits[i]);
		if (sh->orbits[i] == i) ++j;
	    }
	    fprintf(f," [%d]\n",j);
	    if (sh->fixed < 0) break;
	}
    }
}

/************************************************************************/

void
grouporder(int *fix, int nfix, schreier *gp, permnode **ring, 
                           double *grpsize1, int *grpsize2, int n)
/* process the base like in getorbits(), then return the product of the
 * orbits along the base, using the largest orbit at the end if the
 * base is not complete.
*/
{
    schreier *sh;
    int i,j,k,fx;
    int *orb;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n,"grouporder");
#endif

    getorbits(fix,nfix,gp,ring,n);
    expandschreier(gp,ring,n);
    expandschreier(gp,ring,n);
    *grpsize1 = 1.0; *grpsize2 = 0;

    for (i = 0, sh = gp; i < nfix; ++i, sh = sh->next)
    {
	orb = sh->orbits;
	fx = orb[sh->fixed];
	k = 0;
	for (j = fx; j < n; ++j) if (orb[j] == fx) ++k;
	MULTIPLY(*grpsize1,*grpsize2,k);
    }

    orb = sh->orbits;
    k = 1;
    for (i = 0; i < n; ++i)
	if (orb[i] == i)
	    workperm[i] = 1;
	else
	{
	    ++workperm[orb[i]];
	    if (workperm[orb[i]] > k) k = workperm[orb[i]];
	}
    MULTIPLY(*grpsize1,*grpsize2,k);
}

/*****************************************************************************
*                                                                            *
*  schreier_check() checks that this file is compiled compatibly with the    *
*  given parameters.   If not, call exit(1).                                 *
*                                                                            *
*****************************************************************************/

void
schreier_check(int wordsize, int m, int n, int version)
{
        if (wordsize != WORDSIZE)
        {
            fprintf(ERRFILE,"Error: WORDSIZE mismatch in schreier.c\n");
            exit(1);
        }

#if MAXN
        if (m > MAXM)
        {
            fprintf(ERRFILE,"Error: MAXM inadequate in schreier.c\n");
            exit(1);
        }

        if (n > MAXN)
        {
            fprintf(ERRFILE,"Error: MAXN inadequate in schreier.c\n");
            exit(1);
        }
#endif

        if (version < NAUTYREQUIRED)
        {
            fprintf(ERRFILE,"Error: schreier.c version mismatch\n");
            exit(1);
        }
}

/************************************************************************/

void
schreier_freedyn(void)
{
#if !MAXN
    DYNFREE(workperm,workperm_sz);
    DYNFREE(workperm2,workperm2_sz);
    DYNFREE(workpermA,workpermA_sz);
    DYNFREE(workpermB,workpermB_sz);
    DYNFREE(workset,workset_sz);
    DYNFREE(workset2,workset2_sz);
#endif
    clearfreelists();
}
