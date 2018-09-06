/* schreier.h - Version 1.2 (January 2013) */

#ifndef  _SCHREIER_H_    /* only process this file once */
#define  _SCHREIER_H_

#include "nauty.h"
#include "naurng.h"

typedef struct permnodestruct
{
    struct permnodestruct *prev,*next;   /* prev&next in circular list */
    unsigned long refcount;              /* number of references */
    int nalloc;                          /* size of p[] in ints,
                                            <= 0 for a perm marker */
    int mark;                            /* a mark, 0 unless changed */
    int p[2];                            /* actual vector, extended to
                                            nalloc enties */
} permnode;

typedef struct schreierlevel
{
    struct schreierlevel *next;    /* down one level, if any */
    int fixed;                     /* fixed at next level, -1 if none */
                	/* Invariant: next=NULL => fixed = -1 */
    int nalloc;                    /* size of vec[] and orbits[] */
    permnode **vec;                /* vec[i]^pwr[i] is edge label, */
    int *pwr;                      /*  transitive closure maps i->fixed */
    int *orbits;                   /* vector of orbits */
    permnode *marker;              /* points to marker for this level */
} schreier;

#define SCHREIERFAILS 10 
  /* Default number of Schreier failures before giving up. */

#ifdef __cplusplus
extern "C" {
#endif

/* See separate file schreier.txt for a description of usage. */

extern void freeschreier(schreier **gp, permnode **gens);
extern void addpermutation(permnode **ring, int *p, int n); 
extern permnode *findpermutation(permnode *gens, int *p, int n);
extern boolean addgenerator(schreier **gp, permnode **gens, int *p, int n);
extern boolean
	condaddgenerator(schreier **gp, permnode **gens, int *p, int n);
extern boolean expandschreier(schreier *gp, permnode **gens, int n);
extern int *getorbits(int *fix, int nfix,
		 schreier *gp, permnode **gens, int n);
extern int getorbitsmin(int *fix, int nfix, schreier *gp, permnode **gens,
		 int **orbits, int *cell, int ncell, int n, boolean changed);
extern void pruneset(set *fixset, schreier *gp, permnode **gens,
		      set *x, int m, int n);
extern void newgroup(schreier **gp, permnode **gens, int n);
extern void schreier_freedyn(void);
extern int schreier_fails(int nfails);
extern void dumpschreier(FILE *f, schreier *gp, permnode *gens, int n);
extern int schreier_gens(permnode *gens);
extern void deleteunmarked(permnode **gens);
extern void grouporder(int *fix, int nfix,  schreier *gp, permnode **gens,
		double *grpsize1, int *grpsize2, int n);
extern void schreier_check(int wordsize, int m, int n, int version);

#ifdef __cplusplus
}
#endif

#endif  /* _SCHREIER_H_ */
