/* naugroup.h

Procedures for handling groups found by nauty.
*/

#include "nauty.h"

typedef struct perm_struct
{
    struct perm_struct *ptr;   /* general-purpose pointer */
    int p[2];          /* extendable section */
} permrec;

typedef struct coset_struct
{
    int image;       /* image of fixed point */
    permrec *rep;    /* pointer to a representative */
} cosetrec;

typedef struct level_struct
{
    int fixedpt;       /* point that is fixed in this level */
    int orbitsize;     /* the size of the orbit containing fixedpt */
    permrec *gens;     /* pointer to list of generators */
    cosetrec *replist; /* array of orbitsize representatives */
} levelrec;

typedef struct group_struct
{
    int n;                   /* number of points */
    int numorbits;           /* number of orbits */
    int depth;               /* number of points in base */
    levelrec levelinfo[1];   /* extendable section */
} grouprec;

#ifdef __cplusplus
extern "C" {
#endif

extern void freepermrec(permrec*, int);
extern grouprec *groupptr(boolean);
extern permrec *newpermrec(int);
extern void groupautomproc(int,int*,int*,int,int,int);
extern void
   grouplevelproc(int*,int*,int,int*,statsblk*,int,int,int,int,int,int);
extern void makecosetreps(grouprec*);
extern int permcycles(int*,int,int*,boolean);
extern void allgroup(grouprec*,void(*)(int*,int));
extern int allgroup2(grouprec*,void(*)(int*,int,int*));
extern int allgroup3(grouprec*,void(*)(int*,int,int*,void*),void*);
extern void freegroup(grouprec*);

#ifdef __cplusplus
}
#endif
