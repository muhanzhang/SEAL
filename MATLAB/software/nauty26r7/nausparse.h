/* nausparse.h : header file for sparse digraphs, nauty 2.6 */
/* This version allows only simple graphs with loops but
 * contains the data structures for weights on the edges
 * even though they aren't implemented yet. */

/*****************************************************************************
*                                                                            *
*   Copyright (1984-2016) Brendan McKay.  All rights reserved.               *
*   Subject to the waivers and disclaimers in nauty.h.                       *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       10-Nov-09 : removed types shortish and permutation                   *
*       20-May-10 : make some fields type size_t                             *
*       23-May-10 : add sparsenauty()                                        *
*        3-Jun-10 : add *_tr procedures used by Traces                       *
*       30-Jun-10 : add DEFAULTOPTIONS_SPARSEDIGRAPH()                       *
*       18-Aug-12 : fix SG_DECL initialization order                         *
*       18-Jan-13 : add usercanonproc to default options                     *
*       17-Dec-15 : add macros for weighted graphs                           *
*                                                                            *
*****************************************************************************/

#ifndef  _NAUSPARSE_H_    /* only process this file once */
#define  _NAUSPARSE_H_

#include "nauty.h"

#ifndef SG_WEIGHT
#define SG_WEIGHT int
#define SG_WEIGHT_FMT "%d"
#define SG_MINWEIGHT (-NAUTY_INFINITY)
#endif
typedef SG_WEIGHT sg_weight;

#define CHECK_SWG(sg,id) do { if ((sg)->w) { fprintf(stderr, \
 ">E procedure %s does not accept weighted graphs\n",id); exit(1); } } while (0)

typedef struct
{
    size_t nde;  /* Number of directed edges (loops contribute only 1) */
    size_t *v;   /* Array of indexes into e[*] */
    int nv;      /* Number of vertices */
    int *d;      /* Array with out-degree of each vertex */
    int *e;      /* Array to hold lists of neighbours */
    sg_weight *w;      /* Not implemented, should be NULL. */
    size_t vlen,dlen,elen,wlen;  /* Sizes of arrays in units of type */
} sparsegraph;


#define SG_VDE(sgp,vv,dd,ee) do { vv = ((sparsegraph*)(sgp))->v; \
  dd = ((sparsegraph*)(sgp))->d; ee = ((sparsegraph*)(sgp))->e; } while(0)
#define SWG_VDE(sgp,vv,dd,ee,ww) do { vv = ((sparsegraph*)(sgp))->v; \
  dd = ((sparsegraph*)(sgp))->d; ee = ((sparsegraph*)(sgp))->e; \
  ww = ((sparsegraph*)(sgp))->w; } while(0)
#define SG_INIT(sg) do { (sg).v = NULL; (sg).d = (sg).e = (sg).w = NULL; \
   (sg).vlen = (sg).dlen = (sg).elen = (sg).wlen = 0; } while(0)
#define SWG_INIT SG_INIT
#define SG_ALLOC(sg,nlen,ndelen,msg) do { \
   DYNALLOC1(size_t,(sg).v,(sg).vlen,nlen,msg); \
   DYNALLOC1(int,(sg).d,(sg).dlen,nlen,msg); \
   DYNALLOC1(int,(sg).e,(sg).elen,ndelen,msg); \
} while (0)
#define SWG_ALLOC(sg,nlen,ndelen,msg) do { \
   DYNALLOC1(size_t,(sg).v,(sg).vlen,nlen,msg); \
   DYNALLOC1(int,(sg).d,(sg).dlen,nlen,msg); \
   DYNALLOC1(int,(sg).e,(sg).elen,ndelen,msg); \
   DYNALLOC1(sg_weight,(sg).w,(sg).wlen,ndelen,msg); \
} while (0)
#define SG_FREE(sg) do { \
   DYNFREE((sg).v,(sg).vlen); \
   DYNFREE((sg).d,(sg).dlen); \
   DYNFREE((sg).e,(sg).elen); \
   if ((sg).w) DYNFREE((sg).w,(sg).wlen); } while (0)
#define SWG_FREE SG_FREE

#define SG_DECL(sg) sparsegraph sg = {0,NULL,0,NULL,NULL,NULL,0,0,0,0}
#define SWG_DECL SG_DECL

#define DEFAULTOPTIONS_SPARSEGRAPH(options) optionblk options = \
 {0,FALSE,FALSE,FALSE,TRUE,FALSE,CONSOLWIDTH, \
  NULL,NULL,NULL,NULL,NULL,NULL,NULL,100,0,1,0,&dispatch_sparse,FALSE,NULL}
#define DEFAULTOPTIONS_SPARSEDIGRAPH(options) optionblk options = \
 {0,TRUE,FALSE,FALSE,TRUE,FALSE,CONSOLWIDTH, \
  NULL,NULL,NULL,NULL,NULL,NULL,adjacencies_sg,100,0,999,0,&dispatch_sparse,FALSE,NULL}

#ifdef __cplusplus
extern "C" {
#endif

extern dispatchvec dispatch_sparse;

extern int targetcell_sg(graph*,int*,int*,int,int,boolean,int,int,int);
extern boolean cheapautom_sg(int*,int,boolean,int);
extern boolean isautom_sg(graph*,int*,boolean,int,int);
extern void refine_sg(graph*,int*,int*,int,int*,int*,set*,int*,int,int);
extern int testcanlab_sg(graph*,graph*,int*,int*,int,int);
extern void updatecan_sg(graph*,graph*,int*,int,int,int);
extern int testcanlab_tr(sparsegraph*,sparsegraph*,int*,int*,int*);
extern int comparelab_tr(sparsegraph*,int*,int*,int*,int*,int*,int*);
extern void updatecan_tr(sparsegraph*,sparsegraph*,int*,int*,int);
extern void init_sg(graph*,graph**,graph*,graph**,int*,int*,set*,
	                   struct optionstruct*,int*,int,int);
extern void cleanup_sg(graph*,graph**,graph*,graph**,int*,
           int*,optionblk*,statsblk*stats,int,int);
extern void nausparse_freedyn(void);
extern void nausparse_check(int,int,int,int);

extern sparsegraph *nauty_to_sg(graph*,sparsegraph*,int,int);
extern graph* sg_to_nauty(sparsegraph*,graph*,int,int*);
extern void sortlists_sg(sparsegraph*);
extern boolean aresame_sg(sparsegraph*,sparsegraph*);
extern void put_sg(FILE*,sparsegraph*,boolean,int);
extern sparsegraph *copy_sg(sparsegraph*,sparsegraph*);
extern void distvals(sparsegraph*,int,int*,int);

extern void sparsenauty(sparsegraph*g,int*,int*,int*,
                        optionblk*,statsblk*,sparsegraph*);

extern void
   adjacencies_sg(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void
   distances_sg(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);

#ifdef __cplusplus
}
#endif

#endif
