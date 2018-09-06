/*****************************************************************************
* This is the header file for versions 2.5 of nautinv.c.                     *
*                                                                            *
*   Copyright (1984-2013) Brendan McKay.  All rights reserved.               *
*   Subject to the waivers and disclaimers in nauty.h.                       *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       20-Apr-01 : initial creation out of naututil.h                       *
*       10-Nov-10 : remove types shortish and permutation                    *
*                                                                            *
*****************************************************************************/

#include "nauty.h"              /* which includes stdio.h */

#ifdef __cplusplus
extern "C" {
#endif

extern void adjacencies(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void adjtriang(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void cellcliq(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void cellfano(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void cellfano2(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void cellind(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void cellquads(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void cellquins(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void celltrips(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void cellstarts(int*,int,set*,int,int);
extern void cliques(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void distances(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void getbigcells(int*,int,int,int*,int*,int*,int);
extern void indsets(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void nautinv_check(int,int,int,int);
extern void nautinv_freedyn(void);
extern void quadruples(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void refinvar(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void setnbhd(graph*,int,int,set*,set*);
extern void triples(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
extern void twopaths(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
#ifdef __cplusplus
}
#endif
