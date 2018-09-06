/* naugstrings.h : Write graph6 or sparse6 strings into array. */
/* Version 1.1, Jun 2015. */

#include "gtools.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void gtog6string(graph*,char**,int,int);
extern void gtos6string(graph*,char**,int,int);
extern void gtod6string(graph*,char**,int,int);
extern void sgtos6string(sparsegraph*,char**);
extern void sgtog6string(sparsegraph*,char**);
extern void sgtod6string(sparsegraph*,char**);
extern void gtois6string(graph*,graph*,char**,int,int);

#ifdef __cplusplus
}
#endif
