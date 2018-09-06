/*****************************************************************************
*                                                                            *
* This is the header file for version 2.5 of nautaux.c.                      *
*                                                                            *
*   Copyright (1984-2013) Brendan McKay.  All rights reserved.               *
*   Subject to the waivers and disclaimers in nauty.h.                       *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       26-Apr-89 : initial creation for version 1.5.                        *
*       14-Oct-90 : renamed as version 1.6 (no changes to this file)         *
*        5-Jun-93 : renamed as version 1.7+ (no changes to this file)        *
*       18-Aug-93 : renamed as version 1.8 (no changes to this file)         *
*       17-Sep-93 : renamed as version 1.9 (no changes to this file)         *
*       19-Apr-95 : added prototype wrapper for C++                          *
*       16-Nov-00 : made changes listed in nauty.h.                          *
*        3-Nov-04 : declared nautaux_freedyn and nautaux_check               *
*       21-Aug-12 : remove nvector                                           *
*                                                                            *
*****************************************************************************/

#include "nauty.h"           /* which includes stdio.h */

#ifdef __cplusplus
extern "C" {
#endif

extern int component(graph*,int,set*,int,int);
extern boolean equitable(graph*,int*,int*,int,int,int);
extern long ptncode(graph*,int*,int*,int,int,int);
extern void nautaux_freedyn(void);
extern void nautaux_check(int,int,int,int);
#ifdef __cplusplus
}
#endif
