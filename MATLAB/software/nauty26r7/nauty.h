/**************************************************************************
*    This is the header file for Version 2.6 of nauty().                  *
*    nauty.h.  Generated from nauty-h.in by configure.
**************************************************************************/

#ifndef  _NAUTY_H_    /* only process this file once */
#define  _NAUTY_H_

/* The parts between the ==== lines are modified by configure when
creating nauty.h out of nauty-h.in.  If configure is not being used,
it is necessary to check they are correct.
====================================================================*/

/* Check whether various headers or options are available */
#define HAVE_UNISTD_H  1    /* <unistd.h> */
#define HAVE_SYSTYPES_H  1    /* <sys/types.h> */
#define HAVE_STDDEF_H  1     /* <stddef.h> */
#define HAVE_STDLIB_H  1    /* <stdlib.h> */
#define HAVE_STRING_H  1    /* <string.h> */
#define MALLOC_DEC 1  /* 1 = malloc() is declared in stdlib.h, */
				 /* 2 = in malloc.h, 0 = in neither place */
#define HAS_MATH_INF 1 /* INFINITY is defined in math.h or */
				 /* some system header likely to be used */
#define HAS_STDIO_UNLOCK 1  /* Whether there are getc_unlocked, */
		               /* putc_unlocked,flockfile and funlockfile*/

#define DEFAULT_WORDSIZE 0
#define SIZEOF_INT 4
#define SIZEOF_LONG 8
#define SIZEOF_LONG_LONG 8   /* 0 if nonexistent */
#define SIZEOF_INT128 16   /* 0 if nonexistent */

#define HAVE_CONST 1    /* compiler properly supports const */

/* Note that the following is only for running nauty in multiple threads
   and will slow it down a little otherwise. */
#define HAVE_TLS 0   /* have storage attribute for thread-local */
#define TLS_ATTR   /* if so, what it is.  if not, empty */

#define USE_ANSICONTROLS 0 
                          /* whether --enable-ansicontrols is used */

#define _FILE_OFFSET_BITS 0
#if _FILE_OFFSET_BITS == 64
#define _LARGEFILE_SOURCE
#else
#undef _FILE_OFFSET_BITS
#endif

/* Support of gcc extensions __builtin_clz, __builtin_clzl, __builtin_clzll */
#define HAVE_CLZ 1
#define HAVE_CLZL 1
#define HAVE_CLZLL 1

/* Support of gcc extensions
      __builtin_popcount, __builtin_popcountl, __builtin_popcountll
   Note that these may only be fast if the compiler switch -mpopcnt is used.

   Also the intrinsics
      _mm_popcnt_u32, _mm_popcnt_u64
   for the Intel compiler icc.  These need no compiler switch.
*/
#define HAVE_HWPOPCMT 1
#define HAVE_POPCNT 1
#define HAVE_POPCNTL 1
#define HAVE_POPCNTLL 1
#define HAVE_MMPOP32 0
#define HAVE_MMPOP64 0

/*==================================================================*/

/* The following line must be uncommented for compiling into Magma. */
/* #define NAUTY_IN_MAGMA  */

#ifdef NAUTY_IN_MAGMA
#include "defs.h"
#include "system.h"
#include "bs.h"
#define OLDEXTDEFS
#else
#include <stdio.h>
#define P_(x) x
#endif

#if defined(__unix) || defined(__unix__) || defined(unix)
#define SYS_UNIX
#endif

#if !HAVE_CONST
#define const
#endif

/*****************************************************************************
*                                                                            *
*    AUTHOR: Brendan D. McKay                                                *
*            Research School of Computer Science                             *
*            Australian National University                                  *
*            Canberra, ACT 2601, Australia                                   *
*            phone:  +61 2 6125 3845                                         *
*            email:  Brendan.McKay@anu.edu.au                                *
*                                                                            *
*  This software is subject to copyright as detailed in the file COPYRIGHT.  *
*                                                                            *
*   Reference manual:                                                        *
*     B. D. McKay and A. Piperno, nauty User's Guide (Version 2.5),          *
*         http://pallini.di.uniroma1.it                                      *
*         http://cs.anu.edu.au/~bdm/nauty/                                   *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       10-Nov-87 : final changes for version 1.2                            *
*        5-Dec-87 : renamed to version 1.3 (no changes to this file)         *
*       28-Sep-88 : added PC Turbo C support, making version 1.4             *
*       23-Mar-89 : changes for version 1.5 :                                *
*                   - reworked M==1 code                                     *
*                   - defined NAUTYVERSION string                            *
*                   - made NAUTYH_READ to allow this file to be read twice   *
*                   - added optional ANSI function prototypes                *
*                   - added validity check for WORDSIZE                      *
*                   - added new fields to optionblk structure                *
*                   - updated DEFAULTOPTIONS to add invariants fields        *
*                   - added (set*) cast to definition of GRAPHROW            *
*                   - added definition of ALLOCS and FREES                   *
*       25-Mar-89 : - added declaration of new function doref()              *
*                   - added UNION macro                                      *
*       29-Mar-89 : - reduced the default MAXN for small machines            *
*                   - removed OUTOFSPACE (no longer used)                    *
*                   - added SETDIFF and XOR macros                           *
*        2-Apr-89 : - extended statsblk structure                            *
*        4-Apr-89 : - added IS_* macros                                      *
*                   - added ERRFILE definition                               *
*                   - replaced statsblk.outofspace by statsblk.errstatus     *
*        5-Apr-89 : - deleted definition of np2vector (no longer used)       *
*                   - introduced EMPTYSET macro                              *
*       12-Apr-89 : - eliminated MARK, UNMARK and ISMARKED (no longer used)  *
*       18-Apr-89 : - added MTOOBIG and CANONGNIL                            *
*       12-May-89 : - made ISELEM1 and ISELEMENT return 0 or 1               *
*        2-Mar-90 : - added EXTPROC macro and used it                        *
*       12-Mar-90 : - added SYS_CRAY, with help from N. Sloane and A. Grosky *
*                   - added dummy groupopts field to optionblk               *
*                   - select some ANSI things if __STDC__ exists             *
*       20-Mar-90 : - changed default MAXN for Macintosh versions            *
*                   - created SYS_MACTHINK for Macintosh THINK compiler      *
*       27-Mar-90 : - split SYS_MSDOS into SYS_PCMS4 and SYS_PCMS5           *
*       13-Oct-90 : changes for version 1.6:                                 *
*                   - fix definition of setword for WORDSIZE==64             *
*       14-Oct-90 : - added SYS_APOLLO version to avoid compiler bug         *
*       15-Oct-90 : - improve detection of ANSI conformance                  *
*       17-Oct-90 : - changed temp name in EMPTYSET to avoid A/UX bug        *
*       16-Apr-91 : changes for version 1.7:                                 *
*                   - made version SYS_PCTURBO use free(), not cfree()       *
*        2-Sep-91 : - noted that SYS_PCMS5 also works for Quick C            *
*                   - moved MULTIPLY to here from nauty.c                    *
*       12-Jun-92 : - changed the top part of this comment                   *
*       27-Aug-92 : - added version SYS_IBMC, thanks to Ivo Duentsch         *
*        5-Jun-93 : - renamed to version 1.7+, only change in naututil.h     *
*       29-Jul-93 : changes for version 1.8:                                 *
*                   - fixed error in default 64-bit version of FIRSTBIT      *
*                     (not used in any version before ALPHA)                 *
*                   - installed ALPHA version (thanks to Gordon Royle)       *
*                   - defined ALLOCS,FREES for SYS_IBMC                      *
*        3-Sep-93 : - make calloc void* in ALPHA version                     *
*       17-Sep-93 : - renamed to version 1.9,                                *
*                        changed only dreadnaut.c and nautinv.c              *
*       24-Feb-94 : changes for version 1.10:                                *
*                   - added version SYS_AMIGAAZT, thanks to Carsten Saager   *
*                     (making 1.9+)                                          *
*       19-Apr-95 : - added prototype wrapper for C++,                       *
*                     thanks to Daniel Huson                                 *
*        5-Mar-96 : - added SYS_ALPHA32 version (32-bit setwords on Alpha)   *
*       13-Jul-96 : changes for version 2.0:                                 *
*                   - added dynamic allocation                               *
*                   - ERRFILE must be defined                                *
*                   - added FLIPELEM1 and FLIPELEMENT macros                 *
*       13-Aug-96 : - added SWCHUNK? macros                                  *
*                   - added TAKEBIT macro                                    *
*       28-Nov-96 : - include sys/types.h if not ANSI (tentative!)           *
*       24-Jan-97 : - and stdlib.h if ANSI                                   *
*                   - removed use of cfree() from UNIX variants              *
*       25-Jan-97 : - changed options.getcanon from boolean to int           *
*                     Backwards compatibility is ok, as boolean and int      *
*                     are the same.  Now getcanon=2 means to get the label   *
*                     and not care about the group.  Sometimes faster.       *
*        6-Feb-97 : - Put in #undef for FALSE and TRUE to cope with          *
*                     compilers that illegally predefine them.               *
*                   - declared nauty_null and nautil_null                    *
*        2-Jul-98 : - declared ALLBITS                                       *
*       21-Oct-98 : - allow WORDSIZE==64 using unsigned long long            *
*                   - added BIGNAUTY option for really big graphs            *
*       11-Dec-99 : - made bit, leftbit and bytecount static in each file    *
*        9-Jan-00 : - declared nauty_check() and nautil_check()              *
*       12-Feb-00 : - Used #error for compile-time checks                    *
*                   - Added DYNREALLOC                                       *
*        4-Mar-00 : - declared ALLMASK(n)                                    *
*       27-May-00 : - declared CONDYNFREE                                    *
*       28-May-00 : - declared nautil_freedyn()                              *
*       16-Aug-00 : - added OLDNAUTY and changed canonical labelling         *
*       16-Nov-00 : - function prototypes are now default and unavoidable    *
*                   - removed UPROC, now assume all compilers know void      *
*                   - removed nvector, now just int (as it always was)       *
*                   - added extra parameter to targetcell()                  *
*                   - removed old versions which were only to skip around    *
*                     bugs that should have been long fixed:                 *
*                     SYS_APOLLO and SYS_VAXBSD.                             *
*                   - DEFAULTOPIONS now specifies no output                  *
*                   - Removed obsolete SYS_MACLSC version                    *
*       21-Apr-01 : - Added code to satisfy compilation into Magma.  This    *
*                       is activated by defining NAUTY_IN_MAGMA above.       *
*                   - The *_null routines no longer exist                    *
*                   - Default maxinvarlevel is now 1.  (This has no effect   *
*                        unless an invariant is specified.)                  *
*                   - Now labelorg has a concrete declaration in nautil.c    *
*                        and EXTDEFS is not needed                           *
*        5-May-01 : - NILFUNCTION, NILSET, NILGRAPH now obsolete.  Use NULL. *
*       11-Sep-01 : - setword is unsigned int in the event that UINT_MAX     *
*                     is defined and indicates it is big enough              *
*       17-Oct-01 : - major rewrite for 2.1.  SYS_* variables gone!          *
*                     Some modernity assumed, eg size_t                      *
*        8-Aug-02 : - removed MAKEEMPTY  (use EMPTYSET instead)              *
*                   - deleted OLDNAUTY everywhere                            *
*       27-Aug-02 : - converted to use autoconf.  Now the original of this   *
*                     file is nauty-h.in. Run configure to make nauty.h.     *
*       20-Dec-02 : - increased INFINITY                                     *
*                     some reorganization to please Magma                    *
*                   - declared nauty_freedyn()                               *
*       17-Nov-03 : - renamed INFINITY to NAUTY_INFINITY                     *
*       29-May-04 : - added definition of SETWORD_FORMAT                     *
*       14-Sep-04 : - extended prototypes even to recursive functions        *
*       16-Oct-04 : - added DEFAULTOPTIONS_GRAPH                             *
*       24-Oct-04 : Starting 2.3                                             *
*                   - remove register declarations as modern compilers       *
*                     tend to find them a nuisance                           *
*                   - Don't define the obsolete symbol INFINITY if it is     *
*                     defined already                                        *
*       17-Nov-04 : - make 6 counters in statsblk unsigned long              *
*       17-Jan-04 : - add init() and cleanup() to dispatchvec                *
*       12-Nov-05 : - Changed NAUTY_INFINITY to 2^30+2 in BIGNAUTY case      *
*       22-Nov-06 : Starting 2.4                                             *
*                   - removed usertcellproc from options                     *
*                     changed bestcell to targetcell in dispatch vector      *
*                     declare targetcell and maketargetcell                  *
*       29-Nov-06 : - add extraoptions to optionblk                          *
*                   - add declarations of extra_autom and extra_level        *
*       10-Dec-06 : - BIGNAUTY is gone!  Now permutation=shortish=int.       *
*                     NAUTY_INFINITY only depends on whether sizeof(int)=2.  *
*       27-Jun-08 : - define nauty_counter and LONG_LONG_COUNTERS            *
*       30-Jun-08 : - declare version 2.4                                    *
*        8-Nov-09 : - final release of version 2.4;                          *
*       10-Nov-10 : Starting 2.5                                             *
*                   - declare shortish and permutation obsolete, now int     *
*       14-Nov-10 : - SETWORDSNEEDED(n)                                      *
*       23-May-10 : - declare densenauty()                                   *
*       29-Jun-10 : - add PRINT_COUNTER(f,x)                                 *
*                   - add DEFAULTOPTIONS_DIGRAPH()                           *
*       27-Mar-11 : - declare writegroupsize()                               *
*       14-Jan-12 : - add HAVE_TLS and TLS_ATTR                              *
*       21-Feb-12 : - add ENABLE_ANSI                                        *
*       18-Mar-12 : - add COUNTER_FMT                                        *
*       18-Aug-12 : - add ADDONEARC, ADDONEEDGE, EMPTYGRAPH                  *
*       29-Aug-12 : - add CLZ macros and FIRSTBITNZ                          *
*       19-Oct-12 : - add DEFAULT_WORDSIZE                                   *
*        3-Jan-12 : Released 2.5rc1                                          *
*       18-Jan-12 : Froze 2.5                                                *
*       18-Jan-12 : - add NAUABORTED and NAUKILLED                           *
*                   - add nauty_kill_request                                 *
*                   - add usercanonproc                                      *
*        1-Oct-15 : - add COUNTER_FMT_RAW                                    *
*       10-Jan-16 : - defined POPCOUNTMAC, optionally use popcnt             *
*                   - remove SYS_CRAY, let's hope it is long obsolete        *
*                   - add Intel popcount intrinsics for icc                  *
*       12-Jan-16 : - DYNFREE and CONDYNFREE now set the pointer to NULL     *
*       16-Jan-16 : - Change NAUTY_INFINITY to 2 billion + 2                 *
*                 : Froze 2.6                                                *
*       12-Mar-16 : - Add const to alloc_error()                             *
*                                                                            *
* ++++++ This file is automatically generated, don't edit it by hand! ++++++
*                                                                            *
*****************************************************************************/

/*****************************************************************************
*                                                                            *
*   16-bit, 32-bit and 64-bit versions can be selected by defining WORDSIZE. *
*   The largest graph that can be handled has MAXN vertices.                 *
*   Both WORDSIZE and MAXN can be defined on the command line.               *
*   WORDSIZE must be 16, 32 or 64; MAXN must be <= NAUTY_INFINITY-2;         *
*                                                                            *
*   With a very slight loss of efficiency (depending on platform), nauty     *
*   can be compiled to dynamically allocate arrays.  Predefine MAXN=0 to     *
*   achieve this effect, which is default behaviour from version 2.0.        *
*   In that case, graphs of size up to NAUTY_INFINITY-2 can be handled       *
*   if the the memory is available.                                          *
*                                                                            *
*   If only very small graphs need to be processed, use MAXN<=WORDSIZE       *
*   since this causes substantial code optimizations.                        *
*                                                                            *
*   Conventions and Assumptions:                                             *
*                                                                            *
*    A 'setword' is the chunk of memory that is occupied by one part of      *
*    a set.  This is assumed to be >= WORDSIZE bits in size.                 *
*                                                                            *
*    The rightmost (loworder) WORDSIZE bits of setwords are numbered         *
*    0..WORDSIZE-1, left to right.  It is necessary that the 2^WORDSIZE      *
*    setwords with the other bits zero are totally ordered under <,=,>.      *
*    This needs care on a 1's-complement machine.                            *
*                                                                            *
*    The int variables m and n have consistent meanings throughout.          *
*    Graphs have n vertices always, and sets have m setwords always.         *
*                                                                            *
*    A 'set' consists of m contiguous setwords, whose bits are numbered      *
*    0,1,2,... from left (high-order) to right (low-order), using only       *
*    the rightmost WORDSIZE bits of each setword.  It is used to             *
*    represent a subset of {0,1,...,n-1} in the usual way - bit number x     *
*    is 1 iff x is in the subset.  Bits numbered n or greater, and           *
*    unnumbered bits, are assumed permanently zero.                          *
*                                                                            *
*    A 'graph' consists of n contiguous sets.  The i-th set represents       *
*    the vertices adjacent to vertex i, for i = 0,1,...,n-1.                 *
*                                                                            *
*    A 'permutation' is an array of n ints repesenting a permutation of      *
*    the set {0,1,...,n-1}.  The value of the i-th entry is the number to    *
*    which i is mapped.                                                      *
*                                                                            *
*    If g is a graph and p is a permutation, then g^p is the graph in        *
*    which vertex i is adjacent to vertex j iff vertex p[i] is adjacent      *
*    to vertex p[j] in g.                                                    *
*                                                                            *
*    A partition nest is represented by a pair (lab,ptn), where lab and ptn  *
*    are int arrays.  The "partition at level x" is the partition whose      *
*    cells are {lab[i],lab[i+1],...,lab[j]}, where [i,j] is a maximal        *
*    subinterval of [0,n-1] such that ptn[k] > x for i <= k < j and          *
*    ptn[j] <= x.  The partition at level 0 is given to nauty by the user.   *
*    This is  refined for the root of the tree, which has level 1.           *
*                                                                            *
*****************************************************************************/

#ifndef NAUTY_IN_MAGMA
#if HAVE_SYSTYPES_H
#include <sys/types.h>
#endif
#if HAVE_UNISTD_H
#include <unistd.h>
#endif
#if HAVE_STDDEF_H
#include <stddef.h>
#endif
#if HAVE_STDLIB_H
#include <stdlib.h>
#endif
#if HAVE_STRING_H
#include <string.h>
#else
#include <strings.h>
#endif
#endif

/* WORDSIZE is the number of set elements per setword (16, 32 or 64).
   WORDSIZE and setword are defined as follows:

   DEFAULT_WORDSIZE is usually 0 but is set by the configure script
   to NN if --enable-wordsize=NN is used, where NN is 16, 32 or 64.

   If WORDSIZE is not defined, but DEFAULT_WORDSIZE > 0, then set
      WORDSIZE to the same value as DEFAULT_WORDSIZE.
   If WORDSIZE is so far undefined, use 32 unless longs have more 
      than 32 bits, in which case use 64.
   Define setword thus:
      WORDSIZE==16 : unsigned short
      WORDSIZE==32 : unsigned int unless it is too small,
			in which case unsigned long
      WORDSIZE==64 : the first of unsigned int, unsigned long,
                      unsigned long long, which is large enough.
*/

#ifdef NAUTY_IN_MAGMA
#undef WORDSIZE
#define WORDSIZE WORDBITS
#endif

#ifndef WORDSIZE
#if DEFAULT_WORDSIZE > 0
#define WORDSIZE DEFAULT_WORDSIZE
#endif
#endif

#ifdef WORDSIZE

#if  (WORDSIZE != 16) && (WORDSIZE != 32) && (WORDSIZE != 64)
 #error "WORDSIZE must be 16, 32 or 64"
#endif

#else  /* WORDSIZE undefined */

#if SIZEOF_LONG>4
#define WORDSIZE 64
#else
#define WORDSIZE 32
#endif

#endif  /* WORDSIZE */

#ifdef NAUTY_IN_MAGMA
typedef t_uint setword;
#define SETWORD_INT  /* Don't assume this is correct in Magma. */

#else /* NAUTY_IN_MAGMA */

#if WORDSIZE==16
typedef unsigned short setword;
#define SETWORD_SHORT
#endif

#if WORDSIZE==32
#if SIZEOF_INT>=4
typedef unsigned int setword;
#define SETWORD_INT
#else
typedef unsigned long setword;
#define SETWORD_LONG
#endif
#endif

#if WORDSIZE==64
#if SIZEOF_INT>=8
typedef unsigned int setword;
#define SETWORD_INT
#else
#if SIZEOF_LONG>=8
typedef unsigned long setword;
#define SETWORD_LONG
#else
typedef unsigned long long setword;
#define SETWORD_LONGLONG
#endif
#endif
#endif

#endif /* NAUTY_IN_MAGMA else */

#if SIZEOF_LONG_LONG>=8 && SIZEOF_LONG==4
typedef unsigned long long nauty_counter;
#define LONG_LONG_COUNTERS 1
#define COUNTER_FMT "%llu"
#define COUNTER_FMT_RAW "llu"
#else
typedef unsigned long nauty_counter;
#define LONG_LONG_COUNTERS 0
#define COUNTER_FMT "%lu"
#define COUNTER_FMT_RAW "lu"
#endif
#define PRINT_COUNTER(f,x) fprintf(f,COUNTER_FMT,x)

#define NAUTYVERSIONID (26040+HAVE_TLS)  /* 10000*version + HAVE_TLS */
#define NAUTYREQUIRED NAUTYVERSIONID  /* Minimum compatible version */

#if WORDSIZE==16
#define NAUTYVERSION "2.6 (16 bits)"
#endif
#if WORDSIZE==32
#define NAUTYVERSION "2.6 (32 bits)"
#endif
#if WORDSIZE==64
#define NAUTYVERSION "2.6 (64 bits)"
#endif

#ifndef  MAXN  /* maximum allowed n value; use 0 for dynamic sizing. */
#define MAXN 0
#define MAXM 0
#else
#define MAXM ((MAXN+WORDSIZE-1)/WORDSIZE)  /* max setwords in a set */
#endif  /* MAXN */

/* Starting at version 2.2, set operations work for all set sizes unless
   ONE_WORD_SETS is defined.  In the latter case, if MAXM=1, set ops
   work only for single-setword sets.  In any case, macro versions
   ending with 1 work for single-setword sets and versions ending with
   0 work for all set sizes.
*/

#if  WORDSIZE==16
#define SETWD(pos) ((pos)>>4)  /* number of setword containing bit pos */
#define SETBT(pos) ((pos)&0xF) /* position within setword of bit pos */
#define TIMESWORDSIZE(w) ((w)<<4)
#define SETWORDSNEEDED(n) ((((n)-1)>>4)+1)  /* setwords needed for n bits */
#endif

#if  WORDSIZE==32
#define SETWD(pos) ((pos)>>5)
#define SETBT(pos) ((pos)&0x1F)
#define TIMESWORDSIZE(w) ((w)<<5)
#define SETWORDSNEEDED(n) ((((n)-1)>>5)+1)
#endif

#if  WORDSIZE==64
#define SETWD(pos) ((pos)>>6)
#define SETBT(pos) ((pos)&0x3F)
#define TIMESWORDSIZE(w) ((w)<<6)    /* w*WORDSIZE */
#define SETWORDSNEEDED(n) ((((n)-1)>>6)+1)
#endif

#ifdef NAUTY_IN_MAGMA
#define BITT bs_bit
#else
#define BITT bit
#endif

#define ADDELEMENT1(setadd,pos)  (*(setadd) |= BITT[pos])
#define DELELEMENT1(setadd,pos)  (*(setadd) &= ~BITT[pos])
#define FLIPELEMENT1(setadd,pos) (*(setadd) ^= BITT[pos])
#define ISELEMENT1(setadd,pos)   ((*(setadd) & BITT[pos]) != 0)
#define EMPTYSET1(setadd,m)   *(setadd) = 0;
#define GRAPHROW1(g,v,m) ((set*)(g)+(v))
#define ADDONEARC1(g,v,w,m) (g)[v] |= BITT[w]
#define ADDONEEDGE1(g,v,w,m) { ADDONEARC1(g,v,w,m); ADDONEARC1(g,w,v,m); }
#define EMPTYGRAPH1(g,m,n) EMPTYSET0(g,n)  /* really EMPTYSET0 */

#define ADDELEMENT0(setadd,pos)  ((setadd)[SETWD(pos)] |= BITT[SETBT(pos)])
#define DELELEMENT0(setadd,pos)  ((setadd)[SETWD(pos)] &= ~BITT[SETBT(pos)])
#define FLIPELEMENT0(setadd,pos) ((setadd)[SETWD(pos)] ^= BITT[SETBT(pos)])
#define ISELEMENT0(setadd,pos) (((setadd)[SETWD(pos)] & BITT[SETBT(pos)]) != 0)
#define EMPTYSET0(setadd,m) \
    {setword *es; \
    for (es = (setword*)(setadd)+(m); --es >= (setword*)(setadd);) *es=0;}
#define GRAPHROW0(g,v,m) ((set*)(g) + (m)*(size_t)(v))
#define ADDONEARC0(g,v,w,m) ADDELEMENT0(GRAPHROW0(g,v,m),w)
#define ADDONEEDGE0(g,v,w,m) { ADDONEARC0(g,v,w,m); ADDONEARC0(g,w,v,m); }
#define EMPTYGRAPH0(g,m,n) EMPTYSET0(g,(m)*(size_t)(n))

#if  (MAXM==1) && defined(ONE_WORD_SETS)
#define ADDELEMENT ADDELEMENT1
#define DELELEMENT DELELEMENT1
#define FLIPELEMENT FLIPELEMENT1
#define ISELEMENT ISELEMENT1
#define EMPTYSET EMPTYSET1
#define GRAPHROW GRAPHROW1
#define ADDONEARC ADDONEARC1
#define ADDONEEDGE ADDONEEDGE1
#define EMPTYGRAPH EMPTYGRAPH1
#else
#define ADDELEMENT ADDELEMENT0
#define DELELEMENT DELELEMENT0
#define FLIPELEMENT FLIPELEMENT0
#define ISELEMENT ISELEMENT0
#define EMPTYSET EMPTYSET0
#define GRAPHROW GRAPHROW0
#define ADDONEARC ADDONEARC0
#define ADDONEEDGE ADDONEEDGE0
#define EMPTYGRAPH EMPTYGRAPH0
#endif


#ifdef NAUTY_IN_MAGMA
#undef EMPTYSET
#define EMPTYSET(setadd,m) {t_int _i; bsp_makeempty(setadd,m,_i);}
#endif

#define NOTSUBSET(word1,word2) ((word1) & ~(word2))  /* test if the 1-bits
                    in setword word1 do not form a subset of those in word2  */
#define INTERSECT(word1,word2) ((word1) &= (word2))  /* AND word2 into word1 */
#define UNION(word1,word2)     ((word1) |= (word2))  /* OR word2 into word1 */
#define SETDIFF(word1,word2)   ((word1) &= ~(word2)) /* - word2 into word1 */
#define XOR(word1,word2)       ((word1) ^= (word2))  /* XOR word2 into word1 */
#define ZAPBIT(word,x) ((word) &= ~BITT[x])  /* delete bit x in setword */
#define TAKEBIT(iw,w) {(iw) = FIRSTBITNZ(w); (w) ^= BITT[iw];}

#ifdef SETWORD_LONGLONG
#define MSK3232 0xFFFFFFFF00000000ULL
#define MSK1648 0xFFFF000000000000ULL
#define MSK0856 0xFF00000000000000ULL
#define MSK1632 0x0000FFFF00000000ULL
#define MSK0840     0xFF0000000000ULL
#define MSK1616         0xFFFF0000ULL 
#define MSK0824         0xFF000000ULL 
#define MSK0808             0xFF00ULL 
#define MSK63C  0x7FFFFFFFFFFFFFFFULL
#define MSK31C          0x7FFFFFFFULL
#define MSK15C              0x7FFFULL
#define MSK64   0xFFFFFFFFFFFFFFFFULL
#define MSK32           0xFFFFFFFFULL
#define MSK16               0xFFFFULL
#define MSK8                  0xFFULL
#endif

#ifdef SETWORD_LONG
#define MSK3232 0xFFFFFFFF00000000UL
#define MSK1648 0xFFFF000000000000UL
#define MSK0856 0xFF00000000000000UL
#define MSK1632 0x0000FFFF00000000UL
#define MSK0840     0xFF0000000000UL
#define MSK1616         0xFFFF0000UL 
#define MSK0824         0xFF000000UL 
#define MSK0808             0xFF00UL 
#define MSK63C  0x7FFFFFFFFFFFFFFFUL
#define MSK31C          0x7FFFFFFFUL
#define MSK15C              0x7FFFUL
#define MSK64   0xFFFFFFFFFFFFFFFFUL
#define MSK32           0xFFFFFFFFUL
#define MSK16               0xFFFFUL
#define MSK8                  0xFFUL
#endif

#if defined(SETWORD_INT) || defined(SETWORD_SHORT)
#define MSK3232 0xFFFFFFFF00000000U
#define MSK1648 0xFFFF000000000000U
#define MSK0856 0xFF00000000000000U
#define MSK1632 0x0000FFFF00000000U
#define MSK0840     0xFF0000000000U
#define MSK1616         0xFFFF0000U 
#define MSK0824         0xFF000000U
#define MSK0808             0xFF00U 
#define MSK63C  0x7FFFFFFFFFFFFFFFU
#define MSK31C          0x7FFFFFFFU
#define MSK15C              0x7FFFU
#define MSK64   0xFFFFFFFFFFFFFFFFU
#define MSK32           0xFFFFFFFFU
#define MSK16               0xFFFFU
#define MSK8                  0xFFU
#endif

#if defined(SETWORD_LONGLONG)
#if WORDSIZE==16
#define SETWORD_FORMAT "%04llx"
#endif
#if WORDSIZE==32
#define SETWORD_FORMAT "%08llx"
#endif
#if WORDSIZE==64
#define SETWORD_FORMAT "%16llx"
#endif
#endif

#if defined(SETWORD_LONG)
#if WORDSIZE==16
#define SETWORD_FORMAT "%04lx"
#endif
#if WORDSIZE==32
#define SETWORD_FORMAT "%08lx"
#endif
#if WORDSIZE==64
#define SETWORD_FORMAT "%16lx"
#endif
#endif

#if defined(SETWORD_INT)
#if WORDSIZE==16
#define SETWORD_FORMAT "%04x"
#endif
#if WORDSIZE==32
#define SETWORD_FORMAT "%08x"
#endif
#if WORDSIZE==64
#define SETWORD_FORMAT "%16x"
#endif
#endif

#if defined(SETWORD_SHORT)
#if WORDSIZE==16
#define SETWORD_FORMAT "%04hx"
#endif
#if WORDSIZE==32
#define SETWORD_FORMAT "%08hx"
#endif
#if WORDSIZE==64
#define SETWORD_FORMAT "%16hx"
#endif
#endif

/* POPCOUNT(x) = number of 1-bits in a setword x
   POPCOUNTMAC(x) = Macro version of POPCOUNT
   FIRSTBIT(x) = number of first 1-bit in non-zero setword (0..WORDSIZE-1)
                   or WORDSIZE if x == 0
   FIRSTBITNZ(x) = as FIRSTBIT(x) but assumes x is not zero
   BITMASK(x)  = setword whose rightmost WORDSIZE-x-1 (numbered) bits
                 are 1 and the rest 0 (0 <= x < WORDSIZE)
                 (I.e., bits 0..x are unselected and the rest selected.)
   ALLBITS     = all (numbered) bits in a setword  */

#if  WORDSIZE==64
#define POPCOUNTMAC(x) (bytecount[(x)>>56 & 0xFF] + bytecount[(x)>>48 & 0xFF] \
                   + bytecount[(x)>>40 & 0xFF] + bytecount[(x)>>32 & 0xFF] \
                   + bytecount[(x)>>24 & 0xFF] + bytecount[(x)>>16 & 0xFF] \
                   + bytecount[(x)>>8 & 0xFF]  + bytecount[(x) & 0xFF])
#define FIRSTBIT(x) ((x) & MSK3232 ? \
                       (x) &   MSK1648 ? \
                         (x) & MSK0856 ? \
                         0+leftbit[((x)>>56) & MSK8] : \
                         8+leftbit[(x)>>48] \
                       : (x) & MSK0840 ? \
                         16+leftbit[(x)>>40] : \
                         24+leftbit[(x)>>32] \
                     : (x) & MSK1616 ? \
                         (x) & MSK0824 ? \
                         32+leftbit[(x)>>24] : \
                         40+leftbit[(x)>>16] \
                       : (x) & MSK0808 ? \
                         48+leftbit[(x)>>8] : \
                         56+leftbit[x])
#define BITMASK(x)  (MSK63C >> (x))
#define ALLBITS  MSK64
#define SWCHUNK0(w) ((long)((w)>>48)&0xFFFFL)
#define SWCHUNK1(w) ((long)((w)>>32)&0xFFFFL)
#define SWCHUNK2(w) ((long)((w)>>16)&0xFFFFL)
#define SWCHUNK3(w) ((long)(w)&0xFFFFL)
#endif

#if  WORDSIZE==32
#define POPCOUNTMAC(x) (bytecount[(x)>>24 & 0xFF] + bytecount[(x)>>16 & 0xFF] \
                        + bytecount[(x)>>8 & 0xFF] + bytecount[(x) & 0xFF])
#define FIRSTBIT(x) ((x) & MSK1616 ? ((x) & MSK0824 ? \
                     leftbit[((x)>>24) & MSK8] : 8+leftbit[(x)>>16]) \
                    : ((x) & MSK0808 ? 16+leftbit[(x)>>8] : 24+leftbit[x]))
#define BITMASK(x)  (MSK31C >> (x))
#define ALLBITS  MSK32
#define SWCHUNK0(w) ((long)((w)>>16)&0xFFFFL)
#define SWCHUNK1(w) ((long)(w)&0xFFFFL)
#endif

#if  WORDSIZE==16
#define POPCOUNTMAC(x) (bytecount[(x)>>8 & 0xFF] + bytecount[(x) & 0xFF])
#define FIRSTBIT(x) ((x) & MSK0808 ? leftbit[((x)>>8) & MSK8] : 8+leftbit[x])
#define BITMASK(x)  (MSK15C >> (x))
#define ALLBITS  MSK16
#define SWCHUNK0(w) ((long)(w)&0xFFFFL)
#endif

#if defined(SETWORD_LONGLONG) && HAVE_CLZLL
#undef FIRSTBIT
#undef FIRSTBITNZ
#define FIRSTBITNZ(x) __builtin_clzll(x)
#define FIRSTBIT(x) ((x) ? FIRSTBITNZ(x) : WORDSIZE)
#endif
#if defined(SETWORD_LONG) && HAVE_CLZL
#undef FIRSTBIT
#undef FIRSTBITNZ
#define FIRSTBITNZ(x) __builtin_clzl(x)
#define FIRSTBIT(x) ((x) ? FIRSTBITNZ(x) : WORDSIZE)
#endif
#if defined(SETWORD_INT) && HAVE_CLZ
#undef FIRSTBIT
#undef FIRSTBITNZ
#define FIRSTBITNZ(x) __builtin_clz(x)
#define FIRSTBIT(x) ((x) ? FIRSTBITNZ(x) : WORDSIZE)
#endif

#ifndef FIRSTBITNZ
#define FIRSTBITNZ FIRSTBIT
#endif

/* Use popcount instructions if available */
#define POPCOUNT POPCOUNTMAC

#ifdef __INTEL_COMPILER
#include <nmmintrin.h>
#if WORDSIZE==64 && HAVE_MMPOP64
#undef POPCOUNT
#define POPCOUNT(x) _mm_popcnt_u64(x)
#endif
#if WORDSIZE==32 && HAVE_MMPOP32
#undef POPCOUNT
#define POPCOUNT(x) _mm_popcnt_u32(x)
#endif

#else

/* Note that, unlike icc, gcc will not use the POPCNT instruction without
   permission, in which case it defines __POPCNT__ .
*/
#ifdef __POPCNT__
#if defined(SETWORD_LONGLONG) && HAVE_POPCNTLL
#undef POPCOUNT
#define POPCOUNT(x) __builtin_popcountll(x)
#endif
#if defined(SETWORD_LONG) && HAVE_POPCNTL
#undef POPCOUNT
#define POPCOUNT(x) __builtin_popcountl(x)
#endif
#if defined(SETWORD_INT) && HAVE_POPCNT
#undef POPCOUNT
#define POPCOUNT(x) __builtin_popcount(x)
#endif
#endif
#endif

#ifdef NAUTY_IN_MAGMA
#undef POPCOUNT
#undef FIRSTBIT
#undef BITMASK
#define POPCOUNT(x) bs_popcount(x)
#define FIRSTBIT(x) bs_firstbit(x)
#define BITMASK(x)  bs_bitmask(x)
#endif

#define ALLMASK(n) ((n)?~BITMASK((n)-1):(setword)0)  /* First n bits */

    /* various constants: */
#undef FALSE
#undef TRUE
#define FALSE    0
#define TRUE     1

#if SIZEOF_INT>=4
#define NAUTY_INFINITY 2000000002  /* Max graph size is 2 billion */
#else
#define NAUTY_INFINITY 0x7FFF
#endif

/* The following four types are obsolete, use int in new code. */
typedef int shortish;
typedef shortish permutation;
typedef int nvector,np2vector; 

/* For backward compatibility: */
#if !HAS_MATH_INF && !defined(INFINITY)
#define INFINITY NAUTY_INFINITY
#endif

#if MAXN > NAUTY_INFINITY-2
 #error MAXN must be at most NAUTY_INFINITY-2
#endif

    /* typedefs for sets, graphs, permutations, etc.: */

typedef int boolean;    /* boolean MUST be the same as int */

#define UPROC void      /* obsolete */

typedef setword set,graph;
#ifdef NAUTY_IN_MAGMA
typedef graph nauty_graph;
typedef set nauty_set;
#endif

typedef struct
{
    double grpsize1;        /* size of group is */
    int grpsize2;           /*    grpsize1 * 10^grpsize2 */
#define groupsize1 grpsize1     /* for backwards compatibility */
#define groupsize2 grpsize2
    int numorbits;          /* number of orbits in group */
    int numgenerators;      /* number of generators found */
    int errstatus;          /* if non-zero : an error code */
#define outofspace errstatus;   /* for backwards compatibility */
    unsigned long numnodes;      /* total number of nodes */
    unsigned long numbadleaves;  /* number of leaves of no use */
    int maxlevel;                /* maximum depth of search */
    unsigned long tctotal;       /* total size of all target cells */
    unsigned long canupdates;    /* number of updates of best label */
    unsigned long invapplics;    /* number of applications of invarproc */
    unsigned long invsuccesses;  /* number of successful uses of invarproc() */
    int invarsuclevel;      /* least level where invarproc worked */
} statsblk;

/* codes for errstatus field (see nauty.c for more accurate descriptions): */
/* 0 is normal - no error */
#define NTOOBIG      1      /* n > MAXN or n > WORDSIZE*m */
#define MTOOBIG      2      /* m > MAXM */
#define CANONGNIL    3      /* canong = NULL, but getcanon = TRUE */
#define NAUABORTED   4      /* nauty is terminated early under program control */
#define NAUKILLED    5      /* nauty is terminated early by caught signal */

/* manipulation of real approximation to group size */
#define MULTIPLY(s1,s2,i) if ((s1 *= i) >= 1e10) {s1 /= 1e10; s2 += 10;}

struct optionstruct;  /* incomplete definition */

typedef struct
{
    boolean (*isautom)        /* test for automorphism */
            (graph*,int*,boolean,int,int);
    int     (*testcanlab)     /* test for better labelling */
            (graph*,graph*,int*,int*,int,int);
    void    (*updatecan)      /* update canonical object */
            (graph*,graph*,int*,int,int,int);
    void    (*refine)         /* refine partition */
            (graph*,int*,int*,int,int*,int*,set*,int*,int,int);
    void    (*refine1)        /* refine partition, MAXM==1 */
            (graph*,int*,int*,int,int*,int*,set*,int*,int,int);
    boolean (*cheapautom)     /* test for easy automorphism */
            (int*,int,boolean,int);
    int     (*targetcell)     /* decide which cell to split */
            (graph*,int*,int*,int,int,boolean,int,int,int);
    void    (*freedyn)(void); /* free dynamic memory */
    void    (*check)          /* check compilation parameters */
            (int,int,int,int);
    void    (*init)(graph*,graph**,graph*,graph**,int*,int*,set*,
                   struct optionstruct*,int*,int,int);
    void    (*cleanup)(graph*,graph**,graph*,graph**,int*,int*,
                      struct optionstruct*,statsblk*,int,int);
} dispatchvec;

typedef struct optionstruct
{
    int getcanon;             /* make canong and canonlab? */
#define LABELONLY 2   /* new value UNIMPLEMENTED */
    boolean digraph;          /* multiple edges or loops? */
    boolean writeautoms;      /* write automorphisms? */
    boolean writemarkers;     /* write stats on pts fixed, etc.? */
    boolean defaultptn;       /* set lab,ptn,active for single cell? */
    boolean cartesian;        /* use cartesian rep for writing automs? */
    int linelength;           /* max chars/line (excl. '\n') for output */
    FILE *outfile;            /* file for output, if any */
    void (*userrefproc)       /* replacement for usual refine procedure */
         (graph*,int*,int*,int,int*,int*,set*,int*,int,int);
    void (*userautomproc)     /* procedure called for each automorphism */
         (int,int*,int*,int,int,int);
    void (*userlevelproc)     /* procedure called for each level */
         (int*,int*,int,int*,statsblk*,int,int,int,int,int,int);
    void (*usernodeproc)      /* procedure called for each node */
         (graph*,int*,int*,int,int,int,int,int,int);
    int  (*usercanonproc)     /* procedure called for better labellings */
         (graph*,int*,graph*,int,int,int,int);
    void (*invarproc)         /* procedure to compute vertex-invariant */
         (graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
    int tc_level;             /* max level for smart target cell choosing */
    int mininvarlevel;        /* min level for invariant computation */
    int maxinvarlevel;        /* max level for invariant computation */
    int invararg;             /* value passed to (*invarproc)() */
    dispatchvec *dispatch;    /* vector of object-specific routines */
    boolean schreier;         /* use random schreier method */
    void *extra_options;      /* arbitrary extra options */
#ifdef NAUTY_IN_MAGMA
    boolean print_stats;      /* CAYLEY specfic - GYM Sep 1990 */
    char *invarprocname;      /* Magma - no longer global sjc 1994 */
    int lab_h;                /* Magma - no longer global sjc 1994 */
    int ptn_h;                /* Magma - no longer global sjc 1994 */
    int orbitset_h;           /* Magma - no longer global sjc 1994 */
#endif
} optionblk;

#ifndef CONSOLWIDTH
#define CONSOLWIDTH 78
#endif

/* The following are obsolete.  Just use NULL. */
#define NILFUNCTION ((void(*)())NULL)      /* nil pointer to user-function */
#define NILSET      ((set*)NULL)           /* nil pointer to set */
#define NILGRAPH    ((graph*)NULL)         /* nil pointer to graph */

#define DEFAULTOPTIONS_GRAPH(options) optionblk options = \
 {0,FALSE,FALSE,FALSE,TRUE,FALSE,CONSOLWIDTH, \
  NULL,NULL,NULL,NULL,NULL,NULL,NULL,100,0,1,0,&dispatch_graph,FALSE,NULL}
#define DEFAULTOPTIONS_DIGRAPH(options) optionblk options = \
 {0,TRUE,FALSE,FALSE,TRUE,FALSE,CONSOLWIDTH, \
  NULL,NULL,NULL,NULL,NULL,NULL,adjacencies,100,0,999,0,&dispatch_graph,FALSE,NULL}

#ifndef DEFAULTOPTIONS
#define DEFAULTOPTIONS DEFAULTOPTIONS_GRAPH
#endif
#ifndef DEFAULTOPTIONS_DENSEGRAPH
#define DEFAULTOPTIONS_DENSEGRAPH DEFAULTOPTIONS_GRAPH
#endif

#ifdef NAUTY_IN_MAGMA
#define PUTC(c,f) io_putchar(c)
#else
#ifdef IS_JAVA
extern void javastream(FILE* f,char c);
#define PUTC(c,f) javastream(f,c)
#else
#define PUTC(c,f) putc(c,f)
#endif
#endif

/* We hope that malloc, free, realloc are declared either in <stdlib.h>
   or <malloc.h>.  Otherwise we will define them.  We also assume that
   size_t has been defined by the time we get to define malloc(). */
#ifndef NAUTY_IN_MAGMA
#if MALLOC_DEC==2
#include <malloc.h>
#endif
#if MALLOC_DEC==0
extern void *malloc(size_t);
extern void *realloc(void*,size_t);
extern void free(void*);
#endif
#endif

/* ALLOCS(x,y) should return a pointer (any pointer type) to x*y units of new
   storage, not necessarily initialised.  A "unit" of storage is defined by
   the sizeof operator.   x and y are integer values of type int or larger, 
   but x*y may well be too large for an int.  The macro should cast to the
   correct type for the call.  On failure, ALLOCS(x,y) should return a NULL 
   pointer.  FREES(p) should free storage previously allocated by ALLOCS, 
   where p is the value that ALLOCS returned. */

#ifdef NAUTY_IN_MAGMA
#define ALLOCS(x,y) mem_malloc((size_t)(x)*(size_t)(y))
#define REALLOCS(p,x) mem_realloc(p,(size_t)(x))
#define FREES(p) mem_free(p)
#else
#define ALLOCS(x,y) malloc((size_t)(x)*(size_t)(y))
#define REALLOCS(p,x) realloc(p,(size_t)(x)) 
#define FREES(p) free(p)
#endif

/* The following macros are used by nauty if MAXN=0.  They dynamically
   allocate arrays of size dependent on m or n.  For each array there
   should be two static variables:
     type *name;
     size_t name_sz;
   "name" will hold a pointer to an allocated array.  "name_sz" will hold
   the size of the allocated array in units of sizeof(type).  DYNALLSTAT
   declares both variables and initialises name_sz=0.  DYNALLOC1 and
   DYNALLOC2 test if there is enough space allocated, and if not free
   the existing space and allocate a bigger space.  The allocated space
   is not initialised.
   
   In the case of DYNALLOC1, the space is allocated using
       ALLOCS(sz,sizeof(type)).
   In the case of DYNALLOC2, the space is allocated using
       ALLOCS(sz1,sz2*sizeof(type)).

   DYNREALLOC is like DYNALLOC1 except that the old contents are copied
   into the new space. Availability of realloc() is assumed.

   DYNFREE frees any allocated array and sets name_sz back to 0.
   CONDYNFREE does the same, but only if name_sz exceeds some limit.
*/

#define DYNALLSTAT(type,name,name_sz) \
	static TLS_ATTR type *name; static TLS_ATTR size_t name_sz=0
#define DYNALLOC1(type,name,name_sz,sz,msg) \
 if ((size_t)(sz) > name_sz) \
 { if (name_sz) FREES(name); name_sz = (sz); \
 if ((name=(type*)ALLOCS(sz,sizeof(type))) == NULL) {alloc_error(msg);}}
#define DYNALLOC2(type,name,name_sz,sz1,sz2,msg) \
 if ((size_t)(sz1)*(size_t)(sz2) > name_sz) \
 { if (name_sz) FREES(name); name_sz = (size_t)(sz1)*(size_t)(sz2); \
 if ((name=(type*)ALLOCS((sz1),(sz2)*sizeof(type))) == NULL) \
 {alloc_error(msg);}}
#define DYNREALLOC(type,name,name_sz,sz,msg) \
 {if ((size_t)(sz) > name_sz) \
 { if ((name = (type*)REALLOCS(name,(sz)*sizeof(type))) == NULL) \
      {alloc_error(msg);} else name_sz = (sz);}}
#define DYNFREE(name,name_sz) \
  { if (name) FREES(name); name = NULL; name_sz = 0;}
#define CONDYNFREE(name,name_sz,minsz) \
 if (name_sz > (size_t)(minsz)) {DYNFREE(name,name_sz);}

/* File to write error messages to (used as first argument to fprintf()). */
#define ERRFILE stderr

/* Don't use OLDEXTDEFS, it is only still here for Magma. */
#ifdef OLDEXTDEFS   
#define EXTDEF_CLASS
#ifdef EXTDEFS
#define EXTDEF_TYPE 1
#else
#define EXTDEF_TYPE 2
#endif
#else
#define EXTDEF_CLASS static
#define EXTDEF_TYPE 2
#endif

extern int labelorg;   /* Declared in nautil.c */
extern volatile int nauty_kill_request;  /* Also declared in nautil.c */

#ifndef NAUTY_IN_MAGMA
  /* Things equivalent to bit, bytecount, leftbit are defined
     in bs.h for Magma. */
#if  EXTDEF_TYPE==1
extern setword bit[];
extern int bytecount[];
extern int leftbit[];

#else
    /* array giving setwords with single 1-bit */
#if  WORDSIZE==64
#ifdef SETWORD_LONGLONG
EXTDEF_CLASS const
setword bit[] = {01000000000000000000000LL,0400000000000000000000LL,
                 0200000000000000000000LL,0100000000000000000000LL,
                 040000000000000000000LL,020000000000000000000LL,
                 010000000000000000000LL,04000000000000000000LL,
                 02000000000000000000LL,01000000000000000000LL,
                 0400000000000000000LL,0200000000000000000LL,
                 0100000000000000000LL,040000000000000000LL,
                 020000000000000000LL,010000000000000000LL,
                 04000000000000000LL,02000000000000000LL,
                 01000000000000000LL,0400000000000000LL,0200000000000000LL,
                 0100000000000000LL,040000000000000LL,020000000000000LL,
                 010000000000000LL,04000000000000LL,02000000000000LL,
                 01000000000000LL,0400000000000LL,0200000000000LL,
		 0100000000000LL,040000000000LL,020000000000LL,010000000000LL,
		 04000000000LL,02000000000LL,01000000000LL,0400000000LL,
		 0200000000LL,0100000000LL,040000000LL,020000000LL,
		 010000000LL,04000000LL,02000000LL,01000000LL,0400000LL,
		 0200000LL,0100000LL,040000LL,020000LL,010000LL,04000LL,
                 02000LL,01000LL,0400LL,0200LL,0100LL,040LL,020LL,010LL,
		 04LL,02LL,01LL};
#else
EXTDEF_CLASS const
setword bit[] = {01000000000000000000000,0400000000000000000000,
                 0200000000000000000000,0100000000000000000000,
                 040000000000000000000,020000000000000000000,
                 010000000000000000000,04000000000000000000,
                 02000000000000000000,01000000000000000000,
                 0400000000000000000,0200000000000000000,
                 0100000000000000000,040000000000000000,020000000000000000,
                 010000000000000000,04000000000000000,02000000000000000,
                 01000000000000000,0400000000000000,0200000000000000,
                 0100000000000000,040000000000000,020000000000000,
                 010000000000000,04000000000000,02000000000000,
                 01000000000000,0400000000000,0200000000000,0100000000000,
                 040000000000,020000000000,010000000000,04000000000,
                 02000000000,01000000000,0400000000,0200000000,0100000000,
                 040000000,020000000,010000000,04000000,02000000,01000000,
                 0400000,0200000,0100000,040000,020000,010000,04000,
                 02000,01000,0400,0200,0100,040,020,010,04,02,01};
#endif
#endif

#if  WORDSIZE==32
EXTDEF_CLASS const
setword bit[] = {020000000000,010000000000,04000000000,02000000000,
                 01000000000,0400000000,0200000000,0100000000,040000000,
                 020000000,010000000,04000000,02000000,01000000,0400000,
                 0200000,0100000,040000,020000,010000,04000,02000,01000,
                 0400,0200,0100,040,020,010,04,02,01};
#endif

#if WORDSIZE==16
EXTDEF_CLASS const
setword bit[] = {0100000,040000,020000,010000,04000,02000,01000,0400,0200,
                 0100,040,020,010,04,02,01};
#endif

    /*  array giving number of 1-bits in bytes valued 0..255: */
EXTDEF_CLASS const
int bytecount[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                   4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};

    /* array giving position (1..7) of high-order 1-bit in byte: */
EXTDEF_CLASS const
int leftbit[] =   {8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,
                   3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
#endif  /* EXTDEFS */

#endif /* not NAUTY_IN_MAGMA */

#define ANSIPROT 1
#define EXTPROC(func,args) extern func args;     /* obsolete */

/* The following is for C++ programs that read nauty.h.  Compile nauty
   itself using C, not C++.  */

#ifdef __cplusplus
extern "C" {
#endif

extern void alloc_error(const char*);
extern void breakout(int*,int*,int,int,int,set*,int);
extern boolean cheapautom(int*,int,boolean,int);
extern void doref(graph*,int*,int*,int,int*,int*,int*,set*,int*,
  void(*)(graph*,int*,int*,int,int*,int*,set*,int*,int,int),
  void(*)(graph*,int*,int*,int,int,int,int*,int,boolean,int,int),
  int,int,int,boolean,int,int);
extern void extra_autom(int*,int);
extern void extra_level(int,int*,int*,int,int,int,int,int,int);
extern boolean isautom(graph*,int*,boolean,int,int);
extern dispatchvec dispatch_graph;
extern int itos(int,char*);
extern void fmperm(int*,set*,set*,int,int);
extern void fmptn(int*,int*,int,set*,set*,int,int);
extern void longprune(set*,set*,set*,set*,int);
extern void nauty(graph*,int*,int*,set*,int*,optionblk*,
                  statsblk*,set*,int,int,int,graph*);
extern void maketargetcell(graph*,int*,int*,int,set*,int*,int*,int,boolean,
           int,int (*)(graph*,int*,int*,int,int,boolean,int,int,int),int,int);
extern int nextelement(set*,int,int);
extern int orbjoin(int*,int*,int);
extern void permset(set*,set*,int,int*);
extern void putstring(FILE*,char*);
extern void refine(graph*,int*,int*,int,int*,int*,set*,int*,int,int);
extern void refine1(graph*,int*,int*,int,int*,int*,set*,int*,int,int);
extern void shortprune(set*,set*,int);
extern int targetcell(graph*,int*,int*,int,int,boolean,int,int,int);
extern int testcanlab(graph*,graph*,int*,int*,int,int);
extern void updatecan(graph*,graph*,int*,int,int,int);
extern void writeperm(FILE*,int*,boolean,int,int);
extern void nauty_freedyn(void);
extern void nauty_check(int,int,int,int);
extern void naugraph_check(int,int,int,int);
extern void nautil_check(int,int,int,int);
extern void nautil_freedyn(void);
extern void naugraph_freedyn(void);
extern void densenauty(graph*,int*,int*,int*,
			optionblk*,statsblk*,int,int,graph*);
extern void writegroupsize(FILE*,double,int);

#ifdef __cplusplus
}
#endif

/* ++++++ This file is automatically generated, don't edit it by hand! ++++++ */


#endif  /* _NAUTY_H_ */
