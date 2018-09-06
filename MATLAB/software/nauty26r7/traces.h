/******************************************************************************
 *                                                                            *
 * This is the header file for traces() version 2.1, which is included into   *
 *   nauty() version 2.6.                                                     *
 *                                                                            *
 *   nauty is Copyright (1984-2016) Brendan McKay.  All rights reserved.      *
 *   Traces is Copyright Adolfo Piperno, 2008-2016.  All rights reserved.     *
 *   See the file COPYRIGHT for the details of the software license.          *
 *                                                                            *
 *   CHANGE HISTORY                                                           *
 *       28-Dec-12 : final changes for version 2.0                            *
 *       20-Jan-13 : add code for ^C catching in Traces                       *
 *       29-Mar-13 : bug correction in automorphism mode                      *
 *       02-Apr-13 : add preprocessing                                        *
 *       21-May-13 : bug correction (coloured lists)                          *
 *       29-Jun-13 : bug correction (coloured lists and cycles)               *
 *       07-Dec-13 : bug correction in automorphism mode (wrong group size    *
 *                   due to randomness in Schreier-Sims orbit computation)    *
 *                   bug correction (discrete initial partition)              *
 *       15-Feb-14 : CPUDEFS removed (already declared in gtools.h)           *
 *       01-Sep-15 : add weighted edges (not active)                          *
 *       28-Jan-16 : version ready for nauty and Traces v.2.6 distribution    *
 *       12-Jul-16 : bug correction (reaching degree 2 vertices)              *
*****************************************************************************/

#include "gtools.h"
#include "schreier.h" 

typedef struct TracesOptions {
	boolean getcanon;
	boolean writeautoms;
	boolean cartesian;
	boolean digraph;
	boolean defaultptn;
	int linelength;
	FILE* outfile;
	int strategy;    /* Only the value 0 is supported in this version. */
	int verbosity;
	permnode **generators;
    void (*userautomproc)(int,int*,int);
    int  (*usercanonproc)(graph*,int*,graph*,int,int,int,int);
    boolean weighted;
} TracesOptions;

#define DEFAULTOPTIONS_TRACES(opts) TracesOptions opts \
= { FALSE, FALSE, FALSE, FALSE, TRUE, 0, NULL, 0, 0, NULL, NULL, NULL, FALSE }

typedef struct TracesStats {
	double grpsize1;
	int grpsize2;
	int numgenerators;
	int numorbits;
	int treedepth;
	int canupdates;
	int errstatus;
	unsigned long numnodes;
	unsigned long interrupted;
	unsigned long peaknodes;
} TracesStats;

extern void Traces(sparsegraph*,int*,int*,int*,TracesOptions*,
				   TracesStats*,sparsegraph*);									
extern void refine_tr(sparsegraph*,int*,int*,int*,int*,TracesOptions*);		
extern void traces_freedyn(void);
