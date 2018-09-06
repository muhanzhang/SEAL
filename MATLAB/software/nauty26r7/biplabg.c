/* biplabg.c  version 1.1; B D McKay, Nov 10, 2009. */

#define USAGE "biplabg [-q] [infile [outfile]]"

#define HELPTEXT \
" Label bipartite graphs so that the colour classes are contiguous.\n\
  The first vertex of each component is assigned the first colour.\n\
  Vertices in each colour class have the same relative order as before.\n\
  Non-bipartite graphs are rejected.\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 
#include "gutils.h"

/**************************************************************************/

static boolean
biplabel(graph *g, int m, int n, graph *h)
/* h := bipartite labelling of g; else return FALSE */
{
	int i,j;
#if MAXN
	int colour[MAXN];
	int lab[MAXN];
#else
	DYNALLSTAT(int,colour,colour_sz);
	DYNALLSTAT(int,lab,lab_sz);

	DYNALLOC1(int,colour,colour_sz,n,"biplabg");
	DYNALLOC1(int,lab,lab_sz,n,"biplabg");
#endif

	if (!twocolouring(g,colour,m,n)) return FALSE;

	j = 0;
	for (i = 0; i < n; ++i) if (colour[i] == 0) lab[j++] = i;
	for (i = 0; i < n; ++i) if (colour[i] == 1) lab[j++] = i;

	updatecan(g,h,lab,0,m,n);

	return TRUE;
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
        char *infilename,*outfilename;
        FILE *infile,*outfile;
        boolean badargs,quiet;
	int j,m,n,argnum;
	int codetype,outcode;
	graph *g;
	nauty_counter nin,nout;
        char *arg,sw;
	double t;
#if MAXN
	graph h[MAXN*MAXM];
#else
	DYNALLSTAT(graph,h,h_sz);
#endif

	HELP; PUTVERSION;

        infilename = outfilename = NULL;
	quiet = FALSE;

	argnum = 0;
	badargs = FALSE;
	for (j = 1; !badargs && j < argc; ++j)
	{
	    arg = argv[j];
	    if (arg[0] == '-' && arg[1] != '\0')
	    {
		++arg;
		while (*arg != '\0')
		{
		    sw = *arg++;
		         SWBOOLEAN('q',quiet)
		    else badargs = TRUE;
		}
	    }
	    else
	    {
		++argnum;
		if      (argnum == 1) infilename = arg;
	        else if (argnum == 2) outfilename = arg;
		else                  badargs = TRUE;
	    }
	}

	if (badargs)
	{
	    fprintf(stderr,">E Usage: %s\n",USAGE);
	    GETHELP;
	    exit(1);
	}

	if (!quiet)
	{
	    fprintf(stderr,">A biplabg");
	    if (argnum > 0) fprintf(stderr," %s",infilename);
	    if (argnum > 1) fprintf(stderr," %s",outfilename);
	    fprintf(stderr,"\n");
	    fflush(stderr);
	}

	NODIGRAPHSYET(codetype);

	if (infilename && infilename[0] == '-') infilename = NULL;
	infile = opengraphfile(infilename,&codetype,FALSE,1);
	if (!infile) exit(1);
	if (!infilename) infilename = "stdin";

	if (!outfilename || outfilename[0] == '-')
	{
	    outfilename = "stdout";
	    outfile = stdout;
	}
	else if ((outfile = fopen(outfilename,"w")) == NULL)
	{
	    fprintf(stderr,"Can't open output file %s\n",outfilename);
	    gt_abort(NULL);
	}

	if (codetype&SPARSE6) outcode = SPARSE6;
	else                  outcode = GRAPH6;

	if (codetype&HAS_HEADER)
	{
	    if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
	    else    		    writeline(outfile,GRAPH6_HEADER);
	}

	nautil_check(WORDSIZE,1,1,NAUTYVERSIONID);

	nin = nout = 0;
	t = CPUTIME;
	while (TRUE)
	{
	    if ((g = readg(infile,NULL,0,&m,&n)) == NULL) break;
	    ++nin;

#if !MAXN
	    DYNALLOC2(graph,h,h_sz,n,m,"biplabg");
#endif
	    if (biplabel(g,m,n,h))
	    {
		++nout;
	        if (outcode == SPARSE6) writes6(outfile,h,m,n);
	        else                    writeg6(outfile,h,m,n);
	    }
	    FREES(g);
	}
	t = CPUTIME - t;

        if (!quiet)
            fprintf(stderr,
                ">Z  " COUNTER_FMT " graphs read from %s; "
                       COUNTER_FMT " written to %s; %3.2f sec.\n",
                    nin,infilename,nout,outfilename,t);

	exit(0);
}
