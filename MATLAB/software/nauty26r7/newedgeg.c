/* newedgeg.c  version 1.1; B D McKay, Dec 2005. */

#define USAGE "newedgeg [-lq] [infile [outfile]]"

#define HELPTEXT \
" For each pair of non-adjacent edges, output the graph obtained\n\
  by subdividing the edges and joining the new vertices.\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -l  Canonically label outputs\n\
    -q  Suppress auxiliary information\n"

/*************************************************************************/

#include "gtools.h" 
#include "gutils.h"

static FILE *outfile;
static nauty_counter nout;
static int outcode;

/*************************************************************************/

static void
newedge(graph *g1, int m1, int n1,
        int v1, int v2, int w1, int w2,
        graph *g2, int m2) 
/* Make g2 by subdividing edges v1-v2 and w1-w2 in g1
   and adding an edge between them.  Must have m2 >= m1.  */
{
	int i,j;
	setword *s1,*s2;

	s1 = g1;
	s2 = g2;
	for (i = 0; i < n1; ++i)
	{
	    for (j = 0; j < m1; ++j) *(s2++) = *(s1++);
	    for (; j < m2; ++j) *(s2++) = 0;
	}

	s2 = GRAPHROW(g2,v1,m2);
	DELELEMENT(s2,v2);
	ADDELEMENT(s2,n1);

	s2 = GRAPHROW(g2,v2,m2);
        DELELEMENT(s2,v1);
        ADDELEMENT(s2,n1);

	s2 = GRAPHROW(g2,w1,m2);
        DELELEMENT(s2,w2);
        ADDELEMENT(s2,n1+1);

	s2 = GRAPHROW(g2,w2,m2);
        DELELEMENT(s2,w1);
        ADDELEMENT(s2,n1+1);

	s2 = GRAPHROW(g2,n1,m2);
	EMPTYSET(s2,m2);
	ADDELEMENT(s2,v1);
	ADDELEMENT(s2,v2);
	ADDELEMENT(s2,n1+1);

	s2 = GRAPHROW(g2,n1+1,m2);
	EMPTYSET(s2,m2);
	ADDELEMENT(s2,w1);
	ADDELEMENT(s2,w2);
	ADDELEMENT(s2,n1);
}

/*************************************************************************/

static void
na_newedge(graph *g1, int m1, int n1, boolean dolabel)
/* Make all graphs by non-adjacent edge addition. */
{
	int n2,m2;
	int v1,v2,w1,w2;
	set *sv1,*sw1;
	graph *gq;
#if MAXN
	graph h[MAXN*MAXM];
	grapg g2[MAXN*MAXM];
#else
	DYNALLSTAT(graph,h,h_sz);
	DYNALLSTAT(graph,g2,g2_sz);
#endif

	n2 = n1 + 2;
	m2 = (n2 + WORDSIZE - 1) / WORDSIZE;
	if (m2 < m1) m2 = m1;

#if !MAXN
	DYNALLOC2(graph,g2,g2_sz,m2,n2,"newedgeg");
	if (dolabel) DYNALLOC2(graph,h,h_sz,m2,n2,"newedgeg");
#endif

	for (v1 = 0, sv1 = g1; v1 < n1-3; ++v1, sv1 += m1)
	for (w1 = v1+1, sw1 = sv1 + m1; w1 < n1-1; ++w1, sw1 += m1)
	{
	    for (v2 = v1; (v2 = nextelement(sv1,m1,v2)) >= 0; )
	    for (w2 = w1; (w2 = nextelement(sw1,m1,w2)) >= 0; )
	    {
		if (v2 == w1 || v2 == w2) continue;

		newedge(g1,m1,n1,v1,v2,w1,w2,g2,m2);

                gq = g2;

                if (dolabel)
                {
                    fcanonise(g2,m2,n2,h,NULL,FALSE);  /* FIXME (loops) */
                    gq = h;
                }
                if (outcode == SPARSE6) writes6(outfile,gq,m2,n2);
                else                    writeg6(outfile,gq,m2,n2);
                ++nout;
	    }
	}
}

/*************************************************************************/

int
main(int argc, char *argv[])
{
        char *infilename,*outfilename;
        FILE *infile;
        boolean badargs,dolabel,quiet;
	int j,m,n,argnum;
	int codetype;
	graph *g;
	nauty_counter nin;
        char *arg,sw;
	double t;

	HELP; PUTVERSION;

        infilename = outfilename = NULL;
	dolabel = quiet = FALSE;

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
		         SWBOOLEAN('l',dolabel)
		    else SWBOOLEAN('q',quiet)
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
	    fprintf(stderr,">A newedgeg");
	    if (dolabel) fprintf(stderr," -l");
	    if (argnum > 0) fprintf(stderr," %s",infilename);
	    if (argnum > 1) fprintf(stderr," %s",outfilename);
	    fprintf(stderr,"\n");
	    fflush(stderr);
	}

	if (infilename && infilename[0] == '-') infilename = NULL;
	infile = opengraphfile(infilename,&codetype,FALSE,1);
	if (!infile) exit(1);
	if (!infilename) infilename = "stdin";

        NODIGRAPHSYET(codetype);

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

	if (dolabel) nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

	nin = nout = 0;
	t = CPUTIME;
	while (TRUE)
	{
	    if ((g = readg(infile,NULL,0,&m,&n)) == NULL) break;
#if MAXN
	    if (n > MAXN-2)
	    {
		fprintf(stderr,">E newedgeg: input too large\n");
		exit(1);
	    }
#endif

	    ++nin;

	    na_newedge(g,m,n,dolabel);

	    FREES(g);
	}
	t = CPUTIME - t;

        if (!quiet)
            fprintf(stderr,
                ">Z  " COUNTER_FMT " graphs read from %s, "
                       COUNTER_FMT " written to %s; %3.2f sec.\n",
                    nin,infilename,nout,outfilename,t);

	exit(0);
}
