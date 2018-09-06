/* NRswitchg.c  nauty version 2.4; B D McKay, Dec 2005 */

#define USAGE "NRswitchg [-lq] [infile [outfile]]"

#define HELPTEXT \
" For each v, complement the edges from N(v) to V(G)-N(v)-v.\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -l  Canonically label outputs.\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 

/**************************************************************************/

void
NRswitch(graph *g, int m, int n, int v, graph *h)
/* h := g with N-R edge set complemented */
{
	register int i,j;
	register setword *gv,*gi,*hi;
#if MAXN
	set Nset[MAXM];
	set Rset[MAXM];
#else
	DYNALLSTAT(set,Nset,Nset_sz);
	DYNALLSTAT(set,Rset,Rset_sz);
	DYNALLOC1(set,Nset,Nset_sz,m,"NRswitchg");
	DYNALLOC1(set,Rset,Rset_sz,m,"NRswitchg");
#endif

	EMPTYSET(Rset,m);
	for (i = 0; i < n; ++i) ADDELEMENT(Rset,i);

	gv = GRAPHROW(g,v,m);

	for (j = 0; j < m; ++j)
	{
	    Nset[j] = gv[j];
	    Rset[j] ^= gv[j];
	}

	FLIPELEMENT(Rset,v);

	gi = (setword*) g;
	hi = (setword*) h;

	for (i = 0; i < n; ++i)
	{
	    if (i == v)
	        for (j = 0; j < m; ++j) hi[j] = gi[j];
	    else if (ISELEMENT(Nset,i))
		for (j = 0; j < m; ++j) hi[j] = gi[j] ^ Rset[j];
	    else
		for (j = 0; j < m; ++j) hi[j] = gi[j] ^ Nset[j];

	    gi += m;
	    hi += m;
	}
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
        char *infilename,*outfilename;
        FILE *infile,*outfile;
        boolean badargs,dolabel,quiet;
	int j,m,n,v,argnum;
	int codetype,outcode;
	graph *g;
	nauty_counter nin,nout;
        char *arg,sw;
	static graph *gq;
	double t;
#if MAXN
	graph gc[MAXN*MAXM],h[MAXN*MAXM];
#else
	DYNALLSTAT(graph,gc,gc_sz);
	DYNALLSTAT(graph,h,h_sz);
#endif

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
	    fprintf(stderr,">A NRswitchg");
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
	    ++nin;
#if !MAXN
	    DYNALLOC2(graph,gc,gc_sz,n,m,"NRswitchg");
#endif
	    for (v = 0; v < n; ++v)
	    {		
	    	NRswitch(g,m,n,v,gc);
		gq = gc;
	    
	        if (dolabel)
	        {
#if !MAXN
		    DYNALLOC2(graph,h,h_sz,n,m,"compl");
#endif
	 	    fcanonise(gq,m,n,h,NULL,FALSE);
		    gq = h;
	        }
	        if (outcode == SPARSE6) writes6(outfile,gq,m,n);
	        else                    writeg6(outfile,gq,m,n);
		++nout;
	    }
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
