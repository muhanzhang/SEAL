/* delptg.c  version 2.0; B D McKay, Jul 2015. */

#define USAGE "delptg [-lq] [-d#|d#:#] [-n#] [infile [outfile]]"

#define HELPTEXT \
" Delete some vertices from a file of graphs.\n\
\n\
    The output file has a header if and only if the input file does.\n\
    No isomorph reduction is done.\n\
\n\
    -l  Canonically label outputs\n\
    -d# -d#:# Only remove vertices with original degree in the given range\n\
        For digraphs, the out-degree is used.\n\
    -n# The number of vertices to delete (default 1).\n\
        Inputs that are too small are ignored; no empty graphs are output.\n\
    -q  Suppress auxiliary information\n"

/*************************************************************************/

#include "gtools.h" 

static FILE *outfile;
static nauty_counter nout;
static int outcode;
static boolean digraph,dolabel;

/**************************************************************************/

static void
writeone(graph *g, int m, int n, int *del, int ndel)
/* Delete the stated vertices and write */
{
    int i,j,k,nx,mx;
    graph *gi,*gxi,*gq;
#if MAXN
    graph gx[MAXN*MAXM];
    graph hx[MAXN*MAXM];
    int lab[MAXN];
#else
    DYNALLSTAT(graph,gx,gx_sz);
    DYNALLSTAT(graph,hx,hx_sz);
    DYNALLSTAT(int,lab,lab_sz);
#endif

    nx = n - ndel;    /* Always positive */
    mx = SETWORDSNEEDED(nx);

#if !MAXN
    DYNALLOC2(graph,gx,gx_sz,mx,nx,"delptg");
    if (dolabel) DYNALLOC2(graph,hx,hx_sz,mx,nx,"delptg");
    DYNALLOC1(int,lab,lab_sz,nx,"delptg");
#endif

    j = 0;
    for (i = 0; i < del[0]; ++i) lab[j++] = i;
    for (k = 1; k < ndel; ++k)
	for (i = del[k-1]+1; i < del[k]; ++i) lab[j++] = i;
    for (i = del[ndel-1]+1; i < n; ++i) lab[j++] = i;

    EMPTYSET(gx,nx*(size_t)mx);

    for (i = 0, gxi = (set*)gx; i < nx; ++i, gxi += mx)
    {
	gi = GRAPHROW(g,lab[i],m);
	for (j = 0; j < nx; ++j)
	{
	    k = lab[j];
            if (ISELEMENT(gi,k)) ADDELEMENT(gxi,j);
	}
    }

    if (dolabel)
    {
     	fcanonise(gx,mx,nx,hx,NULL,digraph);
    	gq = hx;
    }
    else
	gq = gx;

    if (outcode == DIGRAPH6 || digraph) writed6(outfile,gq,mx,nx);
    else if (outcode == SPARSE6)        writes6(outfile,gq,mx,nx);
    else                                writeg6(outfile,gq,mx,nx);
    ++nout;
}

/**************************************************************************/

static void
search(int level, int ndel, int *del, graph *g, int m, int n, boolean *degok)
{
    int i;

    if (level == ndel)
    {
	writeone(g,m,n,del,ndel);
	return;
    }

    for (i = (level == 0 ? 0 : del[level-1]+1); i <= n-ndel+level; ++i)
	if (degok[i])
	{
	    del[level] = i;
	    search(level+1,ndel,del,g,m,n,degok);
        }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile;
    boolean badargs,quiet,dswitch,nswitch;
    int i,j,m,n,v,argnum;
    int ndel;
    int codetype;
    graph *g;
    nauty_counter nin;
    char *arg,sw;
    setword *gv;
    long mindeg,maxdeg;
    int degv;
    double t;
#if MAXN
    boolean degok[MAXN];
    boolean del[MAXN];
#else
    DYNALLSTAT(boolean,degok,degok_sz);
    DYNALLSTAT(boolean,del,del_sz);
#endif

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    badargs = FALSE;
    dswitch = nswitch = quiet = FALSE;

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
		else SWINT('n',nswitch,ndel,">E delptg -n")
                else SWRANGE('d',":-",dswitch,mindeg,maxdeg,">E delptg -d")
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

    if (nswitch && ndel < 1) gt_abort(">E delptg: bad argument for -n\n");

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!quiet)
    {
        fprintf(stderr,">A delptg");
        if (dolabel) fprintf(stderr," -l");
        if (dswitch) fprintf(stderr," -d%ld:%ld",mindeg,maxdeg);
        if (nswitch) fprintf(stderr," -n%d",ndel);
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    if (!dswitch)
    {
	mindeg = 0;
	maxdeg = NAUTY_INFINITY;
    }

    if (!nswitch) ndel = 1;

    if (dolabel) nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    nin = nout = 0;
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

    if (codetype&SPARSE6)       outcode = SPARSE6;
    else if (codetype&DIGRAPH6) outcode = DIGRAPH6;
    else                        outcode = GRAPH6;

    if (codetype&HAS_HEADER)
    {
        if (outcode == SPARSE6)       writeline(outfile,SPARSE6_HEADER);
        else if (outcode == DIGRAPH6) writeline(outfile,DIGRAPH6_HEADER);
        else    	              writeline(outfile,GRAPH6_HEADER);
    }

    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
        ++nin;

        if (n <= ndel)
        {
	    FREES(g);
	    continue;
        }

#if !MAXN
        DYNALLOC1(boolean,degok,degok_sz,n,"delptg");
        DYNALLOC1(int,del,del_sz,n,"delptg");
#endif

        for (v = 0, gv = g; v < n; ++v, gv += m)
        {
    	    degv = 0;
    	    for (i = 0; i < m; ++i)
    	        degv += POPCOUNT(gv[i]);
    	    degok[v] = (degv >= mindeg) && (degv <= maxdeg);
        }

	search(0,ndel,del,g,m,n,degok);

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
