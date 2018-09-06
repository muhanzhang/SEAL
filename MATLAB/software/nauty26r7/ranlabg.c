/* ranlabg.c  version 1.2; B D McKay, Jun 20, 2015. */

#define USAGE "ranlabg [-q] [-f#] [-S#] [infile [outfile]]"

#define HELPTEXT \
" Randomly relabel graphs.\n\
\n\
    The output file has a header if and only if the input file does.\n\
    Each graph is written in the same format as it is read.\n\
\n\
    -f# Don't relabel the first # vertices.\n\
    -m# Output # randomly labelled copies of each input (default 1).\n\
    -S# Set random number seed (taken from clock otherwise).\n\
\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 

/**************************************************************************/

static void
ranrelabel(graph *g, int fixed, int m, int n, graph *h)
/* h := random labelling of g, fixing some initial vertices */
{
    int i,j,w,imin;
#if MAXN
    int perm[MAXN];
#else
    DYNALLSTAT(int,perm,perm_sz);
    DYNALLOC1(int,perm,perm_sz,n,"ranlabg");
#endif

    for (i = 0 ; i < n; ++i) perm[i] = i;

    if (fixed >= 0 && fixed < n) imin = fixed;
    else                         imin = 0;

    for (i = imin; i < n-1; ++i)
    {
        j = i + KRAN(n-i);
        w = perm[i];
        perm[i] = perm[j];
        perm[j] = w;
    }

    updatecan(g,h,perm,0,m,n);
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,quiet,Sswitch,fswitch,mswitch;
    int j,m,n,argnum,fixed;
    int codetype,outcode;
    graph *g;
    nauty_counter nin;
    int mult;
    boolean digraph;
    char *arg,sw;
    double t;
    long seed;
#if MAXN
    graph h[MAXN*MAXM];
#else
    DYNALLSTAT(graph,h,h_sz);
#endif

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    mswitch = quiet = fswitch = Sswitch = FALSE;

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
                else SWINT('f',fswitch,fixed,"ranlabg -f") 
                else SWINT('m',mswitch,mult,"ranlabg -m") 
                else SWLONG('S',Sswitch,seed,"ranlabg -S") 
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

    if (!Sswitch) INITSEED;
    if (!mswitch) mult = 1;
    if (!fswitch) fixed = 0;

    if (!quiet)
    {
        fprintf(stderr,">A ranlabg");
        fprintf(stderr," -S%ld",seed);
        if (fswitch) fprintf(stderr,"f%d",fixed);
        if (mswitch) fprintf(stderr,"m%d",mult);
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    ran_init(seed);

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
        else                          writeline(outfile,GRAPH6_HEADER);
    }

    nautil_check(WORDSIZE,1,1,NAUTYVERSIONID);

    nin = 0;
    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
        ++nin;

#if !MAXN
        DYNALLOC2(graph,h,h_sz,n,m,"ranlabg");
#endif
        for (j = 0; j < mult; ++j)
        {
            ranrelabel(g,fixed,m,n,h);
            if (readg_code == SPARSE6 || readg_code == INCSPARSE6)
                writes6(outfile,h,m,n);
            else if (readg_code == GRAPH6)
                writeg6(outfile,h,m,n);
            else
                writed6(outfile,h,m,n);
        }
        FREES(g);
    }
    t = CPUTIME - t;

    if (!quiet)
        fprintf(stderr,
            ">Z  " COUNTER_FMT
                 " graphs relabeled from %s to %s; %3.2f sec.\n",
                nin,infilename,outfilename,t);

    exit(0);
}
