/* complg.c  version 2.0; B D McKay, Jun 2015. */

#define USAGE "complg [-lq] [-r|-R] [infile [outfile]]"

#define HELPTEXT \
" Take the complements of a file of graphs.\n\
\n\
    The output file has a header if and only if the input file does.\n\
    The output format is defined by the header or first graph.\n\
\n\
    -r  Only complement if the complement has fewer directed edges.\n\
    -R  Only complement if the complement has fewer directed edges\n\
        or has the same number of directed edges and is canonically\n\
        less than the original.\n\
    -L  Complement the loops too. By default, preserve them.\n\
    -l  Canonically label outputs.\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 

/**************************************************************************/

static void
compl(graph *g, int m, int n, graph *h, boolean comploops)
/* h := complement of g */
{
    int i,j;
    setword *gi,*hi;
#if MAXN
    set all[MAXM];
#else
    DYNALLSTAT(set,all,all_sz);
    DYNALLOC1(set,all,all_sz,m,"complg");
#endif

    EMPTYSET(all,m);
    for (i = 0; i < n; ++i) ADDELEMENT(all,i);

    gi = (setword*) g;
    hi = (setword*) h;

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j) hi[j] = gi[j] ^ all[j];
        if (!comploops) FLIPELEMENT(hi,i);
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
    boolean dolabel,badargs,restricted,Restricted,quiet;
    boolean digraph,Lswitch;
    int j,m,n,argnum;
    int codetype,outcode;
    graph *g;
    size_t ii,ned,nedc,nn,loops,loopsc,gwords;
    nauty_counter nin;
    char *arg,sw;
    static graph *gq;
    double t;
#if MAXN
    graph gc[MAXN*MAXM],h[MAXN*MAXM],hc[MAXN*MAXM];
#else
    DYNALLSTAT(graph,gc,gc_sz);
    DYNALLSTAT(graph,h,h_sz);
    DYNALLSTAT(graph,hc,hc_sz);
#endif

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    dolabel = badargs = FALSE;
    restricted = Restricted = quiet = Lswitch = FALSE;

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
                     SWBOOLEAN('r',restricted)
                else SWBOOLEAN('R',Restricted)
                else SWBOOLEAN('l',dolabel)
                else SWBOOLEAN('L',Lswitch)
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
        fprintf(stderr,">A complg");
        if (restricted || Restricted || dolabel || Lswitch)
            fprintf(stderr," -");
        if (restricted && !Restricted) fprintf(stderr,"r");
        if (Restricted) fprintf(stderr,"R");
        if (dolabel) fprintf(stderr,"l");
        if (Lswitch) fprintf(stderr,"L");
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

    if      (codetype&SPARSE6)  outcode = SPARSE6;
    else if (codetype&DIGRAPH6) outcode = DIGRAPH6;
    else                        outcode = GRAPH6;

    if (codetype&HAS_HEADER)
    {
        if (outcode == SPARSE6)       writeline(outfile,SPARSE6_HEADER);
        else if (outcode == DIGRAPH6) writeline(outfile,DIGRAPH6_HEADER);
        else          	              writeline(outfile,GRAPH6_HEADER);
    }

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    nin = 0;
    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
        ++nin;
#if !MAXN
        DYNALLOC2(graph,gc,gc_sz,n,m,"complg");
#endif

        if (restricted || Restricted)
        {
            ned = loops = 0;
            gwords = m * (size_t)n;
            nn = n * (size_t)(n-1);
            for (ii = 0; ii < gwords; ++ii) ned += POPCOUNT(g[ii]);
            for (ii = 0; ii < n; ++ii) if (ISELEMENT(g+m*ii,ii)) ++loops;
            if (Lswitch)
            {
		loopsc = n - loops;
                nedc = nn + n - (ned-loops)/2;
            }
            else
            {
		loopsc = loops;
                nedc = nn - (ned-3*loops)/2;
            }

            if (ned > nedc || (ned == nedc && !Restricted))
            {
                compl(g,m,n,gc,Lswitch);
                if (dolabel)
                {
#if !MAXN
                    DYNALLOC2(graph,hc,hc_sz,n,m,"complg");
#endif
                    fcanonise(gc,m,n,hc,NULL,digraph||loopsc>0);
            	    gq = hc;
                }
                else
            	    gq = gc;
            }
            else if (ned < nedc)
            {
                if (dolabel)
                {
#if !MAXN
                    DYNALLOC2(graph,h,h_sz,n,m,"complg");
#endif
                    fcanonise(g,m,n,h,NULL,digraph||loops>0);
            	    gq = h;
                }
                else
            	    gq = g;
            }
            else 
            {
                compl(g,m,n,gc,Lswitch);
#if !MAXN
                DYNALLOC2(graph,h,h_sz,n,m,"complg");
                DYNALLOC2(graph,hc,hc_sz,n,m,"complg");
#endif
                fcanonise(g,m,n,h,NULL,digraph||loops>0);
                fcanonise(gc,m,n,hc,NULL,digraph||loopsc>0);
                for (ii = 0; ii < gwords; ++ii)
            	if (h[ii] != hc[ii]) break;
                if (ii == gwords || hc[ii] < h[ii])
                {
            	    if (dolabel) gq = hc; else gq = gc;
                }
                else
                {
            	    if (dolabel) gq = h; else gq = g;
                }
            }
        }
        else   /* Not restricted */
        {		
            compl(g,m,n,gc,Lswitch);
            if (dolabel)
            {
#if !MAXN
                DYNALLOC2(graph,hc,hc_sz,n,m,"complg");
#endif
                fcanonise(gc,m,n,hc,NULL,digraph||loopsc>0);
                gq = hc;
            }
            else
                gq = gc;
        }

        if (outcode == SPARSE6)       writes6(outfile,gq,m,n);
        else if (outcode == DIGRAPH6) writed6(outfile,gq,m,n);
        else                          writeg6(outfile,gq,m,n);
        FREES(g);
    }
    t = CPUTIME - t;

    if (!quiet)
        fprintf(stderr,">Z  " COUNTER_FMT 
                " graphs converted from %s to %s in %3.2f sec.\n",
                nin,infilename,outfilename,t);

    exit(0);
}
