/* converseg.c  version 2.0; B D McKay, Jun 2015. */

#define USAGE "converseg [-q] [infile [outfile]]"

#define HELPTEXT \
" Take the converse digraphs of a file of graphs.\n\
\n\
    The output file has a header if and only if the input file does.\n\
    Undirected graphs are passed through without change, while\n\
    directed graphs are written in digraph6 format.\n\
\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 

/**************************************************************************/

static void
conv(graph *g, int m, int n)
/* Replace g by its converse */
{
    int i,j;
    graph *gi,*gj;

    for (i = 0, gi = g; i < n; ++i, gi += m)
        for (j = i+1, gj = gi+m; j < n; ++j, gj += m)
            if ((ISELEMENT(gi,j)!=0) + (ISELEMENT(gj,i)!=0) == 1)
            {
                FLIPELEMENT(gi,j);
                FLIPELEMENT(gj,i);
            }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,quiet;
    boolean digraph;
    int j,m,n,argnum;
    int codetype,outcode;
    graph *g;
    // size_t ii,ned,nedc,nn,loops,loopsc,gwords;
    nauty_counter nin;
    char *arg,sw;
    double t;

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    badargs = FALSE;
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
        fprintf(stderr,">A converseg");
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

    gtools_check(WORDSIZE,1,1,NAUTYVERSIONID);

    nin = 0;
    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
        ++nin;

        if (!digraph)
	    writelast(outfile);
        else
        {
	    conv(g,m,n);
	    writed6(outfile,g,m,n);
	}
        FREES(g);
    }
    t = CPUTIME - t;

    if (!quiet)
        fprintf(stderr,">Z  " COUNTER_FMT 
                " graphs converted from %s to %s in %3.2f sec.\n",
                nin,infilename,outfilename,t);

    exit(0);
}
