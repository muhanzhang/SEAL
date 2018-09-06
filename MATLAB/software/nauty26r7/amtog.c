/* amtog.c  version 2.0; B D McKay, Jun 2015. */

#define USAGE "amtog [-n#sgzhq] [-o#] [infile [outfile]]"

#define HELPTEXT \
" Read graphs in matrix format.\n\
\n\
    -n#   Set the initial graph order to # (no default).\n\
          This can be overridden in the input.\n\
    -g    Write the output in graph6 format (default).\n\
    -s    Write the output in sparse6 format.\n\
    -z    Write the output in digraph6 format.\n\
    -h    Write a header (according to -g or -s).\n\
    -w    Don't warn about loops (which are suppressed for -g)\n\
    -q    Suppress auxiliary information.\n\
    -o#   Treat digit # as 1 and other digits as 0.\n\
\n\
   Input consists of a sequence of commands restricted to:\n\
\n\
    n=#   set number of vertices (no default)\n\
          The = is optional.\n\
    m     Matrix to follow\n\
          An 'm' is also assumed if a digit is encountered.\n\
    M     Complement of matrix to follow (as m)\n\
    t     Upper triangle of matrix to follow, row by row\n\
           excluding the diagonal.\n\
    T     Complement of upper trangle to follow (as t)\n\
    q     exit (optional)\n"

/*************************************************************************/

#include "gtools.h"  /* which includes nauty.h and stdio.h */
#include <ctype.h>

/**************************************************************************/
/**************************************************************************/

int
main(int argc, char *argv[])
{
    int m,n,outdigit;
    int argnum,i,j,outcode,val;
    char *arg,sw,ochar;
    boolean badargs;
    boolean nswitch,sswitch,gswitch,hswitch,qswitch;
    boolean warn,loop,unsymm,compl,triangle;
    boolean zswitch,oswitch,nowarn;
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    nauty_counter nin;
    char s[10];
#if MAXN
    graph g[MAXN*MAXM];
#else
    DYNALLSTAT(graph,g,g_sz);
#endif

    HELP; PUTVERSION;

    sswitch = gswitch = zswitch = oswitch = FALSE;
    nowarn = qswitch = nswitch = hswitch = FALSE;
    infilename = outfilename = NULL;
    n = -1;

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
                     SWBOOLEAN('s',sswitch)
                else SWBOOLEAN('g',gswitch)
                else SWBOOLEAN('z',zswitch)
                else SWBOOLEAN('h',hswitch)
                else SWBOOLEAN('w',nowarn)
                else SWBOOLEAN('q',qswitch)
                else SWINT('n',nswitch,n,">E amtog -n")
                else SWINT('o',oswitch,outdigit,">E amtog -o")
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

    if ((sswitch!=0) + (gswitch!=0) + (zswitch!=0) > 1)
        gt_abort(">E amtog: -sgz are incompatible\n");

    if (oswitch && (outdigit < 0 || outdigit > 9))
        gt_abort(">E amtog: only digits are allowed for -o\n");

    if (badargs || argnum > 2)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (oswitch) ochar = outdigit + '0';
    else         ochar = '1';

    if (!infilename || infilename[0] == '-')
    {
        infilename = "stdin";
        infile = stdin;
    }
    else if ((infile = fopen(infilename,"r")) == NULL)
    {
        fprintf(stderr,"Can't open input file %s\n",infilename);
        gt_abort(NULL);
    }

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

    if (sswitch)      outcode = SPARSE6;
    else if (zswitch) outcode = DIGRAPH6;
    else              outcode = GRAPH6;

    if (hswitch)
    {
        if (outcode == SPARSE6)       writeline(outfile,SPARSE6_HEADER);
        else if (outcode == DIGRAPH6) writeline(outfile,DIGRAPH6_HEADER);
        else                          writeline(outfile,GRAPH6_HEADER);
    }

#if MAXN
    if (nswitch && n > MAXN)
    {
        gt_abort(">E amtog: value of -n too large\n");
        exit(2);
    }
#else
    if (nswitch)
    {
        m = (n + WORDSIZE - 1) / WORDSIZE;
        DYNALLOC2(graph,g,g_sz,n,m,"amtog");
    }
#endif
    

     /* perform scanning required */

    warn = FALSE;
    nin = 0;
    while (fscanf(infile,"%1s",s) == 1)
    {
        if (s[0] == 'n')
        {
            if (fscanf(infile,"=%d",&n) != 1)
            {
                gt_abort(">E amtog: invalid n=# command\n");
                exit(2);
            }
            m = (n + WORDSIZE - 1) / WORDSIZE;
#if MAXN
            if (n < 1 || n > MAXN || m > MAXM)
                gt_abort(">E amtog: n<0 or n,m too big\n");
#else
            DYNALLOC2(graph,g,g_sz,n,m,"amtog");
#endif
        } 
        else if (s[0] == 'm' || s[0] == 'M' || s[0] == 't' ||
                 s[0] == 'T' || s[0] == '0' || s[0] == '1' ||
                 (oswitch && isdigit(s[0])))
        {
            if (n < 0)
            {
                fprintf(stderr,
                    ">E amtog: matrix found before n is defined\n");
                exit(2);
            }
            if (isdigit(s[0])) ungetc(s[0],infile);
            m = (n + WORDSIZE - 1) / WORDSIZE;
    
            EMPTYSET(g,m*(size_t)n);

            loop = unsymm = FALSE;
            triangle = (s[0] == 't') || (s[0] == 'T');
            compl = (s[0] == 'M') || (s[0] == 'T');

            ++nin;
            for (i = 0; i < n; ++i)
            for (j = (triangle ? i+1 : 0); j < n; ++j)
            {
                if (fscanf(infile,"%1s",s) != 1)
                {
                    fprintf(stderr,">E amtog: incomplete matrix\n");
                    ABORT(">E amtog");
                }
                if (s[0] == '0' || s[0] == '1'
                      || (oswitch && isdigit(s[0])))
                {
                    val = ((i != j) & compl) ^ (s[0] == ochar);
                    if (val == 1)
                    {
                        if (triangle)
                        {
                            ADDELEMENT(GRAPHROW(g,i,m),j);
                            ADDELEMENT(GRAPHROW(g,j,m),i);
                        }
                        else
                        {
                            if (j < i && !ISELEMENT(GRAPHROW(g,j,m),i))
                                unsymm = TRUE;
                            ADDELEMENT(GRAPHROW(g,i,m),j);
                        }
                        if (i == j) loop = TRUE;
                    }
                    else if (j < i && ISELEMENT(GRAPHROW(g,j,m),i))
                        unsymm = TRUE;
                }
                else
                {
                    fprintf(stderr,
                      ">E amtog: illegal character in matrix: \"%c\"\n",
                      s[0]);
                    gt_abort(NULL);
                }
            }

            if (unsymm && outcode != DIGRAPH6) fprintf(stderr,
                   ">W amtog: warning, graph "
                          COUNTER_FMT " is unsymmetric\n",nin);
    
            if (outcode == DIGRAPH6)     writed6(outfile,g,m,n);
            else if (outcode == SPARSE6) writes6(outfile,g,m,n);
            else                         writeg6(outfile,g,m,n);
	    if (loop && outcode == GRAPH6) ++warn;
        }
        else if (s[0] == 'q')
        {
            exit(0);
        }
        else
        {
            fprintf(stderr,">E amtog: invalid command \"%c\"\n",s[0]);
            gt_abort(NULL);
        }
    }

    if (warn) fprintf(stderr,">Z complg: loops were lost\n");

    if (!qswitch)
        fprintf(stderr,">Z  " COUNTER_FMT " graphs converted from %s to %s.\n",
                       nin,infilename,outfilename);

    exit(0);
}
