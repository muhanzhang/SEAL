/* dretog.c  version 1.1; B D McKay, Jan 2013. */

#define USAGE "dretog [-n#o#sghq] [infile [outfile]]"

#define HELPTEXT \
" Read graphs in dreadnaut format.\n\
\n\
   -o#   Label vertices starting at # (default 0).  \n\
         This can be overridden in the input.\n\
   -n#   Set the initial graph order to # (no default).  \n\
         This can be overridden in the input.\n\
   -g    Use graph6 format (default for undirected graphs).\n\
   -z    Use digraph6 format (default for directed graphs).\n\
   -s    Use sparse6 format.\n\
   -h    Write a header (according to -g, -z or -s).\n\
   -q    Suppress auxiliary output.\n\
\n\
  Input consists of a sequence of dreadnaut commands restricted to:\n\
\n\
   n=#   set number of vertices (no default)\n\
         The = is optional.\n\
   $=#   set label of first vertex (default 0)\n\
         The = is optional.\n\
   d     indicate graph will be directed\n\
   $$    return origin to initial value (see -o#)\n\
   \"..\" and !..\\n   comments to ignore\n\
   g     specify graph to follow (as dreadnaut format)\n\
         Can be omitted if first character of graph is a digit or ';'.\n\
   q     exit (optional)\n"

/*************************************************************************/

#include "gtools.h"  /* which includes nauty.h and stdio.h */

/**************************************************************************/
/**************************************************************************/

int
main(int argc, char *argv[])
{
	int m,n,c;
	int argnum,j,outcode,initorg;
	char *arg,sw;
	boolean badargs,prompt,digraph;
	boolean zswitch,sswitch,gswitch,oswitch,nswitch,hswitch,qswitch;
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

	zswitch = sswitch = gswitch = oswitch = FALSE;
	qswitch = nswitch = hswitch = FALSE;
	infilename = outfilename = NULL;
	initorg = 0;
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
		    else SWBOOLEAN('q',qswitch)
		    else SWINT('o',oswitch,initorg,">E dretog -o")
		    else SWINT('n',nswitch,n,">E dretog -n")
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
            gt_abort(">E dretog: -s, -z and -g are incompatible\n");

	if (labelorg < 0) gt_abort(">E dretog: negative origin forbidden\n");

	if (badargs || argnum > 2)
	{
	    fprintf(stderr,">E Usage: %s\n",USAGE);
	    GETHELP;
	    exit(1);
	}

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
	    else    		          writeline(outfile,GRAPH6_HEADER);
	}

#if HAVE_ISATTY
        prompt = isatty(fileno(infile)) && isatty(fileno(outfile));
#else
	prompt = (infile == stdin);
#endif

     /* perform scanning required */

	labelorg = initorg;
	nin = 0;
        digraph = FALSE;

	while (fscanf(infile,"%1s",s) == 1)
	{
	    if (s[0] == 'n')
	    {
		if (fscanf(infile,"%1s",s) == 1 && s[0] != '=')
		    ungetc(s[0],infile);
		if (fscanf(infile,"%d",&n) != 1)
		{
		    fprintf(stderr,">E dretog: invalid n=# command\n");
		    gt_abort(NULL);
		}
		if (n <= 0)
		    gt_abort(">E dretog: n can't be <= 0\n");
	    } 
	    else if (s[0] == 'd')
		digraph = TRUE;
	    else if (s[0] == '"')
	    {
		while ((c = getc(infile)) != '"' && c != EOF) {}
	    }
	    else if (s[0] == '!')
	    {
		while ((c = getc(infile)) != '\n' && c != EOF) {}
	    }
	    else if (s[0] == '$')
	    {
		if ((s[0] = (char)getc(infile)) == '$')
		    labelorg = initorg;
		else
		{
		    if (s[0] != '=') ungetc(s[0],infile);
		    if (fscanf(infile,"%d",&labelorg) != 1)
                        gt_abort(">E dretog: invalid $=# command\n");
                    if (labelorg < 0)
                        gt_abort(">E dretog: must have labelorg >= 0\n");
		}
            }
	    else if (s[0] == 'g' || (s[0] >= '0' && s[0] <= '9')
                                                       || s[0] == ';')
	    {
		if (n < 0)
		    gt_abort(">E dretog: g command before n is defined\n");
		if (s[0] != 'g') ungetc(s[0],infile);
		m = (n + WORDSIZE - 1) / WORDSIZE;
#if MAXN
		if (n > MAXN || m > MAXM)
		    gt_abort(">E n or m too big\n");
#else
		DYNALLOC2(graph,g,g_sz,n,m,"dretog");
#endif
		++nin;
		readgraph(infile,g,digraph,prompt,FALSE,78,m,n);

		if (outcode == DIGRAPH6)
                    writed6(outfile,g,m,n);
                else
	        {
		    if (digraph) fprintf(stderr,
                           ">W writing digraph in undirected format\n");
		    if (outcode == SPARSE6) writes6(outfile,g,m,n);
		    else                    writeg6(outfile,g,m,n);
		}
	    }
	    else if (s[0] == 'q')
		exit(0);
	    else
	    {
	 	fprintf(stderr,">E dretog: invalid command \"%c\"\n",s[0]);
		gt_abort(NULL);
	    }
	}

	if (!qswitch)
	    fprintf(stderr,">Z  " COUNTER_FMT
                               " graphs converted from %s to %s\n",
			    nin,infilename,outfilename);

	exit(0);
}
