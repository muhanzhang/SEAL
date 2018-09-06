/* checks6.c; May 2005 */

#define USAGE "checks6 [-p#:#w] [infile [outfile]]"

#define HELPTEXT \
"  Check a file of graphs, optionally write corrected version\n\
\n\
     -p# -p#:#  \n\
         Specify range of input lines (first is 1)\n\
\n\
     -w  Write corrected graphs (default is not to write)\n\
         A header is written if there is one in the input.\n"

/***********************************************************************/

#include "gtools.h"

/***********************************************************************/

static boolean
seemsbad(char *s)
/* Check graph string for apparent problem, if so, correct it */
{
        int i,j,k,m,n;
        char *p,x,pq;
	set *gj;
	long ii;
	int r,rr,topbit,nb,lastj;
	graph g[16];

	if (s[0] != ':') return FALSE;  /* not sparse6 */

	n = graphsize(s);
	if (n != 2 && n != 4 && n != 8 && n != 16) return FALSE;
	m = 1;

	stringtograph(s,g,m);
	if (g[n-1] != bit[n-1]) return FALSE;
	if (g[n-2] == 0) return FALSE;

	g[n-1] = 0;
	p = s+2;

	for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb) {}
	topbit = 1 << (nb-1);
	k = 6;
	x = 0;

	lastj = 0;
	for (j = 0; j < n; ++j)
	{
	    gj = GRAPHROW(g,j,m);
	    for (i = 0; i <= j; ++i)
	    {
		if (ISELEMENT(gj,i))
		{
		    if (j == lastj)
		    {
		        x <<= 1;
		        if (--k == 0)
		        {
		            p++;
		            k = 6;
		            x = 0;
		        }
		    }
		    else
		    {
			x = (x << 1) | 1;
			if (--k == 0)
			{
		            p++;
			    k = 6;
			    x = 0;
			}
			if (j > lastj+1)
			{
			    for (r = 0, rr = j; r < nb; ++r, rr <<= 1)
			    {
			        if (rr & topbit) x = (x << 1) | 1;
			        else             x <<= 1;
			        if (--k == 0)
			        {
				    p++;
				    k = 6;
				    x = 0;
			        }
			    }
			    x <<= 1;
			    if (--k == 0)
			    {
				p++;
				k = 6;
				x = 0;
			    }
			}
			lastj = j;
		    }
		    for (r = 0, rr = i; r < nb; ++r, rr <<= 1)
		    {
			if (rr & topbit) x = (x << 1) | 1;
			else             x <<= 1;
			if (--k == 0)
			{
			    p++;
			    k = 6;
			    x = 0;
			}
		    }
		}
	    }
	}

        if (k != 6)
        {
            if (k >= nb+1 && lastj == n-2 && n == (1<<nb))
	    {
                *p++ = BIAS6 + ((x << k) | ((1 << (k-1)) - 1));
	    	return TRUE;
	    }
            else
                return FALSE;
        }
}

/***********************************************************************/

int
main(int argc, char *argv[])
{
	int m,n,codetype;
	char *infilename,*outfilename;
	FILE *infile,*outfile;
	int outcode;
	long nin,nerr;
	int argnum,j;
	char *arg,sw,*s;
	boolean pswitch;
	boolean wswitch;
	boolean badargs;
	long pval1,pval2,maxin;

	HELP; PUTVERSION;

	wswitch = pswitch = FALSE;
	infilename = outfilename = NULL;

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
		         SWBOOLEAN('w',wswitch)
		    else SWRANGE('p',":-",pswitch,pval1,pval2,"checks6 -p")
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

	if (badargs || argnum > 2)
	{
	    fprintf(stderr,">E Usage: %s\n",USAGE);
	    GETHELP;
	    exit(1);
	}

	if (infilename && infilename[0] == '-') infilename = NULL;
	infile = opengraphfile(infilename,&codetype,FALSE,
			       pswitch ? pval1 : 1);
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

	if (wswitch && (codetype&HAS_HEADER))
	{
	    if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
	    else    		    writeline(outfile,GRAPH6_HEADER);
	}

	nerr = nin = 0;
	if (!pswitch || pval2 == NOLIMIT)
	    maxin = NOLIMIT;
	else if (pval1 < 1) maxin = pval2;
	else                maxin = pval2 - pval1 + 1;
	while (nin < maxin || maxin == NOLIMIT)
	{
	    if ((s = gtools_getline(infile)) == NULL) break;
	    ++nin;

	    if (seemsbad(s)) ++nerr;
	    if (wswitch) writeline(outfile,s);
	}

	fprintf(stderr,">Z  %ld graphs read",nin);
	if (nerr > 0) fprintf(stderr,"; %ld probable errors",nerr);
	else          fprintf(stderr,"; NO PROBLEMS");
	if (wswitch) fprintf(stderr,"; %ld graphs written",nin);
	fprintf(stderr,"\n");

	exit(0);
}
