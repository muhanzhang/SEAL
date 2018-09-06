/* labelg.c version 2.0; B D McKay, Jun 2015 */

#define USAGE "labelg [-q] [sg | C#W#] [-fxxx] [-S|-t] [-i# -I#:# -K#] [infile [outfile]]"

#define HELPTEXT \
" Canonically label a file of graphs or digraphs.\n\
\n\
    -s  force output to sparse6 format\n\
    -g  force output to graph6 format\n\
    -z  force output to digraph6 format\n\
        If neither -s, -g or -z are given, the output format is\n\
        determined by the header or, if there is none, by the\n\
        format of the first input graph. As an exception, digraphs\n\
        are always written in digraph6 format.\n\
    -S  Use sparse representation internally.\n\
         Note that this changes the canonical labelling.\n\
         Multiple edges are not supported.  One loop per vertex is ok.\n\
    -t  Use Traces.\n\
         Note that this changes the canonical labelling.\n\
         Multiple edges and loops are not supported, nor invariants.\n\
\n\
    -C# Make an invariant in 0..#-1 and output the number of graphs\n\
	with each value of the invariant.  Don't write graphs unless\n\
        -W too.\n\
    -W# (requires -C) Output the graphs with this invariant value,\n\
        in their original labelling.  Don't write the table.\n\
\n\
    The output file will have a header if and only if the input file does.\n\
\n\
    -fxxx  Specify a partition of the point set.  xxx is any\n\
        string of ASCII characters except nul.  This string is\n\
        considered extended to infinity on the right with the\n\
        character 'z'.  One character is associated with each point,\n\
        in the order given.  The labelling used obeys these rules:\n\
         (1) the new order of the points is such that the associated\n\
        characters are in ASCII ascending order\n\
         (2) if two graphs are labelled using the same string xxx,\n\
        the output graphs are identical iff there is an\n\
        associated-character-preserving isomorphism between them.\n\
        No option can be concatenated to the right of -f.\n\
\n\
    -i#  select an invariant (1 = twopaths, 2 = adjtriang(K), 3 = triples,\n\
        4 = quadruples, 5 = celltrips, 6 = cellquads, 7 = cellquins,\n\
        8 = distances(K), 9 = indsets(K), 10 = cliques(K), 11 = cellcliq(K),\n\
       12 = cellind(K), 13 = adjacencies, 14 = cellfano, 15 = cellfano2,\n\
       16 = refinvar(K))\n\
    -I#:#  select mininvarlevel and maxinvarlevel (default 1:1)\n\
    -K#   select invararg (default 3)\n\
\n\
    -q  suppress auxiliary information\n"


/*************************************************************************/

#include "gtools.h"
#include "nautinv.h"
#include "gutils.h"
#include "traces.h"

static struct invarrec
{
    void (*entrypoint)(graph*,int*,int*,int,int,int,int*,
                      int,boolean,int,int);
    void (*entrypoint_sg)(graph*,int*,int*,int,int,int,int*,
                      int,boolean,int,int);
    char *name;
} invarproc[]
    = {{NULL, NULL, "none"},
       {twopaths,   NULL, "twopaths"},
       {adjtriang,  NULL, "adjtriang"},
       {triples,    NULL, "triples"},
       {quadruples, NULL, "quadruples"},
       {celltrips,  NULL, "celltrips"},
       {cellquads,  NULL, "cellquads"},
       {cellquins,  NULL, "cellquins"},
       {distances, distances_sg, "distances"},
       {indsets,    NULL, "indsets"},
       {cliques,    NULL, "cliques"},
       {cellcliq,   NULL, "cellcliq"},
       {cellind,    NULL, "cellind"},
       {adjacencies, adjacencies_sg, "adjacencies"},
       {cellfano,   NULL, "cellfano"},
       {cellfano2,  NULL, "cellfano2"},
       {refinvar,   NULL, "refinvar"}
      };

#define NUMINVARS ((int)(sizeof(invarproc)/sizeof(struct invarrec)))

static nauty_counter orbtotal;
static double unorbtotal;
extern int gt_numorbits;

/**************************************************************************/

int
main(int argc, char *argv[])
{
	graph *g;
	sparsegraph sg,sh;
	int m,n,codetype;
	int argnum,j,outcode;
	char *arg,sw,*fmt;
	boolean badargs,digraph;
	boolean sswitch,gswitch,qswitch,fswitch,Oswitch;
	boolean iswitch,Iswitch,Kswitch,Mswitch,Sswitch;
	boolean uswitch,tswitch,Cswitch,Wswitch,zswitch;
	boolean dooutput;
	int tabsize,outinvar;
	int inv,mininvarlevel,maxinvarlevel,invararg;
	long minil,maxil;
	double t;
	char *infilename,*outfilename;
	FILE *infile,*outfile;
	nauty_counter nin,nout;
	int ii,secret,loops;
	DEFAULTOPTIONS_TRACES(traces_opts);
	TracesStats traces_stats;
#if MAXN
	graph h[MAXN*MAXM];
	int lab[MAXN],ptn[MAXN],orbits[MAXN];
#else
	DYNALLSTAT(graph,h,h_sz);
	DYNALLSTAT(int,lab,lab_sz);
	DYNALLSTAT(int,ptn,ptn_sz);
	DYNALLSTAT(int,orbits,orbits_sz);
#endif
        DYNALLSTAT(nauty_counter,tab,tab_sz);

	HELP; PUTVERSION;

	nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

	sswitch = gswitch = qswitch = FALSE;
	fswitch = Oswitch = Mswitch = FALSE;
	iswitch = Iswitch = Kswitch = FALSE;
	uswitch = Sswitch = tswitch = FALSE;
	zswitch = Cswitch = Wswitch = FALSE;
	infilename = outfilename = NULL;
	inv = 0;

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
		    else SWBOOLEAN('u',uswitch)
		    else SWBOOLEAN('q',qswitch)
		    else SWBOOLEAN('O',Oswitch)
		    else SWBOOLEAN('S',Sswitch)
		    else SWBOOLEAN('t',tswitch)
		    else SWINT('C',Cswitch,tabsize,"labelg -C")
		    else SWINT('W',Wswitch,outinvar,"labelg -W")
		    else SWINT('i',iswitch,inv,"labelg -i")
		    else SWINT('K',Kswitch,invararg,"labelg -K")
		    else SWRANGE('k',":-",Iswitch,minil,maxil,"labelg -k")
		    else SWRANGE('I',":-",Iswitch,minil,maxil,"labelg -I")
		    else if (sw == 'f')
		    {
			fswitch = TRUE;
			fmt = arg;
			break;
		    }
		    else SWINT('M',Mswitch,secret,"labelg -M")
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

	if ((sswitch!=0) + (zswitch!=0) + (gswitch!=0) + (uswitch!=0) > 1)
            gt_abort(">E labelg: -szgu are incompatible\n");

	if (tswitch && (fswitch || Sswitch))
            gt_abort(">E labelg: -t is incompatible with -S and -f \n");

	if (iswitch && inv == 0) iswitch = FALSE;

	if (!Sswitch && iswitch && (inv > NUMINVARS))
	    gt_abort(">E labelg: -i value must be 0..16\n");
	if (tswitch && iswitch)
	    gt_abort(">E labelg: invariants are not available with -t\n");
	if (Wswitch && (sswitch || gswitch || zswitch || uswitch || !Cswitch))
            gt_abort(">E labelg: -W forbids -sgzu and requires -C\n");
	if (Cswitch && tabsize <= 0)
	    gt_abort(">E labelg: bad value for -C\n");
	if (Wswitch && (outinvar < 0 || outinvar >= tabsize)) 
	    gt_abort(">E labelg: value of -W is not valid\n");
        if (Sswitch && iswitch && invarproc[inv].entrypoint_sg == NULL)
	    gt_abort(
                ">E labelg: that invariant is not available in sparse mode\n");


	if (iswitch)
	{
	    if (Iswitch)
	    {
		mininvarlevel = minil;
		maxinvarlevel = maxil;
	    }
	    else
		mininvarlevel = maxinvarlevel = 1;
	    if (!Kswitch) invararg = 3;
	}

	if (!Mswitch) secret = 1;

	if (badargs || argnum > 2)
	{
	    fprintf(stderr,">E Usage: %s\n",USAGE);
	    GETHELP;
	    exit(1);
	}

	if (!qswitch)
	{
	    fprintf(stderr,">A labelg");
	    if (sswitch || gswitch || fswitch || iswitch || zswitch
			|| tswitch || Sswitch || Cswitch || Wswitch)
		fprintf(stderr," -");
	    if (sswitch) fprintf(stderr,"s");
	    if (gswitch) fprintf(stderr,"g");
	    if (zswitch) fprintf(stderr,"z");
	    if (Sswitch) fprintf(stderr,"S");
	    if (tswitch) fprintf(stderr,"t");
	    if (Cswitch) fprintf(stderr,"C%d",tabsize);
	    if (Wswitch) fprintf(stderr,"W%d",outinvar);
	    if (iswitch)
		fprintf(stderr,"i=%s[%d:%d,%d]",invarproc[inv].name,
		        mininvarlevel,maxinvarlevel,invararg);
	    if (fswitch) fprintf(stderr,"f%s",fmt);
	    if (argnum > 0) fprintf(stderr," %s",infilename);
	    if (argnum > 1) fprintf(stderr," %s",outfilename);
	    fprintf(stderr,"\n");
	    fflush(stderr);
	}

	if (Cswitch & !Wswitch)
	{
	    DYNALLOC1(nauty_counter,tab,tab_sz,tabsize,"table");
	    for (j = 0; j < tabsize; ++j) tab[j] = 0;
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

	if (uswitch)
	    outcode = 0;
	else if (Cswitch && !Wswitch)
	    outcode = -1;
	else if (Cswitch && Wswitch)
	    outcode = -2;
        else if (gswitch)
            outcode = GRAPH6;
        else if (sswitch)
            outcode = SPARSE6;
        else if (zswitch)
            outcode = DIGRAPH6;
	else if ((codetype&GRAPH6))
	    outcode = GRAPH6;
	else if ((codetype&SPARSE6))
	    outcode = SPARSE6;
	else if ((codetype&DIGRAPH6))
	    outcode = DIGRAPH6;
	else 
        {
	    outcode = GRAPH6;
	    fprintf(stderr,
                ">W labelg doesn't handle this graph format, writing graph6.\n");
        }

	dooutput = (outcode == GRAPH6 || outcode == DIGRAPH6 || outcode == SPARSE6);

	if (!fswitch) fmt = NULL;

	if (!uswitch && (codetype&HAS_HEADER))
	{
	    if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
	    else    		    writeline(outfile,GRAPH6_HEADER);
	}

	nin = nout = 0;
	orbtotal = 0;
	unorbtotal = 0.0;
	t = CPUTIME;
	if (Sswitch)
	{
	    SG_INIT(sg);
	    SG_INIT(sh);
	    while (TRUE)
	    {
		if (read_sgg_loops(infile,&sg,&loops,&digraph) == NULL) break;
		++nin;
		n = sg.nv;
		m = (n + WORDSIZE - 1) / WORDSIZE;
		SG_ALLOC(sh,n,sg.nde,"labelg");
		for (ii = 0; ii < secret; ++ii)
		    fcanonise_inv_sg(&sg,m,n,&sh,fmt,invarproc[inv].entrypoint_sg,
		                 mininvarlevel,maxinvarlevel,invararg,loops>0||digraph);
		sortlists_sg(&sh);
	        orbtotal += gt_numorbits;
	        unorbtotal += 1.0 / gt_numorbits;
		if (dooutput && (outcode == DIGRAPH6 || digraph)) writed6_sg(outfile,&sh);
	        else if (outcode == SPARSE6) writes6_sg(outfile,&sh);
	        else if (outcode == GRAPH6) writeg6_sg(outfile,&sh);
		else if (outcode == -1)
		    ++tab[hashgraph_sg(&sh,137)%tabsize];
		else if (outcode == -2 && (hashgraph_sg(&sh,137)%tabsize) == outinvar)
		{
		    writelast(outfile);
		    ++nout;
		}
	    }
	}
	else if (tswitch)
	{
	    SG_INIT(sg);
	    SG_INIT(sh);
	    traces_opts.getcanon = TRUE;
	    traces_opts.writeautoms = FALSE;
	    traces_opts.verbosity = 0;
	    traces_opts.outfile = stdout;
	    
	    while (TRUE)
	    {
		if (read_sgg_loops(infile,&sg,&loops,&digraph) == NULL) break;
		if (loops > 0 || digraph) gt_abort(">E Traces does not allow loops or directed edges\n");
		++nin;
		n = sg.nv;
	        DYNALLOC1(int,lab,lab_sz,n,"traces@labelg");
	        DYNALLOC1(int,ptn,ptn_sz,n,"traces@labelg");
	        DYNALLOC1(int,orbits,orbits_sz,n,"traces@labelg");
		SG_ALLOC(sh,n,sg.nde,"labelg");
		for (ii = 0; ii < n; ++ii) { lab[ii] = ii; ptn[ii] = 1; }
		ptn[n-1] = 0;
		for (ii = 0; ii < secret; ++ii)
		    Traces(&sg,lab,ptn,orbits,&traces_opts,&traces_stats,&sh);
		sortlists_sg(&sh);
	        orbtotal += gt_numorbits;
	        unorbtotal += 1.0 / gt_numorbits;
	        if (outcode == SPARSE6) writes6_sg(outfile,&sh);
	        else if (outcode == GRAPH6) writeg6_sg(outfile,&sh);
		else if (outcode == -1)
		    ++tab[hashgraph_sg(&sh,137)%tabsize];
		else if (outcode == -2 && (hashgraph_sg(&sh,137)%tabsize) == outinvar)
		{
		    writelast(outfile);
		    ++nout;
		}
		SG_FREE(sh);
	    }
	}
	else
	{
	    while (TRUE)
	    {
	        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
	        ++nin;
#if !MAXN
	        DYNALLOC2(graph,h,h_sz,n,m,"labelg");
#endif
		loops = loopcount(g,m,n);
	        for (ii = 0; ii < secret; ++ii)
	            fcanonise_inv(g,m,n,h,fmt,invarproc[inv].entrypoint,
			        mininvarlevel,maxinvarlevel,invararg,loops>0||digraph);
	        orbtotal += gt_numorbits;
	        unorbtotal += 1.0 / gt_numorbits;
		if (dooutput && (outcode == DIGRAPH6 || digraph)) writed6(outfile,h,m,n);
	        else if (outcode == SPARSE6) writes6(outfile,h,m,n);
	        else if (outcode == GRAPH6) writeg6(outfile,h,m,n);
		else if (outcode == -1)
		    ++tab[hashgraph(h,m,n,137)%tabsize];
		else if (outcode == -2 && (hashgraph(h,m,n,137)%tabsize) == outinvar)
		{
		    writelast(outfile);
		    ++nout;
		}
	        FREES(g);
	    }
	}
	t = CPUTIME - t;

	if (Oswitch)
	    fprintf(stderr,">C orbit totals = " COUNTER_FMT " %15.8f\n",
			   orbtotal,unorbtotal);
        if (!qswitch)
	{
	    if (outcode == -1)
	    {
		for (j = 0; j < tabsize; ++j)
		    fprintf(stderr,">C %5d " COUNTER_FMT "\n",j,tab[j]);
		fprintf(stderr,">Z " COUNTER_FMT " total from %s in %3.2f sec.\n",
		    nin,infilename,t);
	    }
	    else if (outcode == -2)
	    {
                fprintf(stderr,
                    ">Z " COUNTER_FMT " out of " COUNTER_FMT
                    " graphs copied from %s to %s in %3.2f sec.\n",
                    nout,nin,infilename,outfilename,t);
	    }
	    else
	    {
                fprintf(stderr,
                        ">Z " COUNTER_FMT
                        " graphs labelled from %s to %s in %3.2f sec.\n",
                    nin,infilename,outfilename,t);
	    }
	}

	exit(0);
}
