/* vcolg.c version 1.0; B D McKay, Aug 31, 2013 */

#define USAGE \
"vcolg [-q] [-u|-T|-G|-A|-B] [-e#|-e#:#] \n" \
"       [-m#] [-f#] [-D#|-r#|-l#] [infile [outfile]]"

#define HELPTEXT \
" Read undirected loop-free graphs and colour their vertices in\n\
  in all possible ways with colours 0,1,2,... .\n\
  Isomorphic graphs derived from the same input are suppressed.\n\
  If the input graphs are non-isomorphic then the output graphs are also.\n\
\n\
    -e# | -e#:#  specify a value or range of the total value of the colours\n\
    -m# number of available colours (default 2)\n\
    -f# Use the group that fixes the first # vertices setwise\n\
    -T  use a simple text output format (nv ne {col} {v1 v2})\n\
    -u  no output, just count them\n\
    -q  suppress auxiliary information\n"

/*************************************************************************/

#include "gtools.h"
#include "naugroup.h"

nauty_counter vc_nin,vc_nout;
FILE *outfile;

#define MAXNV 128 
#define MAXNE 1024

static int col[MAXNV];
static boolean first;
static int lastreject[MAXNV];
static boolean lastrejok;
static unsigned long groupsize;
static unsigned long newgroupsize;
static boolean Tswitch;

static int fail_level;

#define GROUPTEST_NOT 
#ifdef GROUPTEST
static long long totallab;
#endif

/* If OUTPROC is defined at compile time, and -u is not used, the
 * procedure OUTPROC is called for each graph.  This must be linked
 * by the compiler.  The arguments are
 * f = open output file
 * g = the input graph
 * col[0..n-1] = the colours
 * m,n = usual nauty meanings
 */

/* SUMMARY feature
 *
 * If SUMMARY is defined, it must expand as the name of a procedure
 * with prototype  void SUMMARY(void).  It is called at the end before
 * the normal summary (which can be suppressed with -q).  The numbers of
 * graphs read and coloured graphs produced are available in the global
 * variables vc_nin and vc_nout (type nauty_counter).
 */

#ifdef OUTPROC
extern void OUTPROC(FILE*,graph*,int*,int,int);
#endif

#ifdef SUMMARY
extern void SUMMARY(void);
#endif

/**************************************************************************/

static void
writeautom(int *p, int n)
/* Called by allgroup. */
{
    int i;

    for (i = 0; i < n; ++i) printf(" %2d",p[i]);
    printf("\n");
}

/**************************************************************************/

static int
ismax(int *p, int n)
/* test if col^p <= col */
{
    int i,k;
    int fail;

    fail = 0;
    for (i = 0; i < n; ++i)
    {
	k = p[i];
	if (k > fail) fail = k;
        if (col[k] > col[i])
        {
            fail_level = fail;
            return FALSE;
        }
        else if (col[k] < col[i]) return TRUE;
    }

    ++newgroupsize;
    return TRUE;
}

/**************************************************************************/

static void
testmax(int *p, int n, int *abort)
/* Called by allgroup2. */
{
    int i;

    if (first)
    {                       /* only the identity */
        first = FALSE;
        return;
    }

    if (!ismax(p,n))
    {
        *abort = 1;
        for (i = 0; i < n; ++i) lastreject[i] = p[i];
        lastrejok = TRUE;
    }
}

/**************************************************************************/

static int
trythisone(grouprec *group, graph *g, int m, int n)
/* Try one solution, accept if maximal. */
/* Return value is level to return to. */
{
    int i,j;
    boolean accept;
    graph *gi;
    size_t ne;

    newgroupsize = 1;

    if (!group || groupsize == 1)
        accept = TRUE;
    else if (lastrejok && !ismax(lastreject,n))
        accept = FALSE;
    else if (lastrejok && groupsize == 2)
        accept = TRUE;
    else
    {
        newgroupsize = 1;
        first = TRUE;

        if (allgroup2(group,testmax) == 0)
            accept = TRUE;
        else
            accept = FALSE;
    }

    if (accept)
    {
#ifdef GROUPTEST
        if (groupsize % newgroupsize != 0)
                    gt_abort("group size error\n");
        totallab += groupsize/newgroupsize;
#endif

        ++vc_nout;

        if (outfile)
        {
#ifdef OUTPROC
            OUTPROC(outfile,g,col,m,n);
#else
	    ne = 0;
	    for (gi = g + m*(size_t)n; --gi >= g; )
		ne += POPCOUNT(*gi);
	    ne /= 2;
            fprintf(outfile,"%d %lu",n,(unsigned long)ne);
    
	    for (i = 0; i < n; ++i) fprintf(outfile," %d",col[i]);
	    fprintf(outfile," ");
	    for (i = 0, gi = g; i < n; ++i, gi += m)
	    {
		for (j = i; (j = nextelement(gi,m,j)) >= 0; )
                    fprintf(outfile," %d %d",i,j);
	    }
            fprintf(outfile,"\n");
#endif
        }
        return n-1;
    }
    else
        return fail_level-1;
}

/**************************************************************************/

static int
scan(int level, graph *g, int *prev, long minedges, long maxedges,
    long sofar, long numcols, grouprec *group, int m, int n)
/* Recursive scan for default case */
/* Returned value is level to return to. */
{
    int left;
    long min,max,k,ret;

    if (level == n)
        return trythisone(group,g,m,n);

    left = n - level - 1;
    min = minedges - sofar - numcols*left;
    if (min < 0) min = 0;
    max = maxedges - sofar;
    if (max >= numcols) max = numcols - 1;
    if (prev[level] >= 0 && col[prev[level]] < max)
	max = col[prev[level]];

    for (k = min; k <= max; ++k)
    {
        col[level] = k;
        ret = scan(level+1,g,prev,minedges,maxedges,sofar+k,numcols,group,m,n);
	if (ret < level) return ret;
    }

    return level-1;
}

/**************************************************************************/

#define SORT_OF_SORT 3
#define SORT_NAME sortwt
#define SORT_TYPE1 int
#define SORT_TYPE2 int
#include "sorttemplates.c"

static void
setlab(int *weight, int *lab, int *ptn, int n)
/* Define (lab,ptn) according to weights. */
{
    int i;

    for (i = 0; i < n; ++i) lab[i] = i;

    sortwt(lab,weight,n);
    for (i = 0; i < n-1; ++i)
    {
        if (weight[lab[i]] != weight[lab[i+1]])
	    ptn[i] = 0;
	else
	    ptn[i] = 1;
    }
    ptn[n-1] = 0;
}

/**************************************************************************/

static void
colourit(graph *g, int nfixed, long minedges, long maxedges,
         long numcols, int m, int n)
{
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[100];
    grouprec *group;
    int i,j,k;
    set *gi,*gj;
    int lab[MAXNV],ptn[MAXNV],orbits[MAXNV],deg[MAXNV];
    int prev[MAXNV]; /* If >= 0, earlier point that must have greater colour */
    int weight[MAXNV];
    int last0,lastn,thisdeg,region,start,stop;

    if (n > MAXNV) gt_abort(">E vcolg: MAXNV exceeded\n");
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    for (region = 0; region < 2; ++region)
    {
	if (region == 0)
	{
	    if (nfixed == 0) continue;
	    start = 0;
	    stop = nfixed;
	    if (stop > n) stop = n;
	}
	else
	{
	    if (nfixed >= n) continue;
	    start = nfixed;
	    stop = n;
	}
	
	for (i = start, gi = g + m*(size_t)start; i < stop; ++i, gi += m)
	{
	    for (j = i-1, gj = gi-m; j >= start; --j, gj -= m)
	    {
		for (k = 0; k < m; ++k) if (gi[k] != gj[k]) break;
		if (k < m)
		{
		    FLIPELEMENT(gi,i); FLIPELEMENT(gj,j);
		    for (k = 0; k < m; ++k) if (gi[k] != gj[k]) break;
		    FLIPELEMENT(gi,i); FLIPELEMENT(gj,j);
		}
		if (k == m) break;
	    }
	    if (j >= start)
	    {
		prev[i] = j;
		weight[i] = weight[j] + 1;
	    }
	    else
	    {
		prev[i] = -1;
		weight[i] = 0;
	    }
	}
    }

    for (i = nfixed; i < n; ++i) weight[i] += nfixed;

    if (maxedges == NOLIMIT || maxedges > n*numcols) maxedges = n*numcols;
    if (n*numcols < minedges) return;

    options.userautomproc = groupautomproc;
    options.userlevelproc = grouplevelproc;
    options.defaultptn = FALSE;

    setlab(weight,lab,ptn,n);
 
    nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,100,m,n,NULL);

    if (stats.grpsize2 == 0)
        groupsize = stats.grpsize1 + 0.1;
    else
        groupsize = 0;

    group = groupptr(FALSE);
    makecosetreps(group);

    if (stats.numorbits < n)
    {
	j = n;
	for (i = 0; i < n; ++i)
	    if (orbits[i] < i && orbits[i] < j) j = orbits[i];

	for (i = j + 1; i < n; ++i)
	    if (orbits[i] == j) prev[i] = j;
    }

    lastrejok = FALSE;
    for (i = 0; i < n; ++i) col[i] = 0;

    scan(0,g,prev,minedges,maxedges,0,numcols,group,m,n);
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    graph *g;
    int i,m,n,codetype;
    int argnum,j,nfixed;
    char *arg,sw;
    boolean badargs;
    boolean fswitch,uswitch,eswitch,qswitch,mswitch;
    long minedges,maxedges,numcols;
    double t;
    char *infilename,*outfilename;
    FILE *infile;
    char msg[201];
    int msglen;

    HELP; PUTVERSION;

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    fswitch = Tswitch = FALSE;
    uswitch = eswitch = mswitch = qswitch = FALSE;
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
                     SWLONG('m',mswitch,numcols,"vcolg -m")
                else SWBOOLEAN('q',qswitch)
                else SWBOOLEAN('u',uswitch)
                else SWBOOLEAN('T',Tswitch)
                else SWINT('f',fswitch,nfixed,"vcolg -f")
                else SWRANGE('e',":-",eswitch,minedges,maxedges,"vcolg -e")
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

    if ((Tswitch!=0) + (uswitch!=0) >= 2)
        gt_abort(">E vcolg: -T and -u are incompatible\n");

#ifndef OUTPROC
    if (!Tswitch && !uswitch)
        gt_abort(">E vcolg: must use -T or -u\n");
#endif

    if (!mswitch) numcols = 2;
    if (!eswitch)
    {
        minedges = 0;
        maxedges = NOLIMIT;
    }
    if (!fswitch) nfixed = 0;

    if (!qswitch)
    {
        msg[0] = '\0';
        CATMSG0(">A vcolg");
        if (eswitch || mswitch || uswitch || (fswitch && nfixed > 0)
              || Tswitch)
            CATMSG0(" -");
        if (mswitch) CATMSG1("m%ld",numcols);
        if (uswitch) CATMSG0("u");
        if (Tswitch) CATMSG0("T");
        if (fswitch) CATMSG1("f%d",nfixed);
        if (eswitch) CATMSG2("e%ld:%ld",minedges,maxedges);
        msglen = strlen(msg);
        if (argnum > 0) msglen += strlen(infilename);
        if (argnum > 1) msglen += strlen(outfilename);
        if (msglen >= 196)
        {
            fputs(msg,stderr);
            if (argnum > 0) fprintf(stderr," %s",infilename);
            if (argnum > 1) fprintf(stderr," %s",outfilename);
            fprintf(stderr,"\n");
        }
        else
        {
            if (argnum > 0) CATMSG1(" %s",infilename);
            if (argnum > 1) CATMSG1(" %s",outfilename);
            CATMSG0("\n");
            fputs(msg,stderr);
        }
        fflush(stderr);
    }

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,FALSE,1);
    if (!infile) exit(1);
    if (!infilename) infilename = "stdin";

    NODIGRAPHSYET(codetype);

    if (uswitch)
        outfile = NULL;
    else
    {
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
    }

    vc_nin = vc_nout = 0;

    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readg(infile,NULL,0,&m,&n)) == NULL) break;
        ++vc_nin;

        colourit(g,nfixed,minedges,maxedges,numcols,m,n);
	if (!uswitch && ferror(outfile))
	    gt_abort(">E vcolg output error\n");
        FREES(g);
    }
    t = CPUTIME - t;

#ifdef SUMMARY
    SUMMARY();
#endif

    if (!qswitch)
    {
        fprintf(stderr,">Z ");
        PRINT_COUNTER(stderr,vc_nin);
        fprintf(stderr," graphs read from %s",infilename);
        fprintf(stderr,"; ");
        PRINT_COUNTER(stderr,vc_nout);
        if (!uswitch)
            fprintf(stderr," coloured graphs written to %s",outfilename);
        else
            fprintf(stderr," coloured graphs generated");
        fprintf(stderr,"; %.2f sec\n",t);
    }

#ifdef GROUPTEST
    fprintf(stderr,"Group test = %lld\n",totallab);
#endif

    exit(0);
}
