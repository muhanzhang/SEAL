/* directg.c version 2.0; B D McKay, June 19, 2015 */

#define USAGE \
  "directg [-q] [-u|-T|-G] [-V] [-o] [-f#] [-e#|-e#:#] [-s#/#] [infile [outfile]]"

#define HELPTEXT \
" Read undirected graphs and orient their edges in all possible ways.\n\
  Edges can be oriented in either or both directions (3 possibilities).\n\
  Isomorphic directed graphs derived from the same input are suppressed.\n\
  If the input graphs are non-isomorphic then the output graphs are also.\n\
\n\
   -e# | -e#:#  specify a value or range of the total number of arcs\n\
    -o     orient each edge in only one direction, never both\n\
   -f#  Use only the subgroup that fixes the first # vertices setwise\n\
\n\
    -T  use a simple text output format (nv ne edges) instead of digraph6\n\
    -G  like -T but includes group size as third item (if less than 10^10)\n\
          The group size does not include exchange of isolated vertices.\n\
    -V  only output graphs with nontrivial groups (including exchange of\n\
          isolated vertices).  The -f option is respected.\n\
    -u  no output, just count them\n\
    -s#/#  Make only a fraction of the orientations: The first integer is\n\
            the part number (first is 0) and the second is the number of\n\
            parts. Splitting is done per input graph independently.\n\
    -q  suppress auxiliary information\n"

/* NOTE: Do not use nauty1.* for this program! */

/*************************************************************************/

#include "gtools.h"
#include "naugroup.h"

nauty_counter dg_nin,dg_ngen,dg_nout,dg_skipped; 

FILE *outfile;

#define MAXNV 128 
#define MAXNE 1024
static int v0[MAXNE],v1[MAXNE];
static int edgeno[MAXNV][MAXNV];

#define MAXME ((2*MAXNE+WORDSIZE-1)/WORDSIZE)

static set x[MAXME];
static int ix[2*MAXNE],nix;
static boolean first;
static int me;
static int lastreject[MAXNV];
static boolean lastrejok;
static int rejectlevel;
static unsigned long groupsize;
static unsigned long newgroupsize;
static boolean Gswitch,Vswitch,ntgroup,ntisol,digraph6;

static int splitlevel,splitmod,splitres,splitcount;
static unsigned long splitcases;

/* DEGPRUNE feature 
 *
 * If DEGPRUNE is defined it must have a value equal to the name of a 
 * procedure to be supplied by the user and linked to this program.
 * The prototype must be
 *     int DEGPRUNE(int *indeg, int outdeg*, int v, int n)
 * Here n is the number of vertices altogether, and v (0..n-1) is the
 * number of one vertex.  At this point in the program, some directed
 * edges have been inserted, and the indegrees and outdegrees have the
 * values given in indeg[] and outdeg[].  Moreover, it is known that
 * no further edges will be added to or from v, so indeg[v] and outdeg[v]
 * are final.  If DEGPRUNE returns a non-zero value, this branch of the
 * search will be abandoned.
 * Before any graph is output, DEGPRUNE will have been called for every
 * vertex, but it cannot be assumed that DEGPRUNE will be called in order
 * of vertex number.
 */

/* INPUTGRAPH feature
 *
 * If INPUTGRAPH is defined, it must expand as the name of a procedure
 * with prototype like  int INPUTGRAPH(graph *g, int m, int n).
 * This procedure will be called for each input graph before it is 
 * processed. The graph will be skipped if the value 0 is returned.
 */

/* PROCESS feature
 *
 * If PROCESS is defined, it must expand as the name of a procedure
 * with prototype like void PROCESS(FILE *f, graph *g, int n). This
 * procedure will be called for each output graph before it is output,
 * with f being the output file (possibly NULL).
 * It is an error if n > WORDSIZE.
 */

/* SUMMARY feature
 *
 * If SUMMARY is defined, it must expand as the name of a procedure
 * with prototype  void SUMMARY(void).  It is called at the end before
 * the normal summary (which can be suppressed with -q).  The numbers of
 * graphs read and digraphs produced are available in the global variables
 * dg_nin and dg_nout (type nauty_counter).
 */

#ifdef DEGPRUNE
static int lastlev[MAXNV],indeg[MAXNV],outdeg[MAXNV];
extern int DEGPRUNE(int*,int*,int,int);
#endif

#ifdef INPUTGRAPH
extern int INPUTGRAPH(graph*,int,int);
#endif

#ifdef PROCESS
extern void PROCESS(FILE*,graph*,int);
#endif

#ifdef SUMMARY
extern void SUMMARY(void);
#endif

/* #define GROUPTEST */
#ifdef GROUPTEST
static long long totallab;
#endif

/**************************************************************************/

void
writeautom(int *p, int n)
/* Called by allgroup. */
{
    int i;

    for (i = 0; i < n; ++i) printf(" %2d",p[i]);
    printf("\n");
}

/**************************************************************************/

static boolean
ismax(int *p, int n)
/* test if x^p <= x */
{
    int i,j,k;
    set px[MAXME];

    EMPTYSET(px,me);

    for (j = 0; j < nix; ++j)
    {
        i = ix[j];
        k = i >> 1;
        if (i & 1) ADDELEMENT(px,edgeno[p[v1[k]]][p[v0[k]]]);
        else       ADDELEMENT(px,edgeno[p[v0[k]]][p[v1[k]]]);

        if (px[0] > x[0])
        {
            rejectlevel = k;
            return FALSE;
        }
    }

    rejectlevel = MAXNE+1;

    if (px[0] < x[0]) return TRUE;

    for (i = 1; i < me; ++i)
        if      (px[i] > x[i]) return FALSE;
        else if (px[i] < x[i]) return TRUE;

    ++newgroupsize;
    ntgroup = TRUE;
    return TRUE;
}

/**************************************************************************/

void
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
trythisone(grouprec *group, int ne, int n)
{
    int i,k;
    boolean accept;
    graph g[WORDSIZE];

    first = TRUE;

    ++dg_ngen;
    nix = ne;
    newgroupsize = 1;
    ntgroup = FALSE;

    if (!group || groupsize == 1)
        accept = TRUE;
    else if (lastrejok && !ismax(lastreject,n))
        accept = FALSE;
    else if (lastrejok && groupsize == 2)
        accept = TRUE;
    else
    {
        newgroupsize = 1;
        ntgroup = FALSE;
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

        if (Vswitch && !ntisol && !ntgroup) return MAXNE+1;

        ++dg_nout;

	if (digraph6)
        {
	    EMPTYSET(g,n);
	    for (i = -1; (i = nextelement(x,me,i)) >= 0; )
            {
                k = i >> 1;
                if (i & 1) g[v1[k]] |= bit[v0[k]];
                else       g[v0[k]] |= bit[v1[k]];
            }
	}

#ifdef PROCESS
	if (!digraph6)
        {
	    EMPTYSET(g,n);
	    for (i = -1; (i = nextelement(x,me,i)) >= 0; )
            {
                k = i >> 1;
                if (i & 1) g[v1[k]] |= bit[v0[k]];
                else       g[v0[k]] |= bit[v1[k]];
            }
	}
	PROCESS(outfile,g,n);
#endif

	if (digraph6)
	    writed6(outfile,g,1,n);
        else if (outfile)
        {
            fprintf(outfile,"%d %d",n,ne);
            if (Gswitch) fprintf(outfile," %lu",newgroupsize);
    
            for (i = -1; (i = nextelement(x,me,i)) >= 0; )
            {
                k = i >> 1;
                if (i & 1) fprintf(outfile," %d %d",v1[k],v0[k]);
                else       fprintf(outfile," %d %d",v0[k],v1[k]);
            }
            fprintf(outfile,"\n");
        }
        return MAXNE+1;
    }
    else
        return (rejectlevel < splitlevel ? splitlevel : rejectlevel);
}

/**************************************************************************/

static int
scan(int level, int ne, int minarcs, int maxarcs, int sofar,
    boolean oriented, grouprec *group, int n)
/* Main recursive scan; returns the level to return to. */
{
    int k,retlev;
#ifdef DEGPRUNE
    int w0,w1;

    w0 = v0[level];
    w1 = v1[level];
#endif

    if (level == splitlevel)
    {
#ifdef SPLITTEST
        ++splitcases;
        return level-1;
#endif
        if (splitcount-- != 0) return level-1;
        splitcount = splitmod - 1;
    }

    if (level == ne)
    {
        retlev = trythisone(group,sofar,n);
        return retlev;
    }

    if (oriented || sofar + 1 + 2*(ne - level - 1) >= minarcs)
    {
        k = 2*level;
        ADDELEMENT(x,k);
        ix[sofar] = k;
#ifdef DEGPRUNE
        ++outdeg[w0]; ++indeg[w1];
        if (lastlev[w0] == level && DEGPRUNE(indeg,outdeg,w0,n)
         || lastlev[w1] == level && DEGPRUNE(indeg,outdeg,w1,n))
            retlev = level;
        else
#endif
        retlev = scan(level+1,ne,minarcs,maxarcs,sofar+1,oriented,group,n);
        DELELEMENT(x,k);
#ifdef DEGPRUNE
        --outdeg[w0]; --indeg[w1];
#endif
        if (retlev < level) return retlev;
        ++k;
        ADDELEMENT(x,k);
        ix[sofar] = k;
#ifdef DEGPRUNE
        ++outdeg[w1]; ++indeg[w0];
        if (lastlev[w0] == level && DEGPRUNE(indeg,outdeg,w0,n)
         || lastlev[w1] == level && DEGPRUNE(indeg,outdeg,w1,n))
            retlev = level;
        else
#endif
        retlev = scan(level+1,ne,minarcs,maxarcs,sofar+1,oriented,group,n);
        DELELEMENT(x,k);
#ifdef DEGPRUNE
        --outdeg[w1]; --indeg[w0];
#endif
        if (retlev < level) return retlev;
    }

    if (!oriented && sofar + 2 + ne - level - 1 <= maxarcs)
    {
        k = 2*level;
        ADDELEMENT(x,k);
        ADDELEMENT(x,k+1);
        ix[sofar] = k;
        ix[sofar+1] = k+1;
#ifdef DEGPRUNE
        ++indeg[w0]; ++indeg[w1];
        ++outdeg[w0]; ++outdeg[w1];
        if (lastlev[w0] == level && DEGPRUNE(indeg,outdeg,w0,n)
         || lastlev[w1] == level && DEGPRUNE(indeg,outdeg,w1,n))
            retlev = level;
        else
#endif
        retlev = scan(level+1,ne,minarcs,maxarcs,sofar+2,oriented,group,n);
        DELELEMENT(x,k+1);
        DELELEMENT(x,k);
#ifdef DEGPRUNE
        --indeg[w0]; --indeg[w1];
        --outdeg[w0]; --outdeg[w1];
#endif
        if (retlev < level) return retlev;
    }

    return level-1;
}

/**************************************************************************/

static void
direct(graph *g, int nfixed, long minarcs, long maxarcs,
       boolean oriented, int m, int n)
{
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[100];
    grouprec *group;
    long ne;
    int i,j,k,j0,j1,deg;
    int isol0,isol1;  /* isolated vertices before and after nfixed */
    set *gi;
    int lab[MAXNV],ptn[MAXNV],orbits[MAXNV];
    set active[(MAXNV+WORDSIZE-1)/WORDSIZE];

    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    j0 = -1;  /* last vertex with degree 0 */
    j1 = n;   /* first vertex with degree > 0 */
    isol0 = isol1 = 0;
 
    ne = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
    {
        deg = 0;
        for (j = 0; j < m; ++j) deg += POPCOUNT(gi[j]);
        if (deg == 0)
        {
            lab[++j0] = i;
            if (i < nfixed) ++isol0; else ++isol1;
        }
        else
            lab[--j1] = i;
        ne += deg;
    }
    ne /= 2;
    ntisol = (isol0 >= 2 || isol1 >= 2);

    me = (2*ne + WORDSIZE - 1) / WORDSIZE;
    if (me == 0) me = 1;
    EMPTYSET(x,me);

    splitcount = splitres;
    if (splitmod == 1)
        splitlevel = -1;
    else
	splitlevel = (ne <= 34 ? ne/2 : 17);

#ifdef DEGPRUNE
    for (i = 0; i < n; ++i) indeg[i] = outdeg[i] = 0;
#endif
#ifndef SPLITTEST
    if (ne == 0 && minarcs <= 0 && (!Vswitch || ntisol) && splitres == 0)
    {
#ifdef DEGPRUNE
	for (i = 0; i < n; ++i)
	    if (DEGPRUNE(indeg,outdeg,i,n)) break;
	if (i == n)
#endif
            trythisone(NULL,0,n);    // only in case 0
        return;
    }
#endif

    if (oriented)
    {
        if (maxarcs < ne || minarcs > ne) return;
    }
    else
    {
        if (maxarcs < ne || minarcs > 2*ne) return;
    }

    if (n > MAXNV || ne > MAXNE)
    {
        fprintf(stderr,">E directg: MAXNV or MAXNE exceeded\n");
        exit(1);
    }

    for (i = 0; i < n; ++i) ptn[i] = 1;
    ptn[n-1] = 0;
    EMPTYSET(active,m);
    ADDELEMENT(active,0);

    for (i = 0; i <= j0; ++i)
    {
        if (i < n-1) ADDELEMENT(active,i+1);
        ptn[i] = 0;
    }

    for (i = j0+1; i < n; ++i)
    if (lab[i] < nfixed) break;

    if (i != j0+1 && i != n)
    {
        ptn[i-1] = 0;
        ADDELEMENT(active,i);
    }

    options.defaultptn = FALSE;
    options.userautomproc = groupautomproc;
    options.userlevelproc = grouplevelproc;

    nauty(g,lab,ptn,active,orbits,&options,&stats,workspace,100,m,n,NULL);

    if (stats.grpsize2 == 0)
        groupsize = stats.grpsize1 + 0.1;
    else
        groupsize = 0;

    if (Vswitch && groupsize == 1 && !ntisol) return;

    group = groupptr(FALSE);
    makecosetreps(group);

#ifdef DEGPRUNE
    for (i = 0; i <= j0; ++i)
        if (DEGPRUNE(indeg,outdeg,lab[i],n)) return;
#endif

    k = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
    {
        for (j = i; (j = nextelement(gi,m,j)) >= 0; )
        {
            v0[k] = i;
            v1[k] = j;
#ifdef DEGPRUNE
            lastlev[i] = lastlev[j] = k;
#endif
            edgeno[i][j] = 2*k;
            edgeno[j][i] = 2*k+1;
            ++k;
        }
    }

    lastrejok = FALSE;

    scan(0,ne,minarcs,maxarcs,0,oriented,group,n);
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    graph *g;
    int m,n,codetype;
    int argnum,j,nfixed;
    char *arg,sw;
    boolean badargs;
    boolean Tswitch,uswitch,eswitch,qswitch,oswitch,fswitch,sswitch;
    long minarcs,maxarcs;
    double t;
    char *infilename,*outfilename;
    FILE *infile;
    char msg[201];
    int msglen;
    long lsplitres,lsplitmod;

    HELP; PUTVERSION;

/* Put 2 in the following to catch incorrect uses of naut*1. */
    nauty_check(WORDSIZE,2,2*WORDSIZE,NAUTYVERSIONID);

    Tswitch = Gswitch = fswitch = Vswitch = FALSE;
    uswitch = eswitch = oswitch = qswitch = FALSE;
    sswitch = digraph6 = FALSE;
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
                     SWBOOLEAN('o',oswitch)
                else SWBOOLEAN('q',qswitch)
                else SWBOOLEAN('u',uswitch)
                else SWBOOLEAN('T',Tswitch)
                else SWBOOLEAN('G',Gswitch)
                else SWBOOLEAN('V',Vswitch)
                else SWINT('f',fswitch,nfixed,"directg -f")
		else SWRANGE('s',"/",sswitch,lsplitres,lsplitmod,"directg -s")
                else SWRANGE('e',":-",eswitch,minarcs,maxarcs,"directg -e")
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

    if ((Gswitch!=0) + (Tswitch!=0) + (uswitch!=0) >= 2)
        gt_abort(">E directg: -G, -T and -u are incompatible\n");

    digraph6 = !(Gswitch || Tswitch || uswitch);

    if (!eswitch)
    {
        minarcs = 0;
        maxarcs = NOLIMIT;
    }
    if (sswitch)
    {
	splitmod = lsplitmod;
	splitres = lsplitres;
    }
    else
    {
	splitmod = 1;
	splitres = 0;
    }
    splitcases = 0;
    if (sswitch && (splitmod < 1 || splitres < 0 || splitres >= splitmod))
    {
        fprintf(stderr,">E -s#/# requires 0 <= res < mod\n");
        exit(1);
    }

    if (!fswitch) nfixed = 0;

    if (!qswitch)
    {
        msg[0] = '\0';
        CATMSG0(">A directg");
        if (sswitch || eswitch || oswitch || uswitch || fswitch) CATMSG0(" -");
        if (oswitch) CATMSG0("o");
        if (uswitch) CATMSG0("u");
        if (Tswitch) CATMSG0("T");
        if (Gswitch) CATMSG0("G");
        if (Vswitch) CATMSG0("V");
        if (fswitch) CATMSG1("f%d",nfixed);
	if (sswitch) CATMSG2("s%d/%d",splitres,splitmod);
        if (eswitch) CATMSG2("e%ld:%ld",minarcs,maxarcs);
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

    dg_nin = 0;
    dg_ngen = 0;
    dg_nout = 0;
    dg_skipped = 0;

    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readg(infile,NULL,0,&m,&n)) == NULL) break;
        ++dg_nin;
#ifdef PROCESS
	if (m > 1)
	{
	    fprintf(stderr,">E m=1 is needed if PROCESS is defined\n");
	    exit(1);
	}
#endif
#ifdef INPUTGRAPH
	if (!INPUTGRAPH(g,m,n))
	{
	   ++dg_skipped;
	   FREES(g);
	   continue;
	}
#endif

        direct(g,nfixed,minarcs,maxarcs,oswitch,m,n);
        if (!uswitch && ferror(outfile))
	    gt_abort(">E directg output error\n");
        FREES(g);
    }
    t = CPUTIME - t;

#ifdef SUMMARY
    SUMMARY();
#endif

#ifdef SPLITTEST
    fprintf(stderr,">S %lu cases at level %d\n",splitcases,splitlevel);
#else
    if (!qswitch)
    {
        fprintf(stderr,">Z " COUNTER_FMT " graphs read from %s",
                       dg_nin,infilename); 
#ifdef INPUTGRAPH
	fprintf(stderr,"; " COUNTER_FMT " skipped",dg_skipped);
#endif
        fprintf(stderr,"; " COUNTER_FMT ,dg_nout);
        if (!uswitch)
            fprintf(stderr," digraphs written to %s",outfilename);
        else
            fprintf(stderr," digraphs generated");
        fprintf(stderr,"; %.2f sec\n",t);
    }
#endif

#ifdef GROUPTEST
    fprintf(stderr,"Group test = %lld\n",totallab);
#endif

    exit(0);
}
