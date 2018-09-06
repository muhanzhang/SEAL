/*****************************************************************************
*                                                                            *
* This is the main file for dreadnaut() version 2.6, which is a test-bed     *
*   for nauty() version 2.6.                                                 *
*                                                                            *
*   Subject to the copyright notice in the file COPYRIGHT.                   *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       10-Nov-87 : final changes for version 1.2                            *
*        5-Dec-87 - replaced all uses of fscanf() by appropriate uses        *
*                   of the new procedures readinteger() and readstring()     *
*                 - changed the '<' command slightly.  If a file of the      *
*                   given name cannot be openned, an attempt is made to      *
*                   open a file with the same name extended by DEFEXT.       *
*                 - improved error processing for 'n' command.               *
*       28-Sep-88 : changes for version 1.4 :                                *
*                 - replaced incorrect %d by %ld in fprintf for ? command    *
*       23-Mar-89 : changes for version 1.5 :                                *
*                 - added optional startup message                           *
*                 - enabled use of refine1 in 'i' command                    *
*                 - implemented $$ command                                   *
*                 - replaced ALLOCS test by DYNALLOC test                    *
*                 - modified @ command and added # command                   *
*                 - declared local procedures static                         *
*       25-Mar-89 - implemented k command                                    *
*       27-Mar-89 - implemented * and I commands                             *
*       29-Mar-89 - implemented K command                                    *
*        2-Apr-89 - added reporting of vertex-invariant statistics           *
*        2-Apr-89 - added ## command                                         *
*        4-Apr-89 - added triples(), quadruples(), adjtriang()               *
*                 - updated error reporting for nauty()                      *
*        5-Apr-89 - removed flushline() from g and e commands                *
*                 - added T command                                          *
*        6-Apr-89 - added cellquads() and distances()                        *
*       26-Apr-89 - modified ? command, added & and && commands              *
*                 - added indsets(), cliques(), cellquins()                  *
*       18-Aug-89 - made g, lab, canong dynamically allocated always         *
*        2-Mar-90 - added celltrips(), cellcliq(), cellind()                 *
*       13-Mar-90 - changed canong and savedg in output to h and h'          *
*       19-Mar-90 - revised help() a little                                  *
*       19-Apr-90 : changes for version 1.6                                  *
*                 - rewrote "*" command to avoid bug in Pyramid C compiler   *
*       20-Apr-90 - rewrote above rewrite to avoid bug in SUN3 gcc           *
*       23-Apr-90 - undid above rewrite and fixed *my* bug <blush> by        *
*                   making NUMINVARS have type int.  Sorry, gcc.             *
*       10-Nov-90 - added calls to null routines (see comment on code)       *
*       27-Aug-92 : renamed to version 1.7, no changes to this file          *
*        5-Jun-93 : renamed to version 1.7+, no changes to this file         *
*       18-Aug-93 : renamed to version 1.8, no changes to this file          *
*       17-Sep-93 : changes for version 1.9 :                                *
*                 - added invariant adjacencies()                            *
*        7-Jun-96 : changes for version 2.0 :                                *
*                 - added invariants cellfano() and cellfano2()              *
*                 - made y=10 the default                                    *
*       11-Jul-96 - added dynamic allocation                                 *
*                 - rewrote h command and added H command                    *
*                 - implemented M and R commands                             *
*       15-Aug-96 - changed z command to use sethash()                       *
*       30-Aug-96 - no need to declare seed; already in naututil.h           *
*       12-Sep-96 - let i and I commands use userrefproc                     *
*        9-Dec-96 - made y=infinity the default                              *
*        6-Sep-97 - allocated arrays before accepting commands               *
*        7-Sep-97 - make g,canong,savedg 1-d arrays even statically          *
*       22-Sep-97 - undid error introduced on 7-Sep (worksize)               *
*        9-Jan-00 - used *_check() instead of *_null()                       *
*       12-Feb-00 - minor code formatting                                    *
*       17-Aug-00 - now use tc_level from DEFAULTOPTIONS                     *
*       16-Nov-00 - made changes listed in nauty.h                           *
*       22-Apr-01 - include nautyinv.h                                       *
*                 - improve worksize processing for MAXN=0                   *
*        5-May-01 - k=0 1 automatic for *, also K=3 or K=0                   *
*        2-Jun-01 - added __ command for digraph converse                    *
*       18-Oct-01 - moved WORKSIZE to here                                   *
*       21-Nov-01 - use NAUTYREQUIRED in *_check() calls                     *
*        1-Sep-02 - Undid the previous change                                *
*       17-Nov-03 - Changed INFINITY to NAUTY_INFINITY                       *
*       15-Nov-04 - Completed all prototypes                                 *
*       23-Nov-06 - no usertcellproc() any more in version 2.4               *
*       10-Nov-09 - removed types shortish and permutation                   *
*       17-Nov-09 - added sparsegraphs, schreier, A and G commands           *
*       19-Nov-09 - added F command, stub for Traces                         *
*        1-Dec-09 - added traces refinement, M also applies to i             *
*       16-Dec-09 - added sr# command for random regular graphs              *
*       19-Dec-09 - added s# command for the sparse case                     *
*                 - w command is in units of 2*m now                         *
*       19-May-10 - Incorporate traces canonical labelling                   *
*        7-Jun-10 - implement %, _ and __ commands for sparse format         *
*        8-Jun-10 - add O and P commands, and sparse && command              *
*       11-Jun-10 - revise command-line parameters                           *
*       14-Jun-10 - digraphs and invariants (top level only) in Traces       *
*       26-Oct-10 - fix u command for sparse graphs; default G=10            *
*       14-Mar-11 - store partition with savedg                              *
*       21-Jul-11 - extend M command                                         *
*       24-Oct-11 - add S and OO commands                                    *
*       15-Jan-12 - use putorbitsplus() if USE_ANSICONTROLS                  *
*       20-Sep-12 - the first argument of ungetc is int, not char            *
*       18-Jan-13 - add code for ^C catching in nauty                        *
*                 - add usercanonproc sample                                 *
*                 - ->> means flush output file                              *
*       14-Nov-14 - fix numcells calculation in 'OO' command                 *
*       16-Mar-15 - add B command                                            *
*       16-Dec-15 - add r& command                                           *
*       22-Jan-16 - commands with short arguments must be all on one line    *
*                 - most errors cause rest of input line to be skipped       *
*       19-Feb-16 - make R command induce a partition if one is defined      *
*                                                                            *
*****************************************************************************/

#include "gtools.h"    /* which includes nauty.h, which includes stdio.h */
#include "nautinv.h"  
#include "schreier.h"
#include "traces.h"

#define USAGE "dreadnaut [-o options]"

#define HELPTEXT \
" Enter nauty+traces test program.\n\
\n\
  -o options  - set initial options.  The parameter value is a string of\n\
                dreadnaut commands from the following set:\n\
                a,c,d,m,p,l,G,P,w,y,$,A,V,M\n\
                The effect is the same as if these commands are entered\n\
                at the beginning of the standard input.\n\
  For help within dreadnaut, use the h command.\n"

#define PM(x) ((x) ? '+' : '-')
#define SS(n,sing,plur)  (n),((n)==1?(sing):(plur))
#define WORKSIZE 60
#define FLUSHANDPROMPT do { flushline(INFILE); if (prompt) fprintf(PROMPTFILE,"> "); } while (0)

#define SORT_OF_SORT 2
#define SORT_NAME sort2ints
#define SORT_TYPE1 int
#define SORT_TYPE2 int
#include "sorttemplates.c"   /* define sort2ints(a,b,n) */

#define INFILE fileptr[curfile]
#define SCHREIER_DEFAULT 10

static long seed;

#if !MAXN
DYNALLSTAT(graph,g,g_sz);
DYNALLSTAT(graph,canong,canong_sz);
DYNALLSTAT(graph,savedg,savedg_sz);
DYNALLSTAT(setword,workspace,workspace_sz);
DYNALLSTAT(int,lab,lab_sz);
DYNALLSTAT(int,ptn,ptn_sz);
DYNALLSTAT(int,orbits,orbits_sz);
DYNALLSTAT(int,templab,templab_sz);
DYNALLSTAT(int,tempptn,tempptn_sz);
DYNALLSTAT(int,perm,perm_sz);
DYNALLSTAT(int,savedlab,savedlab_sz);
DYNALLSTAT(int,savedptn,savedptn_sz);
DYNALLSTAT(set,tempactive,tempactive_sz);
DYNALLSTAT(set,active,active_sz);
#else
static graph g[MAXM*1L*MAXN];
static graph canong[MAXM*1L*MAXN];
static graph savedg[MAXM*1L*MAXN];
static setword workspace[MAXM*2L*WORKSIZE];
static int lab[MAXN];
static int ptn[MAXN];
static int orbits[MAXN];
static int savedlab[MAXN],savedptn[MAXN];
static int perm[MAXN];
static int templab[MAXN];
static int tempptn[MAXN];
static int tempactive[MAXM];
static set active[MAXM];
#endif

static sparsegraph g_sg;
static sparsegraph canong_sg;
static sparsegraph savedg_sg;

static DEFAULTOPTIONS_GRAPH(options);
static DEFAULTOPTIONS_SPARSEGRAPH(options_sg);
static statsblk stats;
static int curfile;
static FILE *fileptr[MAXIFILES];
static FILE *outfile;
static char def_ext[] = DEFEXT;
static boolean firstpath;       /* used in usernode() */

DEFAULTOPTIONS_TRACES(traces_opts);
static TracesStats traces_stats;

#define TMP

#define DENSE_MODE  0
#define SPARSE_MODE 1
#define TRACES_MODE 2
#define SPARSEREP(mode) ((mode)==1||(mode)==2)
#define NOSPARSEYET(c) else if (SPARSEREP(mode)) { fprintf(ERRFILE,\
              "command %s is not implemented in the sparse case\n",c); }
#define NODENSEYET else if (!SPARSEREP(mode)) { fprintf(ERRFILE,\
              "command %c is not implemented in the dense case\n",c); }
#define NOTRACESYET if (mode==TRACES_MODE) {  fprintf(ERRFILE,\
              "command %c is not implemented for Traces\n",c); }

static int mode;

#define U_NODE  1               /* masks for u values */
#define U_AUTOM 2
#define U_LEVEL 4
#define U_TCELL 8     /* At version 2.4, usertcellproc() is gone */
#define U_REF  16
#define U_CANON 32

#ifndef  NODEPROC
#define NODEPROC usernode
#else
extern void NODEPROC(graph*,int*,int*,int,int,int,int,int,int);
#endif

#ifndef  AUTOMPROC
#define AUTOMPROC userautom
#else
extern void AUTOMPROC(int,int*,int*,int,int,int);
#endif

#ifndef  LEVELPROC
#define LEVELPROC userlevel
#else
extern void LEVELPROC(int*,int*,int,int*,statsblk*,int,int,int,int,int,int);
#endif

#ifndef  REFPROC
#define REFPROC NULL
#else
extern void REFPROC(graph*,int*,int*,int,int*,int*,set*,int*,int,int);
#endif

#ifndef  CANONPROC
#define CANONPROC usercanon
#else
extern int CANONPROC(graph*,int*,graph*,int,int,int,int);
#endif

#ifndef  INVARPROC
#define INVARPROC NULL
#define INVARPROCNAME "none"
#else
extern void INVARPROC(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
#define INVARPROCNAME "user-defined"
#endif

#ifndef  INVARPROC_SG
#define INVARPROC_SG NULL
#define INVARPROCNAME_SG "none"
#else
extern void INVARPROC_SG(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
#define INVARPROCNAME_SG "user-defined"
#endif

static struct invarrec
{
    void (*entrypoint)(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
    char *name;
    void (*entrypoint_sg)(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
    char *name_sg;
} invarproc[]
    = {{INVARPROC, INVARPROCNAME, INVARPROC_SG, INVARPROCNAME_SG},
       {NULL,        "none",        NULL,           "none"},
       {twopaths,    "twopaths",    NULL,           "unavailable"},
       {adjtriang,   "adjtriang",   NULL,           "unavailable"},
       {triples,     "triples",     NULL,           "unavailable"},
       {quadruples,  "quadruples",  NULL,           "unavailable"},
       {celltrips,   "celltrips",   NULL,           "unavailable"},
       {cellquads,   "cellquads",   NULL,           "unavailable"},
       {cellquins,   "cellquins",   NULL,           "unavailable"},
       {distances,   "distances",   distances_sg,   "distances_sg"},
       {indsets,     "indsets",     NULL,           "unavailable"},
       {cliques,     "cliques",     NULL,           "unavailable"},
       {cellcliq,    "cellcliq",    NULL,           "unavailable"},
       {cellind,     "cellind",     NULL,           "unavailable"},
       {adjacencies, "adjacencies", adjacencies_sg, "adjacencies_sg"},
       {cellfano,    "cellfano",    NULL,           "unavailable"},
       {cellfano2,   "cellfano2",   NULL,           "unavailable"},
       {refinvar,    "refinvar",    NULL,           "unavailable"}
      };
#define NUMINVARS ((int)(sizeof(invarproc)/sizeof(struct invarrec)))

static void help(FILE*, int);
static void userautom(int,int*,int*,int,int,int);
static void usernode(graph*,int*,int*,int,int,int,int,int,int);
static void userlevel(int*,int*,int,int*,statsblk*,int,int,int,int,int,int);
static int usercanon(graph*,int*,graph*,int,int,int,int);

static boolean options_writeautoms,options_writemarkers,
            options_digraph,options_getcanon,options_linelength;
static int options_invarproc,options_mininvarlevel,options_maxinvarlevel,
	    options_invararg,options_tc_level,options_cartesian;
static int options_schreier,options_keepgroup,options_verbosity,
	   options_strategy;

#if USE_ANSICONTROLS && !DREADTEST
#define PUTORBITS putorbitsplus
#else
#define PUTORBITS putorbits
#endif

#ifdef  EXTRADECLS
EXTRADECLS
#endif

#if !HAVE_SIGACTION
#undef ALLOW_INTERRUPT
#define ALLOW_INTERRUPT 0
#endif

#if ALLOW_INTERRUPT
/*****************************************************************************
*                                                                            *
*  Routines for catching SIGINT                                              *
*                                                                            *
*****************************************************************************/

void
sigintcatcher(int sig)
/* This is the routine called on SIGINT receipt. */
{
    struct sigaction ss;

    nauty_kill_request = 1;
    ss.sa_handler = SIG_DFL;
    sigemptyset(&ss.sa_mask);
    ss.sa_flags = 0;
    sigaction(SIGINT,&ss,0);
}

static void
setsigcatcher(void)
{
    struct sigaction ss;

    nauty_kill_request = 0;
    ss.sa_handler = sigintcatcher;
    sigemptyset(&ss.sa_mask);
    ss.sa_flags = 0;
    sigaction(SIGINT,&ss,0);
}

static void
unsetsigcatcher(void)
{
    struct sigaction ss;

    nauty_kill_request = 0;
    ss.sa_handler = SIG_DFL;
    sigemptyset(&ss.sa_mask);
    ss.sa_flags = 0;
    sigaction(SIGINT,&ss,0);
}
#else
static void
setsigcatcher(void)
{
}

static void
unsetsigcatcher(void)
{
}
#endif

/*****************************************************************************
*                                                                            *
*  This is a program which illustrates the use of nauty.                     *
*  Commands are read from stdin, and may be separated by white space,        *
*  commas or not separated.  Output is written to stdout.                    *
*  For a short description, see the nauty User's Guide.                      *
*                                                                            *
*****************************************************************************/

int
main(int argc, char *argv[])
{
    int m,n,newm,newn;
    boolean gvalid,ovalid,cvalid,pvalid,minus,prompt,doquot;
    boolean gvalid_sg,cvalid_sg;
    int i,j,k,worksize,numcells,savednc,refcode,umask,qinvar;
    int oldorg,oldmode;
    int maxsize,cell1,cell2;
    boolean ranreg,same;
    char *s1,*s2;
    int c,d;
    unsigned long uli;
    size_t sli;
    set *gp;
    double timebefore,timeafter,mintime;
    char filename[515];
    int sgn,sgorg,nperm;
    int multiplicity,actmult;
    long zseed;
    permnode *generators;
    char *ap,*parameters;
    boolean flushing;

    HELP; PUTVERSION;

    if (argc == 3 && strcmp(argv[1],"-o") == 0)
	parameters = argv[2];
    else if (argc != 1)
    {
	fprintf(ERRFILE,USAGE);
	exit(1);
    }
    else
	parameters = "";
 
    mode = DENSE_MODE;
    curfile = 0;
    fileptr[curfile] = stdin;
#ifdef DREADTEST
    prompt = FALSE;
#else
    prompt = DOPROMPT(INFILE);
#endif
    outfile = stdout;
    options_writeautoms = options_writemarkers = TRUE;
    options_digraph = FALSE;
    options_getcanon = options.getcanon;
    options_mininvarlevel = options.mininvarlevel;
    options_maxinvarlevel = options.maxinvarlevel;
    options_invararg = options.invararg;
    options_invarproc = 1; /* index into invarproc[] */
    options_tc_level = options.tc_level;
    options_cartesian = options.cartesian;
    options_linelength = options.linelength;
    options_schreier = SCHREIER_DEFAULT;
    options_keepgroup = FALSE;
    generators = NULL;
    options_verbosity = 1;
    options_strategy = 0;

    n = m = 1;
    worksize = WORKSIZE;

#if !MAXN
    n = WORDSIZE;
    DYNALLOC2(graph,g,g_sz,n,m,"dreadnaut");
    DYNALLOC1(int,lab,lab_sz,n,"dreadnaut");
    DYNALLOC1(int,ptn,ptn_sz,n,"dreadnaut");
    DYNALLOC1(int,orbits,orbits_sz,n,"dreadnaut");
    DYNALLOC1(int,perm,perm_sz,n,"dreadnaut");
    DYNALLOC1(set,active,active_sz,m,"dreadnaut");
    n = 1;
#endif

#ifdef DREADTEST
    seed = 1;
    ran_init(seed);
#else
#ifdef  INITSEED
    INITSEED;
    ran_init(seed);
#endif
#endif

    umask = 0;
    pvalid = FALSE;
    ovalid = FALSE;
    gvalid = gvalid_sg = FALSE;  /* at most one valid */
    cvalid = cvalid_sg = FALSE;  /* at most one valid */
    sgorg = labelorg = oldorg = 0;
    sgn = 0;
    multiplicity = 1;
    mintime = 0.0;
    flushing = FALSE;

#ifdef  INITIALIZE
    INITIALIZE;
#endif

    if (prompt)
    {
        fprintf(PROMPTFILE,"Dreadnaut version %s.\n",NAUTYVERSION);
        fprintf(PROMPTFILE,"> ");
    }

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);
    nautinv_check(WORDSIZE,1,1,NAUTYVERSIONID);
    nautil_check(WORDSIZE,1,1,NAUTYVERSIONID);
    naututil_check(WORDSIZE,1,1,NAUTYVERSIONID);
    nausparse_check(WORDSIZE,1,1,NAUTYVERSIONID);

    SG_INIT(g_sg);
    SG_INIT(canong_sg);
    SG_INIT(savedg_sg);

    minus = FALSE;
    for (ap = parameters; *ap != '\0'; )
    {
	c = *ap++;

	switch (c)
	{
	case ' ':
	case '\t':
	    break;

	case '-':
	    minus = TRUE;
	    break;

	case '+':
	    minus = FALSE;
	    break;

	case 'a': 
	    options_writeautoms = !minus;
            minus = FALSE;
	    break;

	case 'c': 
	    options_getcanon = !minus;
            minus = FALSE;
	    break;

	case 'm': 
	    options_writemarkers = !minus;
            minus = FALSE;
	    break;

	case 'p': 
	    options_cartesian = !minus;
            minus = FALSE;
	    break;

	case 'B': 
	    flushing = !minus;
            minus = FALSE;
	    break;

	case 'd': 
	    options_digraph = !minus;
            minus = FALSE;
	    break;

	case 'P': 
	    options_keepgroup = !minus;
            minus = FALSE;
	    break;

	case 'w':
	    while (*ap == '=' || *ap == ' ') ++ap;
            arg_int(&ap,&worksize,"w");
#if MAXN
            if (worksize > 2*MAXM*WORKSIZE)
            {
                fprintf(ERRFILE,
                   "too big - setting worksize = %d\n",WORKSIZE);
                worksize = WORKSIZE;
            }
#endif
            minus = FALSE;
            break;

	case 'V':
	    if (minus)
            {
                options_verbosity = 0;
                minus = FALSE;
            }
            else
            {
	        while (*ap == '=' || *ap == ' ') ++ap;
                arg_int(&ap,&i,"V");
                if (i < 0) fprintf(ERRFILE,"verbosity must be >= 0\n");
                else       options_verbosity = i;
            }
            break;


	case 'S':
	    if (minus)
            {
                options_strategy = 0;
                minus = FALSE;
            }
            else
            {
	        while (*ap == '=' || *ap == ' ') ++ap;
                arg_int(&ap,&i,"S");
/*
                if (i < 0) fprintf(ERRFILE,"strategy must be >= 0\n");
                else       options_strategy = i;
*/
		if (i != 0) fprintf(ERRFILE,
                      "Only strategy 0 is supported in this version\n");
            }
            break;

        case 'G': 
            if (minus)
            {
                options_schreier = 0;
                minus = FALSE;
            }
            else
            {
	        while (*ap == '=' || *ap == ' ') ++ap;
                arg_int(&ap,&i,"G");
                if (i < 0)
		    fprintf(ERRFILE,"schreierfails must be >= 0\n");
                else
                {
                    options_schreier = i;
                    if (i > 0) schreier_fails(i);
                }
            }
            break;

        case 'y': 
	    while (*ap == '=' || *ap == ' ') ++ap;
            arg_int(&ap,&options_tc_level,"y");
            minus = FALSE;
            break;

        case '$': 
	    while (*ap == '=' || *ap == ' ') ++ap;
            arg_int(&ap,&labelorg,"$");
            minus = FALSE;
            break;

        case 'M': 
	    if (minus)
	    {
	        multiplicity = 1;
		mintime = 0.0;
                minus = FALSE;
	    }
	    else
	    {
		actmult = 0;
	        while (*ap == '=' || *ap == ' ') ++ap;
                arg_int(&ap,&multiplicity,"M");
		if (*ap == '/')
		{
		    ++ap;
		    arg_int(&ap,&actmult,"M/");
		}
	        if (multiplicity < 0) multiplicity = 0;
	        if (actmult < 0) actmult = 0;
		if (multiplicity == 0 && actmult == 0) multiplicity = 1;
		mintime = (double)actmult;
            }
            break;

        case 'l': 
	    while (*ap == '=' || *ap == ' ') ++ap;
            arg_int(&ap,&options_linelength,"l");
            minus = FALSE;
            break;

        case 'A': 
	    d = *ap++;
            if (d == 'n' || d == 'N' || d == 'd' || d == 'D')
		mode = DENSE_MODE;
            else if (d == 's' || d == 'S') mode = SPARSE_MODE;
            else if (d == 't' || d == 'T') mode = TRACES_MODE;
            else
            {
                fprintf(ERRFILE,"Mode %c is unknown\n",(d?d:'0'));
                break;
            }
	    minus = FALSE;
            break;

	default:
	    fprintf(ERRFILE,"Illegal initialization command %c\n",c);
	    exit(1);
	}
    }

    minus = FALSE;
    while (curfile >= 0)
    {
        if ((c = getc(INFILE)) == EOF || c == '\004')
        {
            fclose(INFILE);
            --curfile;
            if (curfile >= 0) prompt = DOPROMPT(INFILE);
        }
        else switch (c)
        {
        case '\n':  /* possibly issue prompt */
            if (prompt) fprintf(PROMPTFILE,"> ");
            minus = FALSE;
            break;

        case ' ':   /* do nothing */
        case '\t':
#ifndef  NLMAP
        case '\r':
#endif
        case '\f':
            break;

        case '-':   /* remember this for next time */
            minus = TRUE;
            break;

        case '+':   /* forget - */
        case ',':
        case ';':
            minus = FALSE;
            break;

        case '<':   /* new input file */
            minus = FALSE;
            if (curfile == MAXIFILES - 1)
	    {
                fprintf(ERRFILE,"exceeded maximum input nesting of %d\n",
                        MAXIFILES);
		FLUSHANDPROMPT;
		break;
	    }
            if (!readstring(INFILE,filename,513))
            {
                fprintf(ERRFILE,
                        "missing file name on '<' command : ignored\n");
                break;
            }
            if ((fileptr[curfile+1] = fopen(filename,"r")) == NULL)
            {
                for (s1 = filename; *s1 != '\0'; ++s1) {}
                for (s2 = def_ext; (*s1 = *s2) != '\0'; ++s1, ++s2) {}
                fileptr[curfile+1] = fopen(filename,"r");
            }
            if (fileptr[curfile+1] != NULL)
            {
                ++curfile;
                prompt = DOPROMPT(INFILE);
                if (prompt)
                    fprintf(PROMPTFILE,"> ");
            }
            else
	    {
                fprintf(ERRFILE,"can't open input file\n");
		FLUSHANDPROMPT;
	    }
            break;

        case '>':   /* new output file, or flush output file */
            if ((d = getc(INFILE)) != '>') ungetc(d,INFILE);
            if (minus)
            {
                minus = FALSE;
		if (d == '>')
		    fflush(outfile);
                else if (outfile != stdout)
                {
                    fclose(outfile);
                    outfile = stdout;
                }
            }
            else
            {
                if (!readstring(INFILE,filename,513))
                {
                    fprintf(ERRFILE,
                        "improper file name, reverting to stdout\n");
                    outfile = stdout;
		    FLUSHANDPROMPT;
                    break;
                }
                OPENOUT(outfile,filename,d=='>');
                if (outfile == NULL)
                {
                    fprintf(ERRFILE,
                        "can't open output file, reverting to stdout\n");
                    outfile = stdout;
		    FLUSHANDPROMPT;
                }
            }
            break;

	case 'B': 
	    flushing = !minus;
            minus = FALSE;
	    break;

        case '!':   /* ignore rest of line */
            do
                c = getc(INFILE);
            while (c != '\n' && c != EOF);
            if (c == '\n') ungetc('\n',INFILE);
            break;

        case 'n':   /* read n value */
            minus = FALSE;
            i = getint_sl(INFILE);
            if (i <= 0 || (MAXN && i > MAXN)
                       || (!MAXN && i > NAUTY_INFINITY-2))
	    {
                fprintf(ERRFILE,
                     " n can't be less than 1 or more than %d\n",
                       MAXN ? MAXN : NAUTY_INFINITY-2);
		FLUSHANDPROMPT;
	    }
            else
            {
                gvalid = FALSE;
                cvalid = FALSE;
                gvalid_sg = FALSE;
                cvalid_sg = FALSE;
                pvalid = FALSE;
                ovalid = FALSE;
                n = i;
                m = SETWORDSNEEDED(n);
                freeschreier(NULL,&generators); 
#if !MAXN
                DYNALLOC1(int,lab,lab_sz,n,"dreadnaut");
                DYNALLOC1(int,ptn,ptn_sz,n,"dreadnaut");
                DYNALLOC1(int,orbits,orbits_sz,n,"dreadnaut");
                DYNALLOC1(int,perm,perm_sz,n,"dreadnaut");
                DYNALLOC1(set,active,active_sz,m,"dreadnaut");
#endif
            }
            break;

        case 'g':   /* read graph */
            minus = FALSE;
            if (SPARSEREP(mode))
            {
                readgraph_sg(INFILE,&g_sg,options_digraph,prompt,
                             options_linelength,n);
                gvalid_sg = TRUE;
                cvalid_sg = FALSE;
            }
            else
            {
#if !MAXN
                DYNALLOC2(graph,g,g_sz,n,m,"dreadnaut");
#endif
                readgraph(INFILE,g,options_digraph,prompt,FALSE,
                          options_linelength,m,n);
                gvalid = TRUE;
                cvalid = FALSE;
            }
            ovalid = FALSE;
            break;

        case 'e':   /* edit graph */
            minus = FALSE;
            if (SPARSEREP(mode))
            {
                fprintf(ERRFILE,"e command only works in dense mode\n");
		FLUSHANDPROMPT;
            }
            else
            {
                readgraph(INFILE,g,options_digraph,prompt,gvalid,
                          options_linelength,m,n);
                gvalid = TRUE;
                cvalid = FALSE;
                ovalid = FALSE;
            }
            break;

        case 'r':   /* relabel graph and current partition */
            minus = FALSE;
	    if ((d = getc(INFILE)) != '&') ungetc(d,INFILE);
            if (gvalid_sg)
            {
                if (d == '&')
		{
		    if (pvalid)
			relabel_sg(&g_sg,lab,lab,&canong_sg);
		}
		else
		{
                    readvperm(INFILE,perm,prompt,n,&nperm);
                    relabel_sg(&g_sg,(pvalid ? lab : NULL),perm,&canong_sg);
		}
                cvalid_sg = FALSE;
                ovalid = FALSE;
            }
            else if (gvalid)
            {
		if (d == '&')
		{
		    if (pvalid)
		    {
#if !MAXN
                        DYNALLOC2(graph,canong,canong_sz,n,m,"dreadnaut");
#endif
			relabel(g,lab,lab,canong,m,n);
		    }
		}
		else
		{
#if !MAXN
                    DYNALLOC2(graph,canong,canong_sz,n,m,"dreadnaut");
#endif
                    readvperm(INFILE,perm,prompt,n,&nperm);
                    relabel(g,(pvalid ? lab : NULL),perm,canong,m,n);
		}
                cvalid = FALSE;
                ovalid = FALSE;
            }
            else
	    {
                fprintf(ERRFILE,"g is not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case 'R':   /* form subgraph */
            if (gvalid)
            {
#if !MAXN
                DYNALLOC2(graph,canong,canong_sz,n,m,"dreadnaut");
#endif
                readvperm(INFILE,perm,prompt,n,&nperm);
                if ((minus && nperm == n) || (!minus && nperm == 0))
		{
                    fprintf(ERRFILE,"can't form null graph\n");
		    FLUSHANDPROMPT;
		}
                else if (minus)
                {
                    sublabel(g,perm+nperm,n-nperm,canong,m,n);
                    if (pvalid) numcells = subpartition(lab,ptn,n,perm+nperm,n-nperm);
                    n = n - nperm;
                }
                else
                {
                    sublabel(g,perm,nperm,canong,m,n);
                    if (pvalid) numcells = subpartition(lab,ptn,n,perm,nperm);
                    n = nperm;
                }
                cvalid = FALSE;
                ovalid = FALSE;
                m = SETWORDSNEEDED(n);
            }
            else if (gvalid_sg)
            {
                readvperm(INFILE,perm,prompt,n,&nperm);
                if ((minus && nperm == n) || (!minus && nperm == 0))
		{
                    fprintf(ERRFILE,"can't form null graph\n");
		    FLUSHANDPROMPT;
		}
                else if (minus)
                {
                    sublabel_sg(&g_sg,perm+nperm,n-nperm,&canong_sg);
                    if (pvalid) numcells = subpartition(lab,ptn,n,perm+nperm,n-nperm);
                    n = n - nperm;
                }
                else
                {
                    sublabel_sg(&g_sg,perm,nperm,&canong_sg);
                    if (pvalid) numcells = subpartition(lab,ptn,n,perm,nperm);
                    n = nperm;
                }
                cvalid_sg = FALSE;
                ovalid = FALSE;
                m = SETWORDSNEEDED(n);
            }
            else
	    {
                fprintf(ERRFILE,"g is not defined\n");
		FLUSHANDPROMPT;
	    }
            minus = FALSE;
            break;

        case '_':   /* complement graph or converse digraph */
            minus = FALSE;
            if ((d = getc(INFILE)) != '_') ungetc(d,INFILE);

            if (!gvalid && !gvalid_sg)
	    {
                fprintf(ERRFILE,"g is not defined\n");
		FLUSHANDPROMPT;
	    }
            else if (gvalid)
            {
                if (d == '_') converse(g,m,n);
                else          complement(g,m,n);
                cvalid = FALSE;
            }
	    else
	    {
		if (d == '_')
		{
		    copy_sg(&g_sg,&canong_sg);
		    converse_sg(&canong_sg,&g_sg);
                    cvalid_sg = FALSE;
		}
	        else
		{
		    copy_sg(&g_sg,&canong_sg);
		    complement_sg(&canong_sg,&g_sg);
                    cvalid_sg = FALSE;
		}
	    }
            break;

        case '@':   /* copy canong into savedg */
            minus = FALSE;
            if (cvalid)
            {
#if !MAXN
                DYNALLOC2(graph,savedg,savedg_sz,n,m,"dreadnaut");
                DYNALLOC1(int,savedlab,savedlab_sz,n,"dreadnaut");
                DYNALLOC1(int,savedptn,savedptn_sz,n,"dreadnaut");
#endif
                sgn = n;
		memcpy(savedg,canong,m*(size_t)n*sizeof(setword));
                for (i = n; --i >= 0;)
		{
		    savedlab[i] = lab[i];
		    savedptn[i] = ptn[i];
		}
                sgorg = labelorg;
            }
            else if (cvalid_sg)
	    {
#if !MAXN
                DYNALLOC1(int,savedlab,savedlab_sz,n,"dreadnaut");
                DYNALLOC1(int,savedptn,savedptn_sz,n,"dreadnaut");
#endif
		sgn = n;
		copy_sg(&canong_sg,&savedg_sg);
                for (i = n; --i >= 0;)
		{
		    savedlab[i] = lab[i];
		    savedptn[i] = ptn[i];
		}
                sgorg = labelorg;
	    }
            else
	    {
                fprintf(ERRFILE,"h is not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case '#':   /* compare canong to savedg */
            if ((d = getc(INFILE)) != '#') ungetc(d,INFILE);

            if (cvalid || cvalid_sg)
            {
                if (sgn > 0)
                {
                    if (sgn != n)
                        fprintf(outfile,
                              "h and h' have different sizes.\n");
                    else
                    {
			if (cvalid)
		 	{
                            for (sli = 0; sli < m*(size_t)n; ++sli)
                                if (savedg[sli] != canong[sli]) break;
			    same = (sli == m*(size_t)n); 
			}
			else
			    same = aresame_sg(&canong_sg,&savedg_sg);

                        if (!same)
                            fprintf(outfile,"h and h' are different.\n");
                        else
                        {
			    for (i = 0; i < n; ++i)
				if ((ptn[i] == 0) != (savedptn[i] == 0))
				    break;
			    if (i < n)
                                fprintf(outfile,
                 "h and h' are identical but have incompatible colourings.\n");
			    else
                                fprintf(outfile,
                                    "h and h' are identical.\n");
                            if (d == '#')
                                putmapping(outfile,savedlab,sgorg,
                                       lab,labelorg,options_linelength,n);
                        }
                    }
                }
                else
		{
                    fprintf(ERRFILE,"h' is not defined\n");
		    FLUSHANDPROMPT;
		}
            }
            else
	    {
                fprintf(ERRFILE,"h is not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case 'j':   /* relabel graph randomly */
            minus = FALSE;
            if (gvalid)
            {
                ranperm(perm,n);
#if !MAXN
                DYNALLOC2(graph,canong,canong_sz,n,m,"dreadnaut");
#endif
                relabel(g,(pvalid?lab:NULL),perm,canong,m,n);
                cvalid = FALSE;
                ovalid = FALSE;
                freeschreier(NULL,&generators);
            }
            else if (gvalid_sg)
            {
                ranperm(perm,n);
                relabel_sg(&g_sg,(pvalid?lab:NULL),perm,&canong_sg);
                cvalid_sg = FALSE;
                ovalid = FALSE;
                freeschreier(NULL,&generators);
            }
            else
	    {
                fprintf(ERRFILE,"g is not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case 'v':   /* write vertex degrees */
            minus = FALSE;
            if ((d = getc(INFILE)) != 'v') ungetc(d,INFILE);

            if (gvalid)
	    {
                if (d == 'v') putdegseq(outfile,g,options_linelength,m,n);
                else          putdegs(outfile,g,options_linelength,m,n);
	    }
            else if (gvalid_sg)
	    {
                if (d == 'v') putdegseq_sg(outfile,&g_sg,options_linelength);
                else          putdegs_sg(outfile,&g_sg,options_linelength);
	    }
            else
	    {
                fprintf(ERRFILE,"g is not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case '%':   /* do Mathon doubling operation */
            minus = FALSE;
            if (gvalid || gvalid_sg)
            {
#if !MAXN
                if (2L * ((long)n + 1L) > NAUTY_INFINITY-2)
                {
                    fprintf(ERRFILE,
                         "n can't be more than %d\n",NAUTY_INFINITY-2);
                    break;
                }
#else
                if (2L * ((long)n + 1L) > MAXN)
                {
                    fprintf(ERRFILE,"n can't be more than %d\n",MAXN);
		    FLUSHANDPROMPT;
                    break;
                }
#endif
                newn = 2 * (n + 1);
                newm = SETWORDSNEEDED(newn);
#if !MAXN
                DYNALLOC1(int,lab,lab_sz,newn,"dreadnaut");
                DYNALLOC1(int,ptn,ptn_sz,newn,"dreadnaut");
                DYNALLOC1(int,orbits,orbits_sz,newn,"dreadnaut");
                DYNALLOC1(int,perm,perm_sz,newn,"dreadnaut");
                DYNALLOC1(set,active,active_sz,newm,"dreadnaut");
#endif
                ovalid = FALSE;
                pvalid = FALSE;
                freeschreier(NULL,&generators);
	    }
            else
	    {
                fprintf(ERRFILE,"g is not defined\n");
		FLUSHANDPROMPT;
	    }

	    if (gvalid)
	    {
#if !MAXN
                DYNALLOC2(graph,canong,canong_sz,n,m,"dreadnaut");
#endif
		memcpy(canong,g,m*(size_t)n*sizeof(setword));

#if !MAXN
                DYNALLOC2(graph,g,g_sz,newn,newm,"dreadnaut");
#endif
                mathon(canong,m,n,g,newm,newn);
                m = newm;
                n = newn;
                cvalid = FALSE;
            }
	    else if (gvalid_sg)
	    {
		copy_sg(&g_sg,&canong_sg);
		mathon_sg(&canong_sg,&g_sg);
                m = newm;
                n = newn;
                cvalid_sg = FALSE;
            }
            break;

        case 's':   /* generate random graph */
            minus = FALSE;
	    d = getc(INFILE);
	    if (d == 'r')
		ranreg = TRUE;
	    else
	    {
		ranreg = FALSE;
		if (d != EOF) ungetc(d,INFILE);
	    }

            i = getint_sl(INFILE);
	    if (ranreg)
	    {
		if (i < 0) i = 3;
		if (i > MAXREG)
		{
		    fprintf(ERRFILE,"sr is limited to degree %d\n",MAXREG);
		    FLUSHANDPROMPT;
		    break;
		}
		if (SPARSEREP(mode))
		{
		    ranreg_sg(&g_sg,i,n);
		    gvalid_sg = TRUE;
		    cvalid = FALSE;
                    ovalid = FALSE;
                    freeschreier(NULL,&generators);
                }
		NODENSEYET
	    }
	    else
	    {
                if (i <= 0) i = 2;
                if (!SPARSEREP(mode))
                {
#if !MAXN
                    DYNALLOC2(graph,g,g_sz,n,m,"dreadnaut");
#endif
                    rangraph(g,options_digraph,i,m,n);
                    gvalid = TRUE;
                    cvalid = FALSE;
                    ovalid = FALSE;
                    freeschreier(NULL,&generators);
                }
                else
		{
		    rangraph2_sg(&g_sg,options_digraph,1,i,n);
		    gvalid_sg = TRUE;
                    cvalid = FALSE;
                    ovalid = FALSE;
                    freeschreier(NULL,&generators);
                }
	    }
            break;

        case 'q':   /* quit */
            EXIT;
            break;

        case '"':   /* copy comment to output */
            minus = FALSE;
            copycomment(INFILE,outfile,'"');
            break;

        case 'I':   /* do refinement and invariants procedure */
	    minus = FALSE;
	    if (!gvalid && !gvalid_sg)
	    {
		fprintf(ERRFILE,"g is not valid\n");
		FLUSHANDPROMPT;
		break;
	    }
            if (!pvalid) unitptn(lab,ptn,&numcells,n);
            cellstarts(ptn,0,active,m,n);
#ifdef  CPUTIME
            timebefore = CPUTIME;
#endif
            if (gvalid)
	    {
                doref(g,lab,ptn,0,&numcells,&qinvar,perm,active,&refcode,
                    options.userrefproc ? options.userrefproc : 
                    (m == 1 ? refine1 : refine),
                    invarproc[options_invarproc].entrypoint,0,0,
                    options_invararg,options_digraph,m,n);
		if (numcells > 1) pvalid = TRUE;
	    }
            else if (gvalid_sg)
	    {
                doref((graph*)&g_sg,lab,ptn,0,&numcells,&qinvar,perm,
                    active,&refcode,
                    options_sg.userrefproc ? options_sg.userrefproc : 
                    refine_sg,
                    invarproc[options_invarproc].entrypoint_sg,0,0,
                    options_invararg,options_digraph,m,n);
		if (numcells > 1) pvalid = TRUE;
	    }
	    else
	    {
		fprintf(ERRFILE,"g is not valid\n");
		FLUSHANDPROMPT;
		break;
	    }
#ifdef  CPUTIME
            timeafter = CPUTIME;
#endif
            fprintf(outfile," %d cell%s; code = %x",
                    SS(numcells,"","s"),refcode);
	    if (mode == SPARSE_MODE)
	    {
		if (invarproc[options_invarproc].entrypoint_sg)
		    fprintf(outfile,
                       " (%s %s)",invarproc[options_invarproc].name_sg,
                       (qinvar == 2 ? "worked" : "failed"));
	    }
	    else if (mode == DENSE_MODE)
	    {
		if (invarproc[options_invarproc].entrypoint)
		    fprintf(outfile,
                       " (%s %s)",invarproc[options_invarproc].name,
                       (qinvar == 2 ? "worked" : "failed"));
	    }
#ifdef  CPUTIME
            fprintf(outfile,"; cpu time = %.2f seconds\n",
                    timeafter-timebefore);
#else
            fprintf(outfile,"\n");
#endif
            if (numcells > 1) pvalid = TRUE;
            break;

        case 'i':   /* do refinement */
	    minus = FALSE;
	    if (!gvalid && !gvalid_sg)
	    {
		fprintf(ERRFILE,"g is not valid\n");
		FLUSHANDPROMPT;
		break;
	    }
            if (!pvalid) unitptn(lab,ptn,&numcells,n);
            cellstarts(ptn,0,active,m,n);

            if (multiplicity != 0 || mintime != 0.0)
	    {
	        savednc = numcells;
#if !MAXN
    		DYNALLOC1(int,tempptn,tempptn_sz,n,"dreadnaut");
    		DYNALLOC1(int,templab,templab_sz,n,"dreadnaut");
    		DYNALLOC1(set,tempactive,tempactive_sz,m,"dreadnaut");
#endif
		memcpy(templab,lab,n*sizeof(int));
		memcpy(tempptn,ptn,n*sizeof(int));
		for (i = 0; i < m; ++i) tempactive[i] = active[i];
	    }

#ifdef  CPUTIME
            timebefore = CPUTIME;
#endif
	    actmult = 0;
	    for (;;)
	    {
		if (actmult > 0)
		{
		    memcpy(lab,templab,n*sizeof(int));
		    memcpy(ptn,tempptn,n*sizeof(int));
		    for (i = 0; i < m; ++i) active[i] = tempactive[i];
		    numcells = savednc;
		}

                if (options.userrefproc)
                    (*options.userrefproc)
                         (g,lab,ptn,0,&numcells,perm,active,&refcode,m,n);
	        else if (gvalid)
	        {
                    if (m == 1)
                        refine1(g,lab,ptn,0,&numcells,perm,active,&refcode,m,n);
                    else
                        refine(g,lab,ptn,0,&numcells,perm,active,&refcode,m,n);
	        }
	        else if (mode == SPARSE_MODE)
		    refine_sg((graph*)&g_sg,lab,ptn,0,&numcells,perm,active,
							         &refcode,m,n);
	        else  /* traces mode */
		    refine_tr(&g_sg,lab,ptn,&numcells,&refcode,&traces_opts);

		++actmult;
		if (multiplicity > 0 && actmult >= multiplicity) break;
#ifdef  CPUTIME
		if (mintime > 0.0 && (actmult < 20 || !(actmult&7)) 
		     	&& CPUTIME >= timebefore+mintime)
		    break;
#endif
	    }
#ifdef  CPUTIME
            timeafter = CPUTIME;
#endif
	    if (numcells > 1) pvalid = TRUE;
            fprintf(outfile," %d cell%s; code = %x",
                    SS(numcells,"","s"),refcode);
#ifdef  CPUTIME
            fprintf(outfile,"; cpu time = %.7f seconds\n",
			      (timeafter-timebefore)/actmult);
#else
	    fprintf(outfile,"\n");
#endif
            break;

        case 'x':   /* execute nauty */
            minus = FALSE;
            if (mode == TRACES_MODE)
	    {
                ovalid = FALSE;
                cvalid_sg = FALSE;
                if (!gvalid_sg)
                {
                    fprintf(ERRFILE,"g is not defined\n");
		    FLUSHANDPROMPT;
                    break;
                }

	        traces_opts.getcanon = options_getcanon;
	        traces_opts.writeautoms = options_writeautoms;
	        traces_opts.cartesian = options_cartesian;
	        traces_opts.linelength = options_linelength;
	        traces_opts.digraph = options_digraph;
	        traces_opts.outfile = outfile;
		traces_opts.verbosity = options_verbosity;
		traces_opts.strategy = options_strategy;
		if (options_keepgroup)
		    traces_opts.generators = &generators;
		else
		    traces_opts.generators = NULL;

#if !MAXN
		DYNALLOC1(int,tempptn,tempptn_sz,n,"dreadnaut"); 
#endif
		if (!pvalid) unitptn(lab,ptn,&numcells,n);
		memcpy(tempptn,ptn,n*sizeof(int));
		savednc = numcells;
		if (options_invarproc != 1 && options_maxinvarlevel > 0)
                {
                    if (options_maxinvarlevel > 1) fprintf(ERRFILE,
		    "Warning: Traces only uses invariants at the top level\n");
                    if (invarproc[options_invarproc].entrypoint_sg)
		    {
#ifdef  CPUTIME
                        timebefore = CPUTIME;
#endif
			cellstarts(tempptn,0,active,m,n);
			doref((graph*)&g_sg,lab,tempptn,0,&savednc,&qinvar,perm,
			    active,&refcode,
                    	    options_sg.userrefproc ? options_sg.userrefproc : 
                    	    refine_sg,
                    	    invarproc[options_invarproc].entrypoint_sg,0,0,
                    	    options_invararg,options_digraph,m,n);
                        fprintf(outfile,"Invariant %s %s; %d cell%s",
			    invarproc[options_invarproc].name_sg,
                            (qinvar == 2 ? "worked" : "failed"),
			    SS(savednc,"","s"));
#ifdef  CPUTIME
                        timeafter = CPUTIME;
			fprintf(outfile,"; cpu time = %.2f seconds\n",
                    	    timeafter-timebefore);
#else
		 	fprintf(outfile,"\n");
#endif
		    }
                }

#ifdef  CPUTIME
                timebefore = CPUTIME;
#endif
		actmult = 0;
                setsigcatcher(); 
                for (;;)
                {
                    traces_opts.defaultptn = !pvalid;
                    Traces(&g_sg,lab,tempptn,orbits,&traces_opts,&traces_stats,
		       &canong_sg);
		    if (traces_stats.errstatus) break;
                    traces_opts.writeautoms = FALSE;
		    traces_opts.verbosity = 0;
		    ++actmult;
		    if (multiplicity > 0 && actmult >= multiplicity) break;
#ifdef  CPUTIME
		    if (mintime > 0.0 && (actmult < 20 || !(actmult&7))
		       	   && CPUTIME >= timebefore+mintime)
			break;
#endif
                }
                unsetsigcatcher(); 
#ifdef  CPUTIME
                timeafter = CPUTIME;
#endif
            if (traces_stats.errstatus)
            {
                if (traces_stats.errstatus == NAUABORTED)
                    fprintf(ERRFILE,"Traces aborted\n");
                else if (traces_stats.errstatus == NAUKILLED)
                    fprintf(ERRFILE,"Traces interrupted\n");
                else
                    fprintf(ERRFILE,
                      "Traces returned error status %d [this can't happen]\n",
                      traces_stats.errstatus);
                cvalid = cvalid_sg = ovalid = FALSE;
            }
            else
            {
                fprintf(outfile,"%d orbit%s",
				SS(traces_stats.numorbits,"","s"));
		fprintf(outfile,"; grpsize=");
		writegroupsize(outfile,
			       traces_stats.grpsize1,traces_stats.grpsize2);
                fprintf(outfile,"; %d gen%s",
                     SS(traces_stats.numgenerators,"","s"));
                fprintf(outfile,
		     "; %lu node%s ", SS(traces_stats.numnodes,"","s"));
		if (traces_stats.interrupted)
		    fprintf(outfile,
			    "(%lu interrupted, ",traces_stats.interrupted);
		else
		    fprintf(outfile,"(");
                fprintf(outfile,"%lu peak); maxlev=%d\n",
		    traces_stats.peaknodes,traces_stats.treedepth);
                if (options_getcanon)
                    fprintf(outfile,
		            "canupdates=%d; ",traces_stats.canupdates);
#ifdef  CPUTIME
                fprintf(outfile,actmult == 1 ?
                              "cpu time = %.2f seconds\n" :
                              "cpu time = %.7f seconds\n",
			      (timeafter-timebefore)/actmult);
#else
		fprintf(outfile,"\n");
#endif
		if (options_getcanon) cvalid_sg = TRUE;
		ovalid = TRUE;
            } 
	    }
            else 
            {
                ovalid = FALSE;
                cvalid = cvalid_sg = FALSE;
                if (!gvalid && !gvalid_sg)
                {
                    fprintf(ERRFILE,"g is not defined\n");
		    FLUSHANDPROMPT;
                    break;
                }
		if (mode == DENSE_MODE)
		{
                    if (pvalid)
                    {
                        fprintf(outfile,"[fixing partition]\n");
                        options.defaultptn = FALSE;
                    }
                    else
                        options.defaultptn = TRUE;

                    options.outfile = outfile;
                    options.digraph = options_digraph;
                    options.cartesian = options_cartesian;
                    options.schreier = (options_schreier > 0);
                    options.getcanon = options_getcanon;
                    options.tc_level = options_tc_level;
                    options.linelength = options_linelength;
		    options.invarproc
				= invarproc[options_invarproc].entrypoint;
		    options.mininvarlevel = options_mininvarlevel;
		    if (options.invarproc)
		       options.maxinvarlevel = options_maxinvarlevel;
		    else
		       options.maxinvarlevel = 0;
		    options.invararg = options_invararg;
		    if (options_schreier > 0)
			schreier_fails(options_schreier);

                    if (umask & U_NODE)  options.usernodeproc = NODEPROC;
                    else                 options.usernodeproc = NULL;
                    if (umask & U_AUTOM) options.userautomproc = AUTOMPROC;
                    else                 options.userautomproc = NULL;
                    if (umask & U_LEVEL) options.userlevelproc = LEVELPROC;
                    else                 options.userlevelproc = NULL;
                    if (umask & U_REF)   options.userrefproc = REFPROC;
                    else                 options.userrefproc = NULL;
                    if (umask & U_CANON) options.usercanonproc = CANONPROC;
                    else                 options.usercanonproc = NULL;
#if !MAXN
                    if (options_getcanon)
                        DYNALLOC2(graph,canong,canong_sz,n,m,"dreadnaut");
                    DYNALLOC1(setword,workspace,workspace_sz,2*m*worksize,
								"dreadnaut");
#endif
                    firstpath = TRUE;
                    options.writeautoms = options_writeautoms;
                    options.writemarkers = options_writemarkers;
#ifdef  CPUTIME
                    timebefore = CPUTIME;
#endif
		    actmult = 0;
		    setsigcatcher();
                    for (;;)
                    {
                        nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,
                             2*m*worksize,m,n,canong);
			if (stats.errstatus) break;
                        options.writeautoms = FALSE;
                        options.writemarkers = FALSE;
			++actmult;
			if (multiplicity > 0 && actmult >= multiplicity)
			    break;
#ifdef  CPUTIME
			if (mintime > 0.0 && (actmult < 20 || !(actmult&7)) 
	  		        && CPUTIME >= timebefore+mintime)
			    break;
#endif
                    }
		    unsetsigcatcher();
#ifdef  CPUTIME
                    timeafter = CPUTIME;
#endif
		}
		if (mode == SPARSE_MODE)
		{
                    if (pvalid)
                    {
                        fprintf(outfile,"[fixing partition]\n");
                        options_sg.defaultptn = FALSE;
                    }
                    else
                        options_sg.defaultptn = TRUE;

                    options_sg.outfile = outfile;
                    options_sg.digraph = options_digraph;
                    options_sg.cartesian = options_cartesian;
                    options_sg.schreier = (options_schreier > 0);
                    options_sg.getcanon = options_getcanon;
                    options_sg.linelength = options_linelength;
		    options_sg.invarproc
                         = invarproc[options_invarproc].entrypoint_sg;
		    options_sg.mininvarlevel = options_mininvarlevel;
		    if (options_sg.invarproc)
		       options_sg.maxinvarlevel = options_maxinvarlevel;
		    else
		       options_sg.maxinvarlevel = 0;
		    options_sg.invararg = options_invararg;
                    options_sg.tc_level = options_tc_level;
		    if (options_schreier > 0)
			schreier_fails(options_schreier);

                    if (umask & U_NODE)  options_sg.usernodeproc = NODEPROC;
                    else                 options_sg.usernodeproc = NULL;
                    if (umask & U_AUTOM) options_sg.userautomproc = AUTOMPROC;
                    else                 options_sg.userautomproc = NULL;
                    if (umask & U_LEVEL) options_sg.userlevelproc = LEVELPROC;
                    else                 options_sg.userlevelproc = NULL;
                    if (umask & U_REF)   options_sg.userrefproc = REFPROC;
                    else                 options_sg.userrefproc = NULL;
                    if (umask & U_CANON) options_sg.usercanonproc = CANONPROC;
                    else                 options_sg.usercanonproc = NULL;
#if !MAXN
                    DYNALLOC1(setword,workspace,workspace_sz,2*m*worksize,
								"dreadnaut");
#endif

                    firstpath = TRUE;
                    options_sg.writeautoms = options_writeautoms;
                    options_sg.writemarkers = options_writemarkers;
#ifdef  CPUTIME
                    timebefore = CPUTIME;
#endif
		    actmult = 0;
		    setsigcatcher();
                    for (;;)
                    {
                        nauty((graph*)&g_sg,lab,ptn,NULL,orbits,&options_sg,
                         &stats,workspace,2*m*worksize,m,n,(graph*)&canong_sg);
			if (stats.errstatus) break;
                        options_sg.writeautoms = FALSE;
                        options_sg.writemarkers = FALSE;
			++actmult;
			if (multiplicity > 0 && actmult >= multiplicity)
			    break;
#ifdef  CPUTIME
			if (mintime > 0.0 && (actmult < 20 || !(actmult&7)) 
	  		        && CPUTIME >= timebefore+mintime)
			    break;
#endif
                    }
		    unsetsigcatcher();
#ifdef  CPUTIME
                    timeafter = CPUTIME;
#endif
		}

		if (stats.errstatus)
		{
		    if (stats.errstatus == NAUABORTED)
		        fprintf(ERRFILE,"nauty aborted\n");
		    else if (stats.errstatus == NAUKILLED)
		        fprintf(ERRFILE,"nauty interrupted\n");
                    else 
                        fprintf(ERRFILE,
                        "nauty returned error status %d [this can't happen]\n",
                           stats.errstatus);
		    cvalid = cvalid_sg = ovalid = FALSE;
		}
                else
                {
                    if (options_getcanon)
		    {
			if (mode == DENSE_MODE) cvalid = TRUE;
			else                    cvalid_sg = TRUE;
		    }
                    ovalid = TRUE;
                    fprintf(outfile,"%d orbit%s",SS(stats.numorbits,"","s"));
		    fprintf(outfile,"; grpsize=");
		    writegroupsize(outfile,stats.grpsize1,stats.grpsize2);
                    fprintf(outfile,"; %d gen%s",
                            SS(stats.numgenerators,"","s"));
                    fprintf(outfile,"; %lu node%s",SS(stats.numnodes,"","s"));
                    if (stats.numbadleaves)
                        fprintf(outfile," (%lu bad lea%s)",
                            SS(stats.numbadleaves,"f","ves"));
                    fprintf(outfile,"; maxlev=%d\n", stats.maxlevel);
                    /* fprintf(outfile,"tctotal=%lu",stats.tctotal); */
                    if (options_getcanon)
                        fprintf(outfile,"canupdates=%lu; ",stats.canupdates);
#ifdef  CPUTIME
                    fprintf(outfile,actmult == 1 ?
                              "cpu time = %.2f seconds\n" :
                              "cpu time = %.7f seconds\n",
                            (timeafter-timebefore)/actmult);
#else
                    fprintf(outfile,"\n");
#endif
		    if (mode == DENSE_MODE && options_maxinvarlevel != 0
		       && invarproc[options_invarproc].entrypoint)
                    {
                        fprintf(outfile,"invarproc \"%s\" succeeded %lu/%lu",
                            invarproc[options_invarproc].name,
			    stats.invsuccesses,stats.invapplics);
                        if (stats.invarsuclevel > 0)
                            fprintf(outfile," beginning at level %d.\n",
                                    stats.invarsuclevel);
                        else
                            fprintf(outfile,".\n");
                    }
		    if (mode == SPARSE_MODE && options_maxinvarlevel != 0
		       && invarproc[options_invarproc].entrypoint_sg)
                    {
                        fprintf(outfile,"invarproc \"%s\" succeeded %lu/%lu",
                            invarproc[options_invarproc].name_sg,
			    stats.invsuccesses,stats.invapplics);
                        if (stats.invarsuclevel > 0)
                            fprintf(outfile," beginning at level %d.\n",
                                    stats.invarsuclevel);
                        else
                            fprintf(outfile,".\n");
                    }
                }
	    }
            break;

        case 'A':   /* change mode, with possible conversion */
	    minus = FALSE;
	    oldmode = mode;
	    d = getc(INFILE);
	    if (d == 'n' || d == 'N' || d == 'd' || d == 'D') mode = DENSE_MODE;
	    else if (d == 's' || d == 'S') mode = SPARSE_MODE;
	    else if (d == 't' || d == 'T') mode = TRACES_MODE;
	    else
	    {
		fprintf(ERRFILE,"Mode %c is unknown\n",(d?d:'0'));
		FLUSHANDPROMPT;
		break;
	    }
	    if ((d = getc(INFILE)) != '+')
	    {
		ungetc(d,INFILE);
		gvalid = gvalid_sg = FALSE;
		pvalid = ovalid = FALSE;
	    }
	    else
	    {
		if (SPARSEREP(oldmode) && !SPARSEREP(mode) && gvalid_sg)
		{
#if !MAXN
		    DYNALLOC2(graph,g,g_sz,n,m,"dreadnaut");
#endif
		    sg_to_nauty(&g_sg,g,m,&m);
		    gvalid_sg = FALSE;
		    gvalid = TRUE;
		}
		if (!SPARSEREP(oldmode) && SPARSEREP(mode) && gvalid)
		{
		    nauty_to_sg(g,&g_sg,m,n);
		    gvalid = FALSE;
		    gvalid_sg = TRUE;
		}
	    }
	    cvalid = cvalid_sg = FALSE;
	    sgn = 0;    /* invalidate saved graph */
	    break;

        case 'f':   /* read initial partition */

            if (minus)
            {
                pvalid = FALSE;
                minus = FALSE;
            }
            else
            {
                readptn(INFILE,lab,ptn,&numcells,prompt,n);
                pvalid = TRUE;
	        freeschreier(NULL,&generators);
            }
            break;

	case 'F':   /* individualise one more vertex */
            if ((d = getc(INFILE)) != 'F') ungetc(d,INFILE);
	    minus = FALSE;
	    if (d != 'F')
	    {
	        i = getint_sl(INFILE);
	        i -= labelorg;
	        if (i < 0 || i >= n)
		{
		    fprintf(ERRFILE,"F argument must be 0..n-1\n");
		    FLUSHANDPROMPT;
		}
	        else
	        {
	            if (!pvalid) unitptn(lab,ptn,&numcells,n);
		    individualise(lab,ptn,0,i,&d,&numcells,n);
		    pvalid = TRUE;
	        }
	    }
	    else
	    {
		if (!gvalid && !gvalid_sg)
		{
		    fprintf(stderr,"g is not defined\n");
		    FLUSHANDPROMPT;
		    break;
		}
		if (!pvalid) unitptn(lab,ptn,&numcells,n);

		if (!SPARSEREP(mode))
		    i = targetcell(g,lab,ptn,0,1,options_digraph,-1,m,n);
		else if (mode == SPARSE_MODE)
		    i = targetcell_sg((graph*)&g_sg,lab,ptn,0,1,
						options_digraph,-1,m,n);
		else /* Traces mode */
	        {
		    maxsize = 0;
		    for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
                    {
                        for (cell2 = cell1; ptn[cell2] > 0; ++cell2) {}
			if (cell2-cell1+1 > maxsize)
			{
			    i = cell1;
			    maxsize = cell2 - cell1 + 1;
			}
		    }
		}
		
		if (ptn[i] > 0)
		{
		    ptn[i] = 0;
		    ++numcells;
		}
		pvalid = TRUE;
	    }
	    freeschreier(NULL,&generators);
	    break;

        case 't':   /* type graph */
            minus = FALSE;
            if (gvalid)
                putgraph(outfile,g,options_linelength,m,n);
            else if (gvalid_sg)
                putgraph_sg(outfile,&g_sg,options_linelength);
            else
	    {
                fprintf(ERRFILE,"g is not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case 'T':   /* type graph preceded by n, $ and g commands */
            minus = FALSE;
            if (gvalid)
            {
                fprintf(outfile,"n=%d $=%d g\n",n,labelorg);
                putgraph(outfile,g,options_linelength,m,n);
                fprintf(outfile,"$$\n");
            }
            else if (gvalid_sg)
            {
                fprintf(outfile,"n=%d $=%d g\n",n,labelorg);
                putgraph_sg(outfile,&g_sg,options_linelength);
                fprintf(outfile,"$$\n");
            }
            else
	    {
                fprintf(ERRFILE,"g is not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case 'u':   /* call user procs */
            if (minus)
            {
                umask = 0;
                minus = FALSE;
            }
            else
            {
                umask = getint_sl(INFILE);
                if (umask < 0) umask = ~U_TCELL;
            }
            if (umask & U_TCELL)
	    {
                fprintf(ERRFILE,"usertcellproc() is gone at version 2.4\n");
		umask &= ~U_TCELL;
	    }
            break;

        case 'o':   /* type orbits */
            minus = FALSE;
            if (ovalid)
                PUTORBITS(outfile,orbits,options_linelength,n);
            else
	    {
                fprintf(ERRFILE,"orbits are not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case 'O':   /* make orbits into a partition*/
            minus = FALSE;
            if ((d = getc(INFILE)) != 'O') ungetc(d,INFILE);

            if (ovalid && d != 'O')
	    {
#if !MAXN
    		DYNALLOC1(int,tempptn,tempptn_sz,n,"dreadnaut");
#endif
		for (i = n; --i >= 0;) tempptn[i] = 0;
                for (i = n; --i >= 0;)
	            if ((j = orbits[i]) < i)
	            {
		        tempptn[i] = tempptn[j];
			tempptn[j] = i;
		    }

                k = 0; 
		numcells = 0;
		for (i = 0; i < n; ++i)
		{
		    if (orbits[i] == i)
		    {
			j = i;
			do
			{
			    lab[k] = j;
			    ptn[k] = 1;
			    ++k;
			    j = tempptn[j];
			} while (j > 0);

			ptn[k-1] = 0;
			++numcells;
		    }
		}
		pvalid = TRUE;
	    }
	    else if (ovalid && d == 'O')
	    {
#if !MAXN
    		DYNALLOC1(int,tempptn,tempptn_sz,n,"dreadnaut");
    		DYNALLOC1(int,templab,templab_sz,n,"dreadnaut");
#endif
		for (i = 0; i < n; ++i) tempptn[i] = 0;
		for (i = 0; i < n; ++i) ++tempptn[orbits[i]];

		j = 0;
		for (i = 0; i < n; ++i)
		    if (tempptn[i] > 0)
		    {
			templab[j] = i;
			tempptn[j] = tempptn[i];
			++j;
		    }
		sort2ints(tempptn,templab,j);

		k = 0;
		for (i = 0; i < j; ++i)
		{ 
		    ptn[templab[i]] = k;
		    k += tempptn[i];
		}
		for (i = 0; i < n; ++i) lab[ptn[orbits[i]]++] = i;

		for (i = 0; i < n; ++i) ptn[i] = 1;
		k = 0;
                numcells = 0;
		for (i = 0; i < j; ++i)
		{ 
		    k += tempptn[i];
		    if (i == j-1 || tempptn[i] != tempptn[i+1])
		    {
			ptn[k-1] = 0;
			++numcells;
		    }
		}
		pvalid = TRUE;
	    }
            else
	    {
                fprintf(ERRFILE,"orbits are not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case 'b':   /* type canonlab and canong */
            minus = FALSE;
            if (cvalid)
                putcanon(outfile,lab,canong,options_linelength,m,n);
            else if (cvalid_sg)
	    {
		sortlists_sg(&canong_sg);
                putcanon_sg(outfile,lab,&canong_sg,options_linelength);
	    }
            else
	    {
                fprintf(ERRFILE,"h is not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case 'z':   /* type hashcode for canong */
            minus = FALSE;
            if (cvalid)
            {
                zseed = hashgraph(canong,m,n,2922320L);
                fprintf(outfile,"[N%07lx",zseed);
                
                zseed = hashgraph(canong,m,n,19883109L);
                fprintf(outfile," %07lx",zseed);
                
                zseed = hashgraph(canong,m,n,489317L); 
                fprintf(outfile," %07lx]\n",zseed);
            }
            else if (cvalid_sg)
            {
                zseed = hashgraph_sg(&canong_sg,2922320L);
                fprintf(outfile,"[%c%07lx",
		        mode==SPARSE_MODE?'S':'T',zseed);
                
                zseed = hashgraph_sg(&canong_sg,19883109L);
                fprintf(outfile," %07lx",zseed);
                
                zseed = hashgraph_sg(&canong_sg,489317L); 
                fprintf(outfile," %07lx]\n",zseed);
            }
            else
	    {
                fprintf(ERRFILE,"h is not defined\n");
		FLUSHANDPROMPT;
	    }
            break;

        case 'c':   /* set getcanon option */
            options_getcanon = !minus;
            minus = FALSE;
            break;

        case 'w':   /* read size of workspace */
            minus = FALSE;
            worksize = getint_sl(INFILE);
#if MAXN
            if (worksize > 2*MAXM*WORKSIZE)
            {
                fprintf(ERRFILE,
                   "too big - setting worksize = %d\n",WORKSIZE);
                worksize = WORKSIZE;
            }
#endif
            break;

        case 'l':   /* read linelength for output */
            options_linelength = getint_sl(INFILE);
            minus = FALSE;
            break;

        case 'y':   /* set tc_level field of options */
            options_tc_level = getint_sl(INFILE);
            minus = FALSE;
            break;

        case 'M':   /* set multiplicity */
	    if (minus)
	    {
	        multiplicity = 1;
		mintime = 0.0;
                minus = FALSE;
	    }
	    else
	    {
                multiplicity = getint_sl(INFILE);
                if (multiplicity < 0) multiplicity = 0;
		if ((d = getc(INFILE)) == '/')
		{
		    actmult = getint_sl(INFILE);
		    if (actmult < 0) actmult = 0;
		}
		else
		{
		    ungetc(d,INFILE);
		    actmult = 0;
		}
		if (multiplicity == 0 && actmult == 0) multiplicity = 1;
		mintime = (double)actmult;
	    }
            break;

        case 'k':   /* set invarlev fields of options */
            options_mininvarlevel = getint_sl(INFILE);
            options_maxinvarlevel = getint_sl(INFILE);
            minus = FALSE;
            break;

        case 'K':   /* set invararg field of options */
            options_invararg = getint_sl(INFILE);
            minus = FALSE;
            break;

        case '*':   /* set invarproc field of options */
            minus = FALSE;
            d = getint_sl(INFILE);
            if (d >= -1 && d <= NUMINVARS-2)
            {
                options_invarproc = d+1;
                options_mininvarlevel = 0;
                options_maxinvarlevel = 1;
                if (options_invarproc >= 10 && options_invarproc <= 13)
                    options.invararg = 3;
                else
                    options.invararg = 0;
            }
            else
	    {
                fprintf(ERRFILE,"no such vertex-invariant\n");
		FLUSHANDPROMPT;
	    }
            break;

        case 'a':   /* set writeautoms option */
            options_writeautoms = !minus;
            minus = FALSE;
            break;

        case 'm':   /* set writemarkers option */
            options_writemarkers = !minus;
            minus = FALSE;
            break;

        case 'V':   /* set verbosity for Traces */
	    if (minus)
	    {
		options_verbosity = 0;
		minus = FALSE;
	    }
	    else
	    {
                i = getint_sl(INFILE);
                if (i < 0)
		{
		    fprintf(ERRFILE,"verbosity must be >= 0\n");
		    FLUSHANDPROMPT;
		}
                else
                    options_verbosity = i;
	    }
            break;

        case 'S':   /* set strategy for Traces */
	    if (minus)
	    {
		options_strategy = 0;
		minus = FALSE;
	    }
	    else
	    {
                i = getint_sl(INFILE);
		if (i != 0) fprintf(ERRFILE,
                      "Only strategy 0 is supported in this version\n");
/*
                if (i < 0)
		{
		    fprintf(ERRFILE,"strategy must be >= 0\n");
		    FLUSHANDPROMPT;
		}
                else
                    options_strategy = i;
*/
	    }
            break;

        case 'G':   /* set schreier option */
	    if (minus)
	    {
		options_schreier = 0;
		minus = FALSE;
	    }
	    else
	    {
                i = getint_sl(INFILE);
                if (i < 0)
		{
		    fprintf(ERRFILE,"schreierfails must be >= 0\n");
		    FLUSHANDPROMPT;
		}
                else
		{
		    options_schreier = i;
		    if (i > 0) schreier_fails(i);
		}
	    }
            break;

        case 'p':   /* set cartesian option */
            options_cartesian = !minus;
            minus = FALSE;
            break;

        case 'd':   /* set digraph option */
            if (options_digraph && minus) gvalid = gvalid_sg = FALSE;
            options_digraph = !minus;
            minus = FALSE;
            break;

        case 'P':   /* set keep-group option */
            if (minus && options_keepgroup)
	    {
		options_keepgroup = FALSE;
		freeschreier(NULL,&generators);
	    }
	    else
	    {
		if ((d = getc(INFILE)) != 'P') ungetc(d,INFILE);
                options_keepgroup = TRUE;
		if (d == 'P')
		{
		    readvperm(INFILE,perm,prompt,n,&nperm);
		    if (nperm != n)
		    {
			fprintf(ERRFILE,"Incomplete permutation\n");
			FLUSHANDPROMPT;
		    }
		    else
			addpermutation(&generators,perm,n);
		}
	    }
            minus = FALSE;
            break;

        case '$':   /* set label origin */
            if ((d = getc(INFILE)) == '$')
                labelorg = oldorg;
            else
            {
                ungetc(d,INFILE);
                oldorg = labelorg;
                i = getint_sl(INFILE);
                if (i < 0)
		{
		    fprintf(ERRFILE,"labelorg must be >= 0\n");
		    FLUSHANDPROMPT;
		}
                else
                    labelorg = i;
            }
            break;

        case '?':   /* type options, etc. */
            minus = FALSE;
	    fprintf(outfile,"Mode=%s ",
		(mode==DENSE_MODE?"dense":
		 mode==SPARSE_MODE?"sparse":"Traces"));
            fprintf(outfile,"m=%d n=%d labelorg=%d",m,n,labelorg);
            if (!gvalid && !gvalid_sg)
                fprintf(outfile," g=undef");
            else if (gvalid)
            {
                uli = 0;
                for (i = 0, gp = g; i < n; ++i, gp += m) uli += setsize(gp,m);
                if (options_digraph) fprintf(outfile," arcs=%lu",uli);
                else                 fprintf(outfile," edges=%lu",uli/2);
            }
            else
            {
                uli = g_sg.nde;
                if (options_digraph) fprintf(outfile," arcs=%lu",uli);
                else                 fprintf(outfile," edges=%lu",uli/2);
            }
            fprintf(outfile," options=(%cc%ca%cm%cp%cd",
                        PM(options_getcanon),PM(options_writeautoms),
                        PM(options_writemarkers),PM(options_cartesian),
                        PM(options_digraph));
	    if (mode == TRACES_MODE)
		fprintf(outfile,"%cP",PM(options_keepgroup));
            if (umask & 31)
                fprintf(outfile," u=%d",umask&31);
            if (options_tc_level > 0)
                fprintf(outfile," y=%d",options_tc_level);
            if (options_mininvarlevel != 0 || options_maxinvarlevel != 0)
                fprintf(outfile," k=(%d,%d)",
                              options_mininvarlevel,options_maxinvarlevel);
            if (options_invararg > 0)
                fprintf(outfile," K=%d",options_invararg);
            if (multiplicity != 1 || mintime != 0.0)
		fprintf(outfile," M=%d/%.0f",multiplicity,mintime);
            fprintf(outfile,")\n");
            fprintf(outfile,"linelen=%d worksize=%d input_depth=%d",
                            options_linelength,worksize,curfile);
	    if (options_schreier > 0)
		fprintf(outfile," G=%d",options_schreier);
	    if (mode == TRACES_MODE)
	    {
		if (options_verbosity != 1)
		    fprintf(outfile," V=%d",options_verbosity);
		if (options_strategy != 0)
		    fprintf(outfile," S=%d",options_strategy);
	    }
            if (options_invarproc != 1)
                fprintf(outfile," invarproc=%s",
		  (mode == DENSE_MODE ?
		    invarproc[options_invarproc].name : 
		    invarproc[options_invarproc].name_sg)); 
            if (pvalid)
                fprintf(outfile,"; %d cell%s",SS(numcells,"","s"));
            else
                fprintf(outfile,"; 1 cell");
            fprintf(outfile,"\n");
            if (outfile != PROMPTFILE)
	    {
	        fprintf(outfile,"Mode=%s ",
		    (mode==DENSE_MODE?"dense":
		     mode==SPARSE_MODE?"sparse":"Traces"));
                fprintf(PROMPTFILE,"n=%d depth=%d labelorg=%d\n",
                     n,curfile,labelorg);
	    }
            break;

        case '&':   /* list the partition and possibly the quotient */
            if ((d = getc(INFILE)) == '&')
                doquot = TRUE;
            else
            {
                ungetc(d,INFILE);
                doquot = FALSE;
            }
            minus = FALSE;
            if (pvalid)
                putptn(outfile,lab,ptn,0,options_linelength,n);
            else
                fprintf(outfile,"unit partition\n");
            if (doquot)
            {
                if (!pvalid) unitptn(lab,ptn,&numcells,n);
		if (SPARSEREP(mode))
		    putquotient_sg(outfile,&g_sg,lab,ptn,0,options_linelength);
		else
                    putquotient(outfile,g,lab,ptn,0,options_linelength,m,n);
            }
            break;

        case 'h':   /* type help information */
        case 'H':
            minus = FALSE;
            help(PROMPTFILE,c == 'H');
            break;

        default:    /* illegal command */
            fprintf(ERRFILE,"'%c' is illegal - type 'h' for help\n",c);
	    FLUSHANDPROMPT;
            break;

        }  /* end of switch */

        if (flushing) fflush(outfile);
    }

    exit(0);
}

/*****************************************************************************
*                                                                            *
*  help(f,i) writes help information to file f (i = 0,1).                    *
*                                                                            *
*****************************************************************************/

static void
help(FILE *f, int i)
{
#define H(ss) fprintf(f," %s\n",ss);

if (i == 0)
{
H("Modes: An = dense, As = sparse, At = Traces; extra + to convert graph")
H("+- a : write automs     v,vv : write degrees    *=# : select invariant:")
H("   b : write canong      w=# : set worksize (units of 2m)")
H("+- c : canonise            x : run nauty          -1 = user-defined")
H("+- d : digraph or loops  y=# : set tc_level        0 = none")
H("   e : edit graph          z : write hashcode      1 = twopaths")
H("-f, f=#, f=[...] : set colours                     2 = adjtriang(K=0,1)")
H("   g : read graph        $=# : set origin          3 = triples")
H(" h,H : help               $$ : restore origin      4 = quadruples")
H("   i : refine              ? : type options        5 = celltrips")
H("   I : refine using invar  _ : compl  __ : conv    6 = cellquads")
H("   j : relabel randomly    % : Mathon doubling     7 = cellquins")
H("k=# # : set invar levels   & : type colouring      8 = distances(K)")
H(" K=# : set invar param    && : + quotient matrix   9 = indsets(K)")
H(" l=# : set line length   >ff : write to file      10 = cliques(K)")
H("+- m : write markers    >>ff : append to file     11 = cellcliq(K)")
H(" n=# : set order      ->/->> : close/flush output 12 = cellind(K)")
H("   o : write orbits      <ff : read from file     13 = adjacencies")
H("+- p : set autom format    @ : save canong        14 = cellfano")
H("   q : quit                # : canong = savedg?   15 = cellfano2")
H(" r,R : relabel/subgraph   ## : + write mapping    16 = refinvar")
H(" s=# : random g (p=1/#) sr=# : random reg   \"...\" : copy comment")
H("-G,G=# : schreier param  F=# : fix extra vertex  FF: fix target vertex")
H(" t,T : type graph          ! : ignore line        O : orbits->partition")
H("+- P : keep group         PP : add automorphism    Type H for more..")
}

if (i == 1)
{
H("Commands for g and e : ")
H("   There is always a \"current vertex\" v, initially first vertex.")
H("   # : add edge v-#       ; : increment v (exit if over limit)")
H("  -# : delete edge v-#   #: : set v := #")
H("   ? : list nbhs of v     . : exit")
H("Mode change:   An = dense nauty, As = sparse nauty, At = Traces")
H("Use An+, As+ or At+ to also convert graph between dense and sparse")
H("Command line argument -o options  allows a,c,d,m,p,l,G,P,w,y,$,A,V,M")
H("Syntax for f :  f=[2 3|4:9|10]  (rest in extra cell at right)")
H("               -f same as f=[], f=# same as f=[#]")
H("Syntax for r :  r 2:4 1 5;    (rest appended in order)")
H("r& relabels the graph and partition in order of the partition")
H("Syntax for R :  R 2:4 1 5;   or  -R 0 3 6:10;")
H("Syntax for PP :  PP 2:4 1 5 0; (must be complete)")
H("Arguments for u : 1=node,2=autom,4=level,16=ref,32=canon (add them)")
H("Accurate times: M=#/# set number of runs and minimum total cpu.")
}

}

/*****************************************************************************
*                                                                            *
*  usernode(g,lab,ptn,level,numcells,tc,code,m,n) is a simple version of the *
*  procedure named by options.usernodeproc.                                  *
*                                                                            *
*****************************************************************************/

static void
usernode(graph *g, int *lab, int *ptn, int level, int numcells,
     int tc, int code, int m, int n)
{
    int i;

    for (i = 0; i < level; ++i) PUTC('.',outfile);
    if (numcells == n)
        fprintf(outfile,"(n/%d)\n",code);
    else if (tc < 0)
        fprintf(outfile,"(%d/%d)\n",numcells,code);
    else
        fprintf(outfile,"(%d/%d/%d)\n",numcells,code,tc);
    if (firstpath) putptn(outfile,lab,ptn,level,options_linelength,n);
    if (numcells == n) firstpath = FALSE;
}

/*****************************************************************************
*                                                                            *
*  userautom(count,perm,orbits,numorbits,stabvertex,n) is a simple           *
*  version of the procedure named by options.userautomproc.                  *
*                                                                            *
*****************************************************************************/

static void
userautom(int count, int *perm, int *orbits,
      int numorbits, int stabvertex, int n)
{
    fprintf(outfile,
         "**userautomproc:  count=%d stabvertex=%d numorbits=%d\n",
         count,stabvertex+labelorg,numorbits);
    PUTORBITS(outfile,orbits,options_linelength,n);
}

/*****************************************************************************
*                                                                            *
*  userlevel(lab,ptn,level,orbits,stats,tv,index,tcellsize,numcells,cc,n)    *
*  is a simple version of the procedure named by options.userlevelproc.      *
*                                                                            *
*****************************************************************************/

static void
userlevel(int *lab, int *ptn, int level, int *orbits, statsblk *stats,
      int tv, int index, int tcellsize, int numcells, int cc, int n)
{
    fprintf(outfile,
        "**userlevelproc:  level=%d tv=%d index=%d tcellsize=%d cc=%d\n",
        level,tv+labelorg,index,tcellsize,cc);
    fprintf(outfile,"    nodes=%lu cells=%d orbits=%d generators=%d\n",
        stats->numnodes,numcells,stats->numorbits,stats->numgenerators);
}

/*****************************************************************************
*                                                                            *
*  usercanon(g,lab,canong,count,code,m,n)                                    *
*  is a simple version of the procedure named by options.userlevelproc.      *
*                                                                            *
*****************************************************************************/

static int
usercanon(graph *g, int *lab, graph *canong, int count, int code,
            int m, int n)
{
    fprintf(outfile,
	"**usercanonproc: count=%d code=%d\n",count,code);
    return 0;
}
