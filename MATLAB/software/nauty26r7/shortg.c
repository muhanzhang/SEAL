/* shortg.c  version 3.1; B D McKay, July 15 2016. */

#define USAGE \
  "shortg [-qvkdu] [-i# -I#:# -K#] [-fxxx] [-S|-t] [-Tdir] [infile [outfile]]"

#define HELPTEXT \
"  Remove isomorphs from a file of graphs.\n\
\n\
    If outfile is omitted, it is taken to be the same as infile\n\
    If both infile and outfile are omitted, input will be taken\n\
            from stdin and written to stdout\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -s  force output to sparse6 format\n\
    -g  force output to graph6 format\n\
    -z  force output to digraph6 format\n\
        If none of -s, -z, -g are given, the output format is determined\n\
        by the header or, if there is none, by the format of the first\n\
        input graph. The output format determines the sorting order too.\n\
        As an exception, digraphs are always written in digraph6 format.\n\
    -S  Use sparse representation internally. Note that this changes the\n\
        canonical labelling. \n\
        Multiple edges are not supported.  One loop per vertex is ok.\n\
    -t  Use Traces.\n\
        Note that this changes the canonical labelling.\n\
        Multiple edges,loops, directed edges are not supported,\n\
        nor invariants.\n\
\n\
    -k  output graphs have the same labelling and format as the inputs.\n\
        Without -k, output graphs have canonical labelling.\n\
        -s, -g, -z are ineffective if -k is given.\n\
\n\
    -v  write to stderr a list of which input graphs correspond to which\n\
        output graphs. The input and output graphs are both numbered\n\
        beginning at 1.  A line like\n\
           23 : 30 154 78\n\
        means that inputs 30, 154 and 78 were isomorphic, and produced\n\
        output 23.\n\
\n\
    -d  include in the output only those inputs which are isomorphic\n\
        to another input.  If -k is specified, all such inputs are\n\
        included in their original labelling.  Without -k, only one\n\
        member of each nontrivial isomorphism class is written,\n\
        with canonical labelling.\n\
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
    -i#  select an invariant (1 = twopaths, 2 = adjtriang(K), 3 = triples,\n\
        4 = quadruples, 5 = celltrips, 6 = cellquads, 7 = cellquins,\n\
        8 = distances(K), 9 = indsets(K), 10 = cliques(K), 11 = cellcliq(K),\n\
       12 = cellind(K), 13 = adjacencies, 14 = cellfano, 15 = cellfano2,\n\
       16 = refinvar(K))\n\
    -I#:#  select mininvarlevel and maxinvarlevel (default 1:1)\n\
    -K#   select invararg (default 3)\n\
\n\
    -u  Write no output, just report how many graphs it would have output.\n\
        In this case, outfile is not permitted.\n\
    -Tdir  Specify that directory \"dir\" will be used for temporary disk\n\
        space by the sort subprocess.  The default is usually /tmp.\n\
    -q  Suppress auxiliary output\n"


/*************************************************************************/

#include "gtools.h" 
#include "nautinv.h"
#include "gutils.h"
#include "traces.h"

#if (HAVE_PIPE==0) || (HAVE_WAIT==0) || (HAVE_FORK==0)
 #error Forget it, either pipe(), wait() or fork() are not available
#endif

#if !HAVE_PID_T
typedef int pid_t;
#endif

#if !FDOPEN_DEC
FILE *fdopen(int, const char*);
#endif

#if SORT_NEWKEY == 0
#define SORTCOMMAND  SORTPROG,SORTPROG,"-u","+0","-1"
#define VSORTCOMMAND1  SORTPROG,SORTPROG
#define VSORTCOMMAND2  SORTPROG,SORTPROG,"+0","-1","+2"

#define SORTCOMMANDT  SORTPROG,SORTPROG,"-T",tempdir,"-u","+0","-1"
#define VSORTCOMMANDT1  SORTPROG,SORTPROG,"-T",tempdir
#define VSORTCOMMANDT2  SORTPROG,SORTPROG,"-T",tempdir,"+0","-1","+2"
#else
#define SORTCOMMAND  SORTPROG,SORTPROG,"-u","-k","1,1"
#define VSORTCOMMAND1  SORTPROG,SORTPROG
#define VSORTCOMMAND2  SORTPROG,SORTPROG,"-k","1,1","-k","3"

#define SORTCOMMANDT  SORTPROG,SORTPROG,"-T",tempdir,"-u","-k","1,1"
#define VSORTCOMMANDT1  SORTPROG,SORTPROG,"-T",tempdir
#define VSORTCOMMANDT2  SORTPROG,SORTPROG,"-T",tempdir,"-k","1,1","-k","3"
#endif

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

/**************************************************************************/

static pid_t
beginsort(FILE **sortin, FILE **sortout, char *tempdir,
          boolean vdswitch, boolean keep)
/* begin sort process, open streams for i/o to it, and return its pid */
{
    int pid;
    int inpipe[2],outpipe[2];

    if (pipe(inpipe) < 0 || pipe(outpipe) < 0)
        gt_abort(">E shortg: can't create pipes to sort process\n");

    if ((pid = fork()) < 0) gt_abort(">E shortg: can't fork\n");

    if (pid > 0)            /* parent */
    {
        close(inpipe[0]);
        close(outpipe[1]);
        if ((*sortin = fdopen(inpipe[1],"w")) == NULL)
            gt_abort(">E shortg: can't open stream to sort process\n");
        if ((*sortout = fdopen(outpipe[0],"r")) == NULL)
            gt_abort(">E shortg: can't open stream from sort process\n");
    }
    else                   /* child */
    {
        SET_C_COLLATION;

        close(inpipe[1]);
        close(outpipe[0]);
        if (dup2(inpipe[0],0) < 0 || dup2(outpipe[1],1) < 0)
            gt_abort(">E shortg: dup2 failed\n");

        if (tempdir == NULL)
        {
            if (vdswitch)
                if (keep) execlp(VSORTCOMMAND2,NULL);
                else      execlp(VSORTCOMMAND1,NULL);
            else
                execlp(SORTCOMMAND,NULL);
        }
        else
        {
            if (vdswitch)
                if (keep) execlp(VSORTCOMMANDT2,NULL);
                else      execlp(VSORTCOMMANDT1,NULL);
            else
                execlp(SORTCOMMANDT,NULL);
        }
        gt_abort(">E shortg: can't start sort process\n");
    }

    return pid;
}

/**************************************************************************/

static void
tosort(FILE *f, char *cdstr, char *dstr, nauty_counter index)
/* write one graph to sort process 
   cdstr = canonical string 
   dstr = optional original string
   index = optional index number */
{
    int i;
    char buff[30];

    for (i = 0; cdstr[i] != '\n'; ++i) {}
    cdstr[i] = '\0';
    writeline(f,cdstr);

    if (dstr != NULL)
    {
        writeline(f," ");
        for (i = 0; dstr[i] != '\n'; ++i) {}
        dstr[i] = '\0';
        writeline(f,dstr);
    }

    if (index > 0)
    {
#if LONG_LONG_COUNTERS
        sprintf(buff,"\t%11llu\n",index);
#else
        sprintf(buff,"\t%10lu\n",index);
#endif
        writeline(f,buff);
    }
    else
        writeline(f,"\n");
}

/**************************************************************************/

static boolean
fromsort(FILE *f, char **cdstr, char **dstr, nauty_counter *index)
/* read one graph from sort process */
{
    int j;
    char *s;

    if ((s = gtools_getline(f)) == NULL) return FALSE;

    *cdstr = s;
    for (j = 0; s[j] != ' ' && s[j] != '\t' && s[j] != '\n'; ++j) {}

    if (s[j] == ' ')
    {
        s[j] = '\0';
        *dstr = &s[j+1];
        for (++j; s[j] != '\t' && s[j] != '\n'; ++j) {}
    }
    else
        *dstr = NULL;

    if (s[j] == '\t')
    {
#if LONG_LONG_COUNTERS
        if (sscanf(&s[j+1],"%llu",index) != 1)
#else
        if (sscanf(&s[j+1],"%lu",index) != 1)
#endif
            gt_abort(">E shortg: index field corrupted\n");
    }
    else
        *index = 0;
    s[j] = '\0';

    return TRUE;
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    FILE *sortin,*sortout;
    int status,loops;
    char *dstr,*cdstr,*prevdstr,*prevcdstr;
    char sw,*fmt;
    boolean badargs,quiet,vswitch,dswitch,kswitch,format,uswitch;
    boolean iswitch,Iswitch,Kswitch,Tswitch;
    boolean zswitch,sswitch,gswitch,Sswitch,tswitch;
    boolean digraph;
    nauty_counter numread,prevnumread,numwritten,classsize;
    int m,n,i,ii,argnum,line;
    int outcode,codetype;
    int inv,mininvarlevel,maxinvarlevel,invararg;
    long minil,maxil;
    pid_t sortpid;
    graph *g;
    char *arg,*tempdir;
    sparsegraph sg,sh;
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

    HELP; PUTVERSION;

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    infilename = outfilename = NULL;
    dswitch = format = quiet = vswitch = kswitch = uswitch = FALSE;
    zswitch = sswitch = gswitch = Sswitch = Tswitch = FALSE;
    tswitch = iswitch = Iswitch = Kswitch = FALSE;
    tempdir = NULL;
    inv = 0;

     /* parse argument list */

    argnum = 0;
    badargs = FALSE;

    for (i = 1; !badargs && i < argc; ++i)
    {
        arg = argv[i];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                     SWBOOLEAN('q',quiet)
                else SWBOOLEAN('v',vswitch)
                else SWBOOLEAN('k',kswitch)
                else SWBOOLEAN('d',dswitch)
                else SWBOOLEAN('u',uswitch)
                else SWBOOLEAN('s',sswitch)
                else SWBOOLEAN('g',gswitch)
                else SWBOOLEAN('z',zswitch)
                else SWBOOLEAN('S',Sswitch)
                else SWBOOLEAN('t',tswitch)
                else SWINT('i',iswitch,inv,"shortg -i")
                else SWINT('K',Kswitch,invararg,"shortg -K")
                else SWRANGE('I',":-",Iswitch,minil,maxil,"shortg -I")
                else if (sw == 'f')
                {
                    format = TRUE;
                    fmt = arg;
                    break;
                }
                else if (sw == 'T')
                {
                    Tswitch = TRUE;
                    tempdir = arg;
                    break;
                }
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

    if (iswitch && inv == 0) iswitch = FALSE;

    if (strcmp(SORTPROG,"no_sort_found") == 0)
        gt_abort(">E shortg: no sort program known\n");
    if (uswitch && outfilename != NULL)
        gt_abort(">E shortg: -u and outfile are incompatible\n");
    if ((sswitch!=0) + (gswitch!=0) + (zswitch!=0) + (kswitch!=0) > 1)
        gt_abort(">E shortg: -sgzk are incompatible\n");
    if (Tswitch && *tempdir == '\0')
        gt_abort(">E shortg: -T needs a non-empty argument\n");
    if (tswitch && (format || Sswitch))
        gt_abort(">E shortg: -t is incompatible with -S and -f \n");
    if (iswitch && (inv > 16))
        gt_abort(">E shortg: -i value must be 0..16\n");
    if (Sswitch && iswitch && invarproc[inv].entrypoint_sg == NULL)
        gt_abort(
          ">E shortg: that invariant is not available in sparse mode\n");

    if (argnum == 1 && !uswitch) outfilename = infilename;

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

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!quiet)
    {
        fprintf(stderr,">A shortg");
        if (uswitch || kswitch || vswitch || format || Tswitch || zswitch
                  || tswitch  || sswitch || gswitch || Sswitch || iswitch)
        fprintf(stderr," -");
        if (sswitch) fprintf(stderr,"s");
        if (gswitch) fprintf(stderr,"g");
        if (zswitch) fprintf(stderr,"z");
        if (Sswitch) fprintf(stderr,"S");
        if (tswitch) fprintf(stderr,"t");
        if (kswitch) fprintf(stderr,"k");
        if (vswitch) fprintf(stderr,"v");
        if (dswitch) fprintf(stderr,"d");
        if (uswitch) fprintf(stderr,"u");
        if (iswitch)
            fprintf(stderr,"i=%s[%d:%d,%d]",invarproc[inv].name,
                    mininvarlevel,maxinvarlevel,invararg);
        if (format) fprintf(stderr,"f%s",fmt);
        if (format && Tswitch) fprintf(stderr," -");
        if (Tswitch) fprintf(stderr,"T%s",tempdir);
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
    }

     /* open input file */

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,FALSE,1);
    if (!infile) exit(1);
    if (!infilename) infilename = "stdin";

    if (gswitch)      outcode = GRAPH6;
    else if (sswitch) outcode = SPARSE6;
    else if (zswitch) outcode = DIGRAPH6;
    else if ((codetype&GRAPH6))   outcode = GRAPH6;
    else if ((codetype&SPARSE6))  outcode = SPARSE6;
    else if ((codetype&DIGRAPH6)) outcode = DIGRAPH6;
    else
    {
        outcode = GRAPH6;
        fprintf(stderr,
             ">W shortg doesn't handle this graph format, writing graph6.\n");
    }

    if (sswitch || (!gswitch && (codetype&SPARSE6))) outcode = SPARSE6;
    else                                             outcode = GRAPH6;

#if defined(SIG_IGN) && defined(SIGPIPE)
    signal(SIGPIPE,SIG_IGN);        /* process pipe errors ourselves */
#endif

     /* begin sort in a subprocess */

    sortpid = beginsort(&sortin,&sortout,(Tswitch?tempdir:NULL),
                                            dswitch||vswitch,kswitch);

  /* feed input graphs, possibly relabelled, to sort process */

    numread = 0;

    if (Sswitch)
    {
	SG_INIT(sg);
	SG_INIT(sh);
        while (TRUE)
        {
	    if (read_sgg_loops(infile,&sg,&loops,&digraph) == NULL) break;
            dstr = readg_line;
            ++numread;
	    n = sg.nv;
#if MAXN
	    if (n > MAXN) gt_abort(">E shortg: graph larger than MAXN read\n");
#endif
            m = SETWORDSNEEDED(n);
            SG_ALLOC(sh,n,sg.nde,"shortg");
            fcanonise_inv_sg(&sg,m,n,&sh,format?fmt:NULL,
		invarproc[inv].entrypoint_sg,
                mininvarlevel,maxinvarlevel,invararg,loops>0||digraph);
            sortlists_sg(&sh);
	    if (outcode == DIGRAPH6 || digraph) cdstr = sgtod6(&sh);
            else if (outcode == SPARSE6)        cdstr = sgtos6(&sh);
            else                                cdstr = sgtog6(&sh);

            tosort(sortin,cdstr,kswitch ? dstr : NULL,vswitch ? numread : 0);
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
            if (loops > 0) gt_abort(
               ">E shortg: Traces does not allow loops or directed edges\n");
            ++numread;
            n = sg.nv;
#if MAXN
	    if (n > MAXN) gt_abort(">E shortg: graph larger than MAXN read\n");
#endif
            DYNALLOC1(int,lab,lab_sz,n,"traces@shortg");
            DYNALLOC1(int,ptn,ptn_sz,n,"traces@shortg");
            DYNALLOC1(int,orbits,orbits_sz,n,"traces@shortg");
            SG_ALLOC(sh,n,sg.nde,"labelg");
            for (ii = 0; ii < n; ++ii) { lab[ii] = ii; ptn[ii] = 1; }
            ptn[n-1] = 0;
            Traces(&sg,lab,ptn,orbits,&traces_opts,&traces_stats,&sh);
            sortlists_sg(&sh);
	    if (outcode == DIGRAPH6 || digraph) cdstr = sgtod6(&sh);
            else if (outcode == SPARSE6)        cdstr = sgtos6(&sh);
            else                                cdstr = sgtog6(&sh);
	    tosort(sortin,cdstr,kswitch ? dstr : NULL,vswitch ? numread : 0);
        }
    }
    else
    {
        while (TRUE)
        {
            if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
            dstr = readg_line;
            ++numread;
            loops = loopcount(g,m,n);
#if !MAXN
            DYNALLOC2(graph,h,h_sz,n,m,"shortg");
#endif
            fcanonise_inv(g,m,n,h,format?fmt:NULL,
                invarproc[inv].entrypoint,mininvarlevel,maxinvarlevel,
                invararg, loops>0||digraph);
	    if (outcode == DIGRAPH6 || digraph) cdstr = ntod6(h,m,n);
            else if (outcode == SPARSE6)        cdstr = ntos6(h,m,n);
            else                                cdstr = ntog6(h,m,n);

            tosort(sortin,cdstr,kswitch ? dstr : NULL,vswitch ? numread : 0);
            FREES(g);
        }
    }
    fclose(sortin);
    fclose(infile);

     /* open output file */

    if (uswitch)
        outfilename = "<none>";
    else if (outfilename == NULL || outfilename[0] == '-' || is_pipe)
    {
        outfile = stdout;
        outfilename = "stdout";
    }
    else
    {
        if ((outfile = fopen(outfilename,"w")) == NULL)
        {
            fprintf(stderr,
                ">E shortg: can't open %s for writing\n",outfilename);
            gt_abort(NULL);
        }
    }

    if (!uswitch && (codetype&HAS_HEADER))
    {
        if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
        else                    writeline(outfile,GRAPH6_HEADER);
    }

    if (!quiet)
        fprintf(stderr,">Z " COUNTER_FMT " graphs read from %s\n",
                numread,infilename);

     /* collect output from sort process and write to output file */

    prevcdstr = prevdstr = NULL;
    numwritten = 0;
    if (dswitch)
    {
        classsize = 0;
        while (fromsort(sortout,&cdstr,&dstr,&numread))
        {
            if (classsize == 0 || strcmp(cdstr,prevcdstr) != 0)
                classsize = 1;
            else
            {
                ++classsize;
                if (classsize == 2)
                {
                    ++numwritten;
                    if (!uswitch)
                    {
                        writeline(outfile,kswitch ? prevdstr : prevcdstr);
                        writeline(outfile,"\n");
                    }
                    if (kswitch)
                    {
                        ++numwritten;
                        if (!uswitch)
                        {
                            writeline(outfile,kswitch ? dstr : cdstr);
                            writeline(outfile,"\n");
                        }
                    }
                    if (vswitch)
                    {
                        fprintf(stderr,"\n");
#if LONG_LONG_COUNTERS
                        fprintf(stderr,"%3llu : %3llu %3llu",
                            numwritten,prevnumread,numread);
#else
                        fprintf(stderr,"%3lu : %3lu %3lu",
                            numwritten,prevnumread,numread);
#endif
                        line = 1;
                    }
                }
                else
                {
                    if (kswitch)
                    {
                        ++numwritten;
                        if (!uswitch)
                        {
                            writeline(outfile,kswitch ? dstr : cdstr);
                            writeline(outfile,"\n");
                        }
                    }
                    if (vswitch)
                    {
                        if (line == 15)
                        {
                            line = 0;
                            fprintf(stderr,"\n     ");
                        }
#if LONG_LONG_COUNTERS
                        fprintf(stderr," %3llu",numread);
#else
                        fprintf(stderr," %3lu",numread);
#endif
                        ++line;
                    }
                }
            }
            if (prevcdstr) FREES(prevcdstr);
            prevcdstr = stringcopy(cdstr);
            if (prevdstr) FREES(prevdstr);
            if (kswitch) prevdstr = stringcopy(dstr);
            prevnumread = numread;
        }
        if (vswitch) fprintf(stderr,"\n\n");
    }
    else if (vswitch)
    {
        while (fromsort(sortout,&cdstr,&dstr,&numread))
        {
            if (numwritten == 0 || strcmp(cdstr,prevcdstr) != 0)
            {
                ++numwritten;
                if (!uswitch)
                {
                    writeline(outfile,kswitch ? dstr : cdstr);
                    writeline(outfile,"\n");
                }
                fprintf(stderr,"\n");
#if LONG_LONG_COUNTERS
                fprintf(stderr,"%3llu : %3llu",numwritten,numread);
#else
                fprintf(stderr,"%3lu : %3lu",numwritten,numread);
#endif
                line = 1;
            }
            else
            {
                if (line == 15)
                {
                    line = 0;
                    fprintf(stderr,"\n     ");
                }
#if LONG_LONG_COUNTERS
                fprintf(stderr," %3llu",numread);
#else
                fprintf(stderr," %3lu",numread);
#endif
                ++line;
            }
            if (prevcdstr) FREES(prevcdstr);
            prevcdstr = stringcopy(cdstr);
        }
        fprintf(stderr,"\n\n");
    }
    else
    {
        while (fromsort(sortout,&cdstr,&dstr,&numread))
        {
            ++numwritten;
            if (!uswitch)
            {
                writeline(outfile,kswitch ? dstr : cdstr);
                writeline(outfile,"\n");
            }
        }
    }

    fclose(sortout);
    if (!uswitch) fclose(outfile);

    if (!quiet)
    {
        if (uswitch)
            fprintf(stderr,">Z " COUNTER_FMT " graphs produced\n",numwritten);
        else
            fprintf(stderr,
                  ">Z " COUNTER_FMT " graphs written to %s\n",numwritten,outfilename);
    }

     /* check that the subprocess exitted properly */

    while (wait(&status) != sortpid) {}

#if (defined(WIFSIGNALED) || defined(WTERMSIG)) && defined(WEXITSTATUS)
#ifdef WIFSIGNALED
    if (WIFSIGNALED(status) && WTERMSIG(status) != 0)
#else
    if (WTERMSIG(status) != 0)
#endif
    {
        fprintf(stderr,">E shortg: sort process killed (signal %d)\n",
                      WTERMSIG(status)); 
        gt_abort(NULL);
    }   
    else if (WEXITSTATUS(status) != 0)
    {
        fprintf(stderr,
                ">E shortg: sort process exited abnormally (code %d)\n",
                WEXITSTATUS(status));
        gt_abort(NULL);
    }
#endif

    exit(0);
}
