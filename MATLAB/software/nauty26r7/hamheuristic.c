/* hamheuristic.c  */

#define USAGE "hamheuristic [-sgu] [-vq] [-L#] [-t#] [infile [outfile]]"

#define HELPTEXT \
" Apply a heuristic for finding hamiltonian cycles.\n\
  Output those which are unsuccessful.\n\
\n\
    -s  force output to sparse6 format\n\
    -g  force output to graph6 format\n\
        If neither -s or -g are given, the output format is\n\
        determined by the header or, if there is none, by the\n\
        format of the first input graph.\n\
    -u  Suppress output to outfile, give statistics instead.\n\
\n\
    The output file will have a header if and only if the input file does.\n\
\n\
    -p  Be content with a hamiltonian path\n\
    -v  Give a cycle or path if one is found.\n\
    -L# Limit number of sideways steps (default 1000+5*n)\n\
    -t# Try # times (default 1)\n\
\n\
    -q  suppress auxiliary information\n"

#define DEBUG(X) X

/*************************************************************************/

#include "gtools.h"

#define NO 0
#define YES 1
#define TIMEOUT 2

static long seed = 314159265;

/**************************************************************************/

static int
hamheur(sparsegraph *sg, boolean pathok, unsigned long limit, int *cyc)
/* Try up to limit, fill cyc if YES and cyc!=NULL */
/* For payjok = TRUE, return when a hamiltonian path is found */
{
    DYNALLSTAT(int,path,path_sz);
    DYNALLSTAT(int,work,work_sz);
    DYNALLSTAT(int,pos,pos_sz);  /* position on path or -1 */
    int end0,end1,v0,v1;
    size_t *v;
    int n,*e,*d,len,d0,*e0,d1,*e1,dx,*ex;
    int i,ix,x,j,w,vext,exts;
    unsigned long count;
    boolean left,cycle;
    long ran;
 
    SG_VDE(sg,v,d,e);
    n = sg->nv;
    if (n < 3) return NO;

    DYNALLOC1(int,path,path_sz,2*n+4,"malloc");
    DYNALLOC1(int,pos,pos_sz,n,"malloc");
    DYNALLOC1(int,work,work_sz,n,"malloc");

    for (i = 0; i < n; ++i)
    {
	if (d[i] < 2) return NO;
	pos[i] = -1;
    }
    count = 0;
 
    v0 = KRAN(n);
    v1 = e[v[v0]+KRAN(d[v0])];
    end0 = n+1;
    end1 = n+2;
    path[end0] = v0; pos[v0] = end0;
    path[end1] = v1; pos[v1] = end1;
    len = 2;

    for (; count < limit;)
    {
     /* First look for an extension, but note if cycle */

	if (pathok && len == n)
	{
	    if (cyc)
	    {
	        j = 0;
	        for (i = end0; i <= end1; ++i)
		    cyc[j++] = path[i];
	    }
	    return YES;
	}

	v0 = path[end0];
	v1 = path[end1];

	exts = 0;
	d0 = d[v0];
        e0 = e + v[v0];
	cycle = FALSE;
        ran = NEXTRAN;
	for (i = 0; i < d0; ++i)
	{
	    w = e0[i];
	    if (pos[w] < 0)
	    {
		++exts;
		if (ran%exts == 0) {left = TRUE; vext = w;}
	    }
	    else if (w == v1) cycle = TRUE;
        }
	d1 = d[v1];
        e1 = e + v[v1];
        ran = NEXTRAN;
	for (i = 0; i < d1; ++i)
	{
	    w = e1[i];
	    if (pos[w] < 0) 
            {
                ++exts;
                if (ran%exts == 0) {left = FALSE; vext = w;}
            }
	}

	if (exts > 0)
	{
	    if (left)
	    {
		--end0;
		path[end0] = vext;
		pos[vext] = end0;
	    }
	    else
	    {
		++end1;
		path[end1] = vext;
		pos[vext] = end1;
	    }
	    ++len;
	    continue;
	}

     /* Can't extend but a cycle was found */

	if (cycle)
        {
	    if (len == n)
	    {
		if (cyc)
		{
		    j = 0;
		    for (i = pos[0]; i <= end1; ++i) cyc[j++] = path[i];
		    for (i = end0; i < pos[0]; ++i) cyc[j++] = path[i];
		}
		return YES;
	    }

	    ix = end0 + KRAN(len);
	    for (i = 0; i < len; ++i)      
            {
	        x = path[ix];

	        exts = 0;
	        dx = d[x];
	        ex = e + v[x];
		ran = NEXTRAN;
	        for (j = 0; j < dx; ++j)
	        {
		    w = ex[j];
		    if (pos[w] < 0) 
                    {
                        ++exts;
                        if (ran%exts == 0) vext = w;
                    }
	        }
		if (exts > 0) break;
		if (ix == end1) ix = end0; else ++ix;
            }
	    if (i == len) return NO;    /* isolated component */

	    work[0] = vext;
	    j = 1;
	    for (i = ix; i <= end1; ++i) work[j++] = path[i];
	    for (i = end0; i < ix; ++i) work[j++] = path[i];

	    for (i = 0; i < j; ++i)
	    {
		path[end0+i] = work[i];
                pos[work[i]] = end0+i;
	    }
	    ++end1;
	    ++len;
	    continue;
        }

     /* In the last resort, do a sideways move. */
    
	exts = 0;
	d0 = d[v0];
        e0 = e + v[v0];
        ran = NEXTRAN;
	for (i = 0; i < d0; ++i)
	{
	    w = e0[i];
	    if (pos[w] >= 0 && w != path[end0+1])
	    {
		++exts;
		if (ran%exts == 0) {left = TRUE; vext = w;}
	    }
        }
	d1 = d[v1];
        e1 = e + v[v1];
        ran = NEXTRAN;
	for (i = 0; i < d1; ++i)
	{
	    w = e1[i];
	    if (pos[w] >= 0 && w != path[end1-1]) 
            {
                ++exts;
                if (ran%exts == 0) {left = FALSE; vext = w;}
            }
	}

	if (left)
	{
	    i = end0;
            j = pos[vext]-1;
	}
	else
	{
	    i = pos[vext]+1;
	    j = end1;
	}
	for (; i < j; ++i, --j)
	{
	    w = path[i];
	    path[i] = path[j];
	    path[j] = w;
	    pos[path[i]] = i;
	    pos[path[j]] = j;
	}

	++count;
    }

    return TIMEOUT;
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    sparsegraph sg;
    int n,codetype;
    int argnum,i,j,outcode,tvalue;
    char *arg,sw;
    boolean badargs;
    boolean pswitch,sswitch,gswitch,qswitch,Lswitch,tswitch,vswitch,uswitch;
    long Lvalue;
    double t;
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    nauty_counter nin,nout,nNO,nYES,nTIMEOUT;
    int status;
    DYNALLSTAT(int,cyc,cyc_sz);

    HELP; PUTVERSION;

    INITSEED;
    ran_init(seed);

    uswitch = sswitch = gswitch = Lswitch = qswitch = FALSE;
    pswitch = vswitch = tswitch = FALSE;
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
                     SWBOOLEAN('s',sswitch)
                else SWBOOLEAN('g',gswitch)
                else SWBOOLEAN('u',uswitch)
                else SWBOOLEAN('p',pswitch)
                else SWBOOLEAN('q',qswitch)
                else SWBOOLEAN('v',vswitch)
		else SWLONG('L',Lswitch,Lvalue,"-L")
		else SWINT('t',tswitch,tvalue,"-t")
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

    if ((sswitch != 0) + (gswitch != 0) + (uswitch != 0) > 1)
        gt_abort(">E hamheuristic: -s, -g and -u are incompatible\n");

    if (!tswitch) tvalue = 1;

    if (badargs || argnum > 2)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!qswitch)
    {
        fprintf(stderr,">A hamheuristic");
        if (pswitch || sswitch || gswitch || vswitch || tswitch || Lswitch)
            fprintf(stderr," -");
        if (sswitch) fprintf(stderr,"s");
        if (gswitch) fprintf(stderr,"g");
        if (uswitch) fprintf(stderr,"u");
        if (vswitch) fprintf(stderr,"v");
        if (pswitch) fprintf(stderr,"p");
	if (Lswitch) fprintf(stderr,"L%ld",Lvalue);
	if (tswitch) fprintf(stderr,"t%d",tvalue);
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,FALSE,1);
    if (!infile) exit(1);
    if (!infilename) infilename = "stdin";

    NODIGRAPHSYET(codetype);

    if (!uswitch)
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

        if (sswitch || (!gswitch && (codetype&SPARSE6)))
            outcode = SPARSE6;
        else 
            outcode = GRAPH6;

        if (codetype&HAS_HEADER)
        {
            if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
            else                    writeline(outfile,GRAPH6_HEADER);
        }
    }

    nin = nout = 0;
    nNO = nYES = nTIMEOUT = 0;
    t = CPUTIME;
    SG_INIT(sg);

    while (read_sg(infile,&sg) != NULL)
    {
        ++nin;

	n = sg.nv;

        if (vswitch) { DYNALLOC1(int,cyc,cyc_sz,n,"malloc"); }
	else         cyc = NULL;

	status = TIMEOUT;
	for (i = 0; i < tvalue; ++i)
	{
	    status = hamheur(&sg,pswitch,(Lswitch?Lvalue:1000+5*n),cyc);
	    if (status != TIMEOUT) break;
	}

	if (status == NO) ++nNO;
        else if (status == YES) ++nYES;
        else ++nTIMEOUT;

        if (status != YES && !uswitch)
        {
            if (outcode == SPARSE6) writes6_sg(outfile,&sg);
            else if (outcode == GRAPH6) writeg6_sg(outfile,&sg);
	    ++nout;
	}

	if (vswitch)
	{
	    fprintf(stderr,">H " COUNTER_FMT ":",nin);
	    if (status == TIMEOUT)
		fprintf(stderr," timed out\n");
            else if (status == NO)
                    fprintf(stderr," disconnected or low degree\n");
	    else
	    {
		for (i = 0; i < n; ++i) fprintf(stderr," %d",cyc[i]);
		fprintf(stderr,"\n");
	    }
	}
    }
    t = CPUTIME - t;

    if (uswitch)
    {
	fprintf(stderr,">Z " COUNTER_FMT " graphs read from %s; "
           COUNTER_FMT " hamiltonian, " COUNTER_FMT " not, "
           COUNTER_FMT " timed out; %3.2f sec.\n",
           nin,infilename,nYES,nNO,nTIMEOUT,t);
    }
    else if (!qswitch)
    {
        fprintf(stderr,
                ">Z " COUNTER_FMT " graphs read from %s, "
                COUNTER_FMT " written to %s; %3.2f sec.\n",
                nin,infilename,nout,outfilename,t);
    }

    exit(0);
}
