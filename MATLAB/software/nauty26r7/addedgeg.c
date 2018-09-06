/* addedgeg.c   nauty version 2.6; B D McKay, Jan 2013. */

#define USAGE "addedgeg [-lq] [-D#] [-btfF] [-z#] [infile [outfile]]"

#define HELPTEXT \
" For each edge nonedge e, output G+e if it satisfies certain conditions\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -l  Canonically label outputs\n\
    -D# Specify an upper bound on the maximum degree of the output\n\
    -b  Output has no new cycles of odd length\n\
    -t  Output has no new 3-cycle if input doesn't\n\
    -f  Output has no new 4-cycle if input doesn't\n\
    -F  Output has no new 5-cycle if input doesn't\n\
    -z# Output has no new cycles of length less than #\n\
        -btfFz can be used in arbitrary combinations\n\
    -q  Suppress auxiliary information\n"

/*************************************************************************/

#include "gtools.h" 
#include "gutils.h"

/*************************************************************************/

static void
no3path(graph *g, int m, int n, int v, int *dist)
/* For each i, set dist[i]=0 if there is a 3-path from v to i */
{
	set *gv,*gv1,*gv2;
	int v1,v2,v3;

	gv = GRAPHROW(g,v,m);
	for (v1 = -1; (v1 = nextelement(gv,m,v1)) >= 0; )
	{
	    gv1 = GRAPHROW(g,v1,m);
	    for (v2 = -1; (v2 = nextelement(gv1,m,v2)) >= 0; )
            {
		if (v2 == v) continue;
                gv2 = GRAPHROW(g,v2,m);
		for (v3 = -1; (v3 = nextelement(gv2,m,v3)) >= 0; )
                    if (v3 != v && v3 != v1) dist[v3] = 0;
	    }
	}
}

/*************************************************************************/

static void
no4path(graph *g, int m, int n, int v, int *dist)
/* For each i, set dist[i]=0 if there is a 4-path from v to i */
{
        set *gv,*gv1,*gv2,*gv3;
        int v1,v2,v3,v4;

        gv = GRAPHROW(g,v,m);
        for (v1 = -1; (v1 = nextelement(gv,m,v1)) >= 0; )
        {   
            gv1 = GRAPHROW(g,v1,m);
            for (v2 = -1; (v2 = nextelement(gv1,m,v2)) >= 0; )
            {   
                if (v2 == v) continue;
                gv2 = GRAPHROW(g,v2,m);
                for (v3 = -1; (v3 = nextelement(gv2,m,v3)) >= 0; )
		{
		    if (v3 == v || v3 == v1) continue;
		    gv3 = GRAPHROW(g,v3,m);
		    for (v4 = -1; (v4 = nextelement(gv3,m,v4)) >= 0; )
			if (v4 != v && v4 != v1 && v4 != v2) dist[v4] = 0;
		}
            }
        }
}

/*************************************************************************/

int
main(int argc, char *argv[])
{
        char *infilename,*outfilename;
        FILE *infile,*outfile;
        boolean badargs,dolabel,quiet,Dswitch;
	boolean bswitch,tswitch,fswitch,Fswitch,zswitch;
	int mincycle;
	int i,j,m,n,v,w,argnum;
	int codetype,outcode;
	graph *g,*gq;
	nauty_counter nin,nout;
        char *arg,sw;
	setword *gv,*gw;
	int maxdeg,actmaxdeg,degv;
	double t;
#if MAXN
	graph h[MAXN*MAXM];
	int deg[MAXN],dist[MAXN];
	boolean okdist[MAXN+1];
#else
	DYNALLSTAT(graph,h,h_sz);
	DYNALLSTAT(int,deg,deg_sz);
	DYNALLSTAT(boolean,okdist,okdist_sz);
	DYNALLSTAT(int,dist,dist_sz);
#endif

	HELP; PUTVERSION;

        infilename = outfilename = NULL;
	Dswitch = dolabel = quiet = zswitch = FALSE;
	bswitch = tswitch = fswitch = Fswitch = FALSE;

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
		         SWBOOLEAN('l',dolabel)
		    else SWBOOLEAN('q',quiet)
		    else SWBOOLEAN('b',bswitch)
		    else SWBOOLEAN('t',tswitch)
		    else SWBOOLEAN('f',fswitch)
		    else SWBOOLEAN('F',Fswitch)
                    else SWINT('z',zswitch,mincycle,">E addedgeg -z")
                    else SWINT('D',Dswitch,maxdeg,">E addedgeg -D")
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
	    fprintf(stderr,">A addedgeg");
	    if (dolabel) fprintf(stderr," -l");
	    if (Dswitch) fprintf(stderr," -D%d",maxdeg);
	    if (bswitch || tswitch || fswitch || Fswitch || zswitch)
	    {
		fprintf(stderr," -");
		if (bswitch) fprintf(stderr,"b");
		if (tswitch) fprintf(stderr,"t");
		if (fswitch) fprintf(stderr,"f");
		if (Fswitch) fprintf(stderr,"F");
		if (zswitch) fprintf(stderr,"z%d",mincycle);
	    }
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

	if (codetype&SPARSE6) outcode = SPARSE6;
	else                  outcode = GRAPH6;

	if (codetype&HAS_HEADER)
	{
	    if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
	    else    		    writeline(outfile,GRAPH6_HEADER);
	}

	if (!Dswitch) maxdeg = NAUTY_INFINITY;

	if (dolabel) nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

	nin = nout = 0;
	t = CPUTIME;
	while (TRUE)
	{
	    if ((g = readg(infile,NULL,0,&m,&n)) == NULL) break;
	    ++nin;

#if !MAXN
	    DYNALLOC1(int,deg,deg_sz,n,"addedgeg");
	    DYNALLOC1(boolean,okdist,okdist_sz,n+1,"addedgeg");
	    DYNALLOC1(int,dist,dist_sz,n,"addedgeg");
#endif

	    actmaxdeg = n;
	    for (v = 0, gv = g; v < n; ++v, gv += m)
	    {
		degv = 0;
		for (i = 0; i < m; ++i)
		    degv += POPCOUNT(gv[i]);
		if (degv < actmaxdeg) actmaxdeg = degv;
		deg[v] = degv;
	    }

	    if (actmaxdeg > maxdeg) continue;

	    okdist[0] = okdist[1] = FALSE;
	    for (i = 2; i <= n; ++i) okdist[i] = TRUE;

	    if (bswitch)
	        for (i = 2; i < n; i += 2) okdist[i] = FALSE;

	    if (tswitch && n >= 3) okdist[2] = FALSE;
	    if (fswitch && n >= 4) okdist[3] = FALSE;
	    if (Fswitch && n >= 5) okdist[4] = FALSE;
	    if (zswitch)
		for (i = 2; i < mincycle-1 && i < n; ++i) okdist[i] = FALSE;

	    for (v = 0, gv = g; v < n-1; ++v, gv += m)
	    {
		if (deg[v] >= maxdeg) continue;

		find_dist(g,m,n,v,dist);

		for (w = v+1; w < n; ++w)
		    dist[w] = okdist[dist[w]] && (deg[w] < maxdeg);

		if (fswitch) no3path(g,m,n,v,dist);
		if (Fswitch) no4path(g,m,n,v,dist);

	    	for (w = v+1; w < n; ++w)
	        {
		    if (!dist[w]) continue;

		    gw = GRAPHROW(g,w,m);
		    ADDELEMENT(gv,w);
		    ADDELEMENT(gw,v);
		    gq = g;
	    
	            if (dolabel)
	            {
#if !MAXN
		        DYNALLOC2(graph,h,h_sz,n,m,"addedgeg");
#endif
	 	        fcanonise(g,m,n,h,NULL,FALSE);  /*FIXME (loops)*/
		        gq = h;
	            }
	            if (outcode == SPARSE6) writes6(outfile,gq,m,n);
	            else                    writeg6(outfile,gq,m,n);
		    ++nout;
		    DELELEMENT(gv,w);
                    DELELEMENT(gw,v);
		}
	    }
	    FREES(g);
	}
	t = CPUTIME - t;

        if (!quiet)
            fprintf(stderr,
              ">Z  " COUNTER_FMT " graphs read from %s, "
                           COUNTER_FMT " written to %s; %3.2f sec.\n",
                    nin,infilename,nout,outfilename,t);

	 exit(0);
}
