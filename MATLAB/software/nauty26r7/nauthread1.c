/* nauthread1.c

   This program tests dense nauty running in multiple threads.
   It must be linked with nauty as configured with 
   --enable-tls and will only run on systems which support
   thread-local storage.
*/

#include <pthread.h>
#include "nauty.h"   
/* MAXN=0 is defined by nauty.h, which implies dynamic allocation */

#define THREADS 1000    /* Total number of threads to run */
#define ATONCE 20       /* Number of threads to run at once */
#define GRAPHSIZE 200   /* Least graph size to use */

typedef struct
{
    int n;
    boolean writeautoms;
} params;   /* Used to pass parameters to the thread */

static void*
runit(void * threadarg)          /* Main routine for one thread */
{
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    set *gv;

    int n,m,v;

    n = ((params*)threadarg)->n;

 /* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
    Here we change those options that we want to be different from the
    defaults.  writeautoms=TRUE causes automorphisms to be written.     */

    options.writeautoms = ((params*)threadarg)->writeautoms;

    m = SETWORDSNEEDED(n);

 /* The following optional call verifies that we are linking
    to compatible versions of the nauty routines.            */

    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");

 /* Now we will make a polygon of n vertices */

    EMPTYGRAPH(g,m,n);
    for (v = 0; v < n; ++v) ADDONEEDGE(g,v,(v+1)%n,m);

    if (options.writeautoms)
         printf("Generators for Aut(C[%d]):\n",n);
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);

    if (options.writeautoms)
    {
        printf("order = ");
        writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
        printf("\n");
    }
    if (stats.numorbits != 1 || stats.grpsize1 != 2*n)
        fprintf(stderr,">E group error\n");

 /* If we are using multiple threads, we need to free all the dynamic
    memory we have allocated.  We don't have to do this after each 
    call to nauty, just once before the thread finishes. */

    DYNFREE(g,g_sz);
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    nauty_freedyn();
    nautil_freedyn();
    naugraph_freedyn();  /* Use nausparse_freedyn() instead if
                            sparse format is being used. */

    return NULL;
}


int
main(int argc, char *argv[])
{
    int n,ret;
    pthread_t thread[THREADS];
    params par[THREADS];
    int started,finished;

#if !HAVE_TLS
    fprintf(stderr,">E This program needs to be linked with a version\n");
    fprintf(stderr,"  of nauty successfully configured with --enable-tls.\n");
    exit(1);
#endif

    for (started = finished = 0; finished < THREADS; )
    {
	if (started == THREADS || started-finished == ATONCE)
	{
	    if ((ret = pthread_join(thread[finished],NULL)) != 0)
	    {
                fprintf(stderr,">E Thread joining failed, code=%d\n",ret);    
                exit(1);
            }
	    ++finished;
	}
	else
	{
		/* We vary the graph size a bit as it tests the
                   thread independence better. */
	    par[started].n = GRAPHSIZE + (started % 17);
	    par[started].writeautoms = FALSE;

            if ((ret = pthread_create(&thread[started],NULL,
					runit,&par[started])) != 0)
            {
                fprintf(stderr,">E Thread creation failed, code=%d\n",ret);       
                exit(1);
            }
	    ++started;
	}
    }

    exit(0);
}
