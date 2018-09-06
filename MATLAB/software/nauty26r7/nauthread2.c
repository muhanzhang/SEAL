/* nauthread2.c

   This program tests sparse nauty running in multiple threads.
   It must be linked with nauty as configured with 
   --enable-tls and will only run on systems which support
   thread-local storage.
*/

#include <pthread.h>
#include "nauty.h"   
#include "nausparse.h"
/* MAXN=0 is defined by nauty.h, which implies dynamic allocation */

#define THREADS 1000    /* Total number of threads to run */
#define ATONCE 20       /* Number of threads to run at once */
#define GRAPHSIZE 500   /* Least graph size to use */

typedef struct
{
    int n;
    boolean writeautoms;
} params;   /* Used to pass parameters to the thread */

static void*
runit(void * threadarg)          /* Main routine for one thread */
{
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DEFAULTOPTIONS_SPARSEGRAPH(options);
    statsblk stats;
    sparsegraph sg;   /* Declare sparse graph structure */

    int n,m,i;

    n = ((params*)threadarg)->n;
    options.writeautoms = ((params*)threadarg)->writeautoms;

 /* Initialise sparse graph structure. */

    SG_INIT(sg);

    m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");

 /* SG_ALLOC makes sure that the v,d,e fields of a sparse graph
    structure point to arrays that are large enough.  This only
    works if the structure has been initialised. */

    SG_ALLOC(sg,n,2*n,"malloc");

    sg.nv = n;              /* Number of vertices */
    sg.nde = 2*n;           /* Number of directed edges */

    for (i = 0; i < n; ++i)
    {
        sg.v[i] = 2*i;
        sg.d[i] = 2;
        sg.e[2*i] = (i+n-1)%n;      /* edge i->i-1 */
        sg.e[2*i+1] = (i+n+1)%n;    /* edge i->i+1 */
    }

    if (options.writeautoms)
        printf("Generators for Aut(C[%d]):\n",n);
    sparsenauty(&sg,lab,ptn,orbits,&options,&stats,NULL);

    if (options.writeautoms)
    {
        printf("Automorphism group size = ");
        writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
        printf("\n");
    }
    if (stats.numorbits != 1 || stats.grpsize1 != 2*n)
        fprintf(stderr,">E group error\n");

 /* If we are using multiple threads, we need to free all the dynamic
    memory we have allocated.  We don't have to do this after each 
    call to nauty, just once before the thread finishes. */

    SG_FREE(sg);
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    nauty_freedyn();
    nautil_freedyn();
    nausparse_freedyn();  /* Use naugraph_freedyn() instead if
                            dense format is being used. */

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
    fprintf(stderr,"  of nauty successfully configured using --enable-tls.\n");
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
