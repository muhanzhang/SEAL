/* This program prints generators for the automorphism group of an
   n-vertex polygon, where n is a number supplied by the user.

   This version uses dynamic allocation.
*/

#include "nauty.h"   
/* MAXN=0 is defined by nauty.h, which implies dynamic allocation */

int
main(int argc, char *argv[])
{
  /* DYNALLSTAT declares a pointer variable (to hold an array when it
     is allocated) and a size variable to remember how big the array is.
     Nothing is allocated yet.  */
 
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    int n,m,v;
    set *gv;

/* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
   Here we change those options that we want to be different from the
   defaults.  writeautoms=TRUE causes automorphisms to be written. */

    options.writeautoms = TRUE;

    while (1)
    {
        printf("\nenter n : ");
        if (scanf("%d",&n) == 1 && n > 0)
        {

     /* The nauty parameter m is a value such that an array of
        m setwords is sufficient to hold n bits.  The type setword
        is defined in nauty.h.  The number of bits in a setword is
        WORDSIZE, which is 16, 32 or 64.  Here we calculate
        m = ceiling(n/WORDSIZE). */

            m = SETWORDSNEEDED(n);

         /* The following optional call verifies that we are linking
            to compatible versions of the nauty routines. */

            nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

         /* Now that we know how big the graph will be, we allocate
          * space for the graph and the other arrays we need. */

            DYNALLOC2(graph,g,g_sz,m,n,"malloc");
            DYNALLOC1(int,lab,lab_sz,n,"malloc");
            DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
            DYNALLOC1(int,orbits,orbits_sz,n,"malloc");

            EMPTYGRAPH(g,m,n);
            for (v = 0; v < n; ++v) ADDONEEDGE(g,v,(v+1)%n,m);

            printf("Generators for Aut(C[%d]):\n",n);
            densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);

            printf("order = ");
            writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
            printf("\n");
        }
        else
            break;
    }

    exit(0);
}
