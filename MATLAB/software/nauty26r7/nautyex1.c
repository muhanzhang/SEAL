/* This program prints generators for the automorphism group of an
   n-vertex polygon, where n is a number supplied by the user.

   This version uses a fixed limit for MAXN.
*/

#define MAXN 1000    /* Define this before including nauty.h */
#include "nauty.h"   /* which includes <stdio.h> and other system files */

int
main(int argc, char *argv[])
{
    graph g[MAXN*MAXM];
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    int n,m,v;

 /* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
    Here we change those options that we want to be different from the
    defaults.  writeautoms=TRUE causes automorphisms to be written.     */

    options.writeautoms = TRUE;

    while (1)
    {
        printf("\nenter n : ");
        if (scanf("%d",&n) != 1 || n <= 0)    /* Exit if EOF or bad number */
            break;

        if (n > MAXN)
        {
            printf("n must be in the range 1..%d\n",MAXN);
            exit(1);
        }

     /* The nauty parameter m is a value such that an array of
        m setwords is sufficient to hold n bits.  The type setword
        is defined in nauty.h.  The number of bits in a setword is
        WORDSIZE, which is 16, 32 or 64.  Here we calculate
        m = ceiling(n/WORDSIZE).                                  */

        m = SETWORDSNEEDED(n);

     /* The following optional call verifies that we are linking
        to compatible versions of the nauty routines.            */

        nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

     /* Now we create the cycle.  First we zero the graph, than for
        each v, we add the edge (v,v+1), where values are mod n. */

        EMPTYGRAPH(g,m,n);
        for (v = 0; v < n; ++v)  ADDONEEDGE(g,v,(v+1)%n,m);

        printf("Generators for Aut(C[%d]):\n",n);

     /* Since we are not requiring a canonical labelling, the last
        parameter to densenauty() is not required and can be NULL. */

        densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);

     /* The size of the group is returned in stats.grpsize1 and
        stats.grpsize2. */

        printf("Automorphism group size = ");
        writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
        printf("\n");
    }

    exit(0);
}
