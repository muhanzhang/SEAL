/* This program demonstrates how known automorphisms can be given
   to Traces.  We compute the automorphism group of the circulant
   graph of order n with i is adjacent to j iff j-i is a square
   mod n.  We need that -1 is a square so that the graph is
   undirected, which means that the prime factors of n must be 
   congruent to 1 mod 4.  (This is the Paley graph in the event
   that p is a prime.)
*/

#include "traces.h"

int
main(int argc, char *argv[])
{
    DYNALLSTAT(int,p,p_sz);
    DYNALLSTAT(boolean,issquare,issquare_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_TRACES(options);
    TracesStats stats;
 /* Declare and initialize sparse graph structures */
    SG_DECL(sg);

    int deg,n,m,i,j;
    size_t k;
    permnode *gens;

 /* Select option for passing generators to Traces */

    options.generators = &gens;
 
 /* Read a number of vertices and process it */

    while (1)
    {
        printf("\nenter n : ");
        if (scanf("%d",&n) == 1 && n > 2)
        {
            m = SETWORDSNEEDED(n);
            nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

            DYNALLOC1(int,lab,lab_sz,n,"malloc");
            DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
            DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
            DYNALLOC1(int,p,p_sz,n,"malloc");
            DYNALLOC1(boolean,issquare,issquare_sz,n,"malloc");

         /* Initialise list of automorphisms */

            gens = NULL;

         /* Find the squares and the degree */

            for (i = 0; i < n; ++i) issquare[i] = FALSE;
            for (i = 0; i < n; ++i) issquare[(i*i)%n] = TRUE;
            if (!issquare[n-1])
            {
                printf("-1 must be a square mod n; try again\n");
                continue;
            }

            deg = 0;
            for (i = 1; i < n; ++i) if (issquare[i]) ++deg;

         /* Now make the graph */

            SG_ALLOC(sg,n,n*deg,"malloc");
            sg.nv = n;              /* Number of vertices */
            sg.nde = n*deg;           /* Number of directed edges */

            for (i = 0; i < n; ++i)
            {
                sg.v[i] = i*deg;     /* Position of vertex i in v array */
                sg.d[i] = deg;       /* Degree of vertex i */
            }
             
            for (i = 0; i < n; ++i)   /* Edges */
            {
                k = sg.v[i];
                for (j = 1; j < n; ++j)
                    if (issquare[j]) sg.e[k++] = (i + j) % n;
            }

         /* Add known automorphisms */

            /* We wouldn't need freeschreier() if we were only
               processing one graph, but it doesn't hurt.  This
               is how to properly dispose of previous generators. */

            freeschreier(NULL,&gens);

            /* Cyclic rotation */
            for (i = 0; i < n; ++i) p[i] = (i + 1) % n;
            addpermutation(&gens,p,n);

            /* Reflection about 0 */
            for (i = 0; i < n; ++i) p[i] = (n - i) % n;
            addpermutation(&gens,p,n);

         /* Call Traces */
            
            Traces(&sg,lab,ptn,orbits,&options,&stats,NULL);

            printf("Automorphism group size = ");
            writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
            printf("\n");

        /* Traces left the automorphims we gave it, augmented by
           any extra automorphims it found, in a circular list
           pointed to by gens.  See schreier.txt for documentation. */
        }
        else
            break;
    }

    exit(0);
}
