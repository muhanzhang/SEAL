/* This program demonstrates how an isomorphism is found between
   two graphs, using the Moebius graphs as an example.
   This version uses dense form with dynamic allocation.
*/

#include "nauty.h" 

int
main(int argc, char *argv[])
{
    DYNALLSTAT(int,lab1,lab1_sz);
    DYNALLSTAT(int,lab2,lab2_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,map,map_sz);
    DYNALLSTAT(graph,g1,g1_sz);
    DYNALLSTAT(graph,g2,g2_sz);
    DYNALLSTAT(graph,cg1,cg1_sz);
    DYNALLSTAT(graph,cg2,cg2_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    int n,m,i;
    size_t k;

 /* Select option for canonical labelling */

    options.getcanon = TRUE;
 
    while (1)
    {
        printf("\nenter n : ");
        if (scanf("%d",&n) == 1 && n > 0)
        {
            if (n%2 != 0)
            {
                fprintf(stderr,"Sorry, n must be even\n");
                continue;
            }

            m = SETWORDSNEEDED(n);
            nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

            DYNALLOC1(int,lab1,lab1_sz,n,"malloc");
            DYNALLOC1(int,lab2,lab2_sz,n,"malloc");
            DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
            DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
            DYNALLOC1(int,map,map_sz,n,"malloc");

         /* Now make the first graph */

            DYNALLOC2(graph,g1,g1_sz,n,m,"malloc");
            EMPTYGRAPH(g1,m,n);

            for (i = 0; i < n; i += 2)   /* Spokes */
                ADDONEEDGE(g1,i,i+1,m);

            for (i = 0; i < n-2; ++i)  /* Cycle */
                ADDONEEDGE(g1,i,i+2,m);
            ADDONEEDGE(g1,1,n-2,m);
            ADDONEEDGE(g1,0,n-1,m);

         /* Now make the second graph */

            DYNALLOC2(graph,g2,g2_sz,n,m,"malloc");
            EMPTYGRAPH(g2,m,n);

            for (i = 0; i < n; ++i)
            {
                ADDONEEDGE(g2,i,(i+1)%n,m);   /* Rim */
                ADDONEEDGE(g2,i,(i+n/2)%n,m); /* Diagonals */
            }

         /* Create canonical graphs */

            DYNALLOC2(graph,cg1,cg1_sz,n,m,"malloc");
            DYNALLOC2(graph,cg2,cg2_sz,n,m,"malloc");

            densenauty(g1,lab1,ptn,orbits,&options,&stats,m,n,cg1);
            densenauty(g2,lab2,ptn,orbits,&options,&stats,m,n,cg2);

         /* Compare canonically labelled graphs */

            for (k = 0; k < m*(size_t)n; ++k)
                if (cg1[k] != cg2[k]) break;

            if (k == m*(size_t)n)
            {
                printf("Isomorphic.\n");
                if (n <= 1000)
                {
                 /* Write the isomorphism.  For each i, vertex lab1[i]
                    of sg1 maps onto vertex lab2[i] of sg2.  We compute
                    the map in order of labelling because it looks better. */

                    for (i = 0; i < n; ++i) map[lab1[i]] = lab2[i];
                    for (i = 0; i < n; ++i) printf(" %d-%d",i,map[i]);
                    printf("\n");
                }
            }
            else
                printf("Not isomorphic.\n");
        }
        else
            break;
    }

    exit(0);
}
