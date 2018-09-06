/* This program demonstrates how an isomorphism is found between
   two graphs, using the Moebius graph as an example.
   This version uses Traces and demonstrates how to compute the
   automorphism group separately before computing the canonical
   labelling.  Although this is slower for easy graphs like
   those here, it can be faster for some very difficult graphs.
*/

#include "traces.h"

int
main(int argc, char *argv[])
{
    DYNALLSTAT(int,lab1,lab1_sz);
    DYNALLSTAT(int,lab2,lab2_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,map,map_sz);
    static DEFAULTOPTIONS_TRACES(options);
    TracesStats stats;
    permnode *generators;
 /* Declare and initialize sparse graph structures */
    SG_DECL(sg1); SG_DECL(sg2);
    SG_DECL(cg1); SG_DECL(cg2);

    int n,m,i;

 /* Read a number of vertices and process */

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

            SG_ALLOC(sg1,n,3*n,"malloc");
            sg1.nv = n;              /* Number of vertices */
            sg1.nde = 3*n;           /* Number of directed edges */

            for (i = 0; i < n; ++i)
            {
                sg1.v[i] = 3*i;     /* Position of vertex i in v array */
                sg1.d[i] = 3;       /* Degree of vertex i */
            }
             
            for (i = 0; i < n; i += 2)   /* Spokes */
            {
                sg1.e[sg1.v[i]] = i+1;
                sg1.e[sg1.v[i+1]] = i;
            }

            for (i = 0; i < n-2; ++i)  /* Clockwise edges */
                sg1.e[sg1.v[i]+1] = i+2;
            sg1.e[sg1.v[n-2]+1] = 1;
            sg1.e[sg1.v[n-1]+1] = 0;

            for (i = 2; i < n; ++i)  /* Anticlockwise edges */
                sg1.e[sg1.v[i]+2] = i-2;
            sg1.e[sg1.v[1]+2] = n-2;
            sg1.e[sg1.v[0]+2] = n-1;
                
         /* Now make the second graph */

            SG_ALLOC(sg2,n,3*n,"malloc");
            sg2.nv = n;              /* Number of vertices */
            sg2.nde = 3*n;           /* Number of directed edges */

            for (i = 0; i < n; ++i)
            {
                sg2.v[i] = 3*i;
                sg2.d[i] = 3;
            }

            for (i = 0; i < n; ++i)
            {
                sg2.v[i] = 3*i;
                sg2.d[i] = 3;
                sg2.e[sg2.v[i]] = (i+1) % n;      /* Clockwise */
                sg2.e[sg2.v[i]+1] = (i+n-1) % n;  /* Anti-clockwise */
                sg2.e[sg2.v[i]+2] = (i+n/2) % n;  /* Diagonals */
            }

         /* Now we make the canonically labelled graphs by a two-step
            process.  The first call to Traces computes the
            automorphism group.  The second call computes the
            canonical labelling, using the automorphism group from
            the first call.

            We have declared a variable "generators" that will be
            used to hold the group generators between the two calls.
            It has to be initialised to NULL and its address has to 
            be given to Traces using options.generators.  After the
            second call, we need to discard the generators with a
            call to freeschreier(), which also initializes it again. */
            
            generators = NULL;
            options.generators = &generators;

            options.getcanon = FALSE;
            Traces(&sg1,lab1,ptn,orbits,&options,&stats,NULL);
            options.getcanon = TRUE;
            Traces(&sg1,lab1,ptn,orbits,&options,&stats,&cg1);
            freeschreier(NULL,&generators);

            options.getcanon = FALSE;
            Traces(&sg2,lab1,ptn,orbits,&options,&stats,NULL);
            options.getcanon = TRUE;
            Traces(&sg2,lab1,ptn,orbits,&options,&stats,&cg2);
            freeschreier(NULL,&generators);

         /* Compare canonically labelled graphs */

            if (aresame_sg(&cg1,&cg2))
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
