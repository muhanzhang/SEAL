/* This program prints the entire automorphism group of an n-vertex
   polygon, where n is a number supplied by the user. 
*/

#include "nauty.h"    /* which includes <stdio.h> */
#include "naugroup.h"

/**************************************************************************/

void
writeautom(int *p, int n)
/* Called by allgroup.  Just writes the permutation p. */
{
    int i;

    for (i = 0; i < n; ++i) printf(" %2d",p[i]); printf("\n");
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    int n,m,v;
    grouprec *group;

 /* The following cause nauty to call two procedures which
        store the group information as nauty runs. */
        
    options.userautomproc = groupautomproc;
    options.userlevelproc = grouplevelproc;

    while (1)
    {
        printf("\nenter n : ");
        if (scanf("%d",&n) == 1 && n > 0)
        {
            m = SETWORDSNEEDED(n);
            nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

            DYNALLOC2(graph,g,g_sz,m,n,"malloc");
            DYNALLOC1(int,lab,lab_sz,n,"malloc");
            DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
            DYNALLOC1(int,orbits,orbits_sz,n,"malloc");

            EMPTYGRAPH(g,m,n);
            for (v = 0; v < n; ++v) ADDONEEDGE(g,v,(v+1)%n,m);

            printf("Automorphisms of C[%d]:\n",n);
            densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);

         /* Get a pointer to the structure in which the group information
            has been stored.  If you use TRUE as an argument, the
            structure will be "cut loose" so that it won't be used
            again the next time nauty() is called.  Otherwise, as
            here, the same structure is used repeatedly. */
                
            group = groupptr(FALSE);

         /* Expand the group structure to include a full set of coset
            representatives at every level.  This step is necessary
            if allgroup() is to be called. */
                
            makecosetreps(group);

         /* Call the procedure writeautom() for every element of the group.
            The first call is always for the identity. */
                
            allgroup(group,writeautom);
        }
        else
            break;
    }
    exit(0);
}
