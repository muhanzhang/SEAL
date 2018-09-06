#include <stdio.h>
#include <stdlib.h>
/* Reads a graph in Bliss format from stdin and writes it in dreadnaut
 * format to stdout.  If there is an argument, it is written to the
 * output before the graph.  If there are two arguments, they are
 * written to the output before and after the graph.
 */

typedef struct 
{
   int v,w;
} vpair;

int comparedge(const void *p1, const void *p2)
{
    vpair *e1,*e2;

    e1 = (vpair*)p1;
    e2 = (vpair*)p2;
    if      (e1->v < e2->v) return -1;
    else if (e1->v > e2->v) return 1;
    else if (e1->w < e2->w) return -1;
    else if (e1->w > e2->w) return 1;
    else                    return 0;
}

static int
nextchar(void)
{
    char s[2];

    if (scanf("%1s",s) != 1) return EOF;
    else                     return s[0];
}

int
main(int argc, char *argv[])
{
    int n,c;
    unsigned long ne,j;
    vpair *elist,*vlist;
    int haven;
    int i,v,w;
    int haveptn;

    haven = 0;
    j = 0;
    while ((c = nextchar()) >= 0)
    {
	switch (c)
	{
	case 'c':
	    putchar('"');
	    while ((c = getchar()) != '\n' && c != EOF) putchar(c);
	    printf("\\n\"\n");
	    break;

	case 'p':
	    if (haven)
	    {
		fprintf(stderr,"Duplicate p line\n");
		exit(1);
	    }
	    if (scanf(" edge %d %lu",&n,&ne) != 2)
	    {
		fprintf(stderr,"Bad p line\n");
		exit(1);
	    }
	    if ((elist = (vpair*)malloc(ne*sizeof(vpair))) == NULL
               || (vlist = (vpair*)malloc(n*sizeof(vpair))) == NULL)
	    {
		fprintf(stderr,"Malloc failed\n");
		exit(1);
	    }
	    haven = 1;
	    for (i = 0; i < n; ++i)
	    {
		vlist[i].v = 0;  /* default colour */
		vlist[i].w = i+1;
	    }
	    break;

	case 'n':
	    if (!haven)
	    {
                fprintf(stderr,"Missing p line\n");
                exit(1);
            }  
            if (scanf("%d%d",&w,&v) != 2 || w < 1 || w > n)
            {
                fprintf(stderr,"Bad n line\n");
                exit(1);
            }
	    vlist[w-1].v = v;
	    break;

	case 'e':
	    if (!haven || j == ne)
	    {
		fprintf(stderr,"Missing p line or too many e lines\n");
		exit(1);
	    }
	    if (scanf("%d%d",&v,&w) != 2 || v < 1 || w < 1 || v > n || w > n)
	    {
		fprintf(stderr,"Bad e line\n");
		exit(1);
	    }
	    elist[j].v = v; elist[j].w = w;
	    ++j;
	    break;

	default:
	    fprintf(stderr,"Unknown line %c\n",c);
	    exit(1);
	}
    }

    if (j != ne)
    {
        fprintf(stderr,"Wrong number of e lines\n");
        exit(1);
    }

    if (argc > 1) printf("%s\n",argv[1]);

    printf("$=1 n=%d g\n",n);
    qsort(elist,ne,sizeof(vpair),comparedge);
    
    v = -1;
    for (j = 0; j < ne; ++j)
    {
        if (elist[j].v != v)
        {
    	    v = elist[j].v;
    	    if (j > 0) printf("\n");
	    printf("%d:",v);
	}
	printf(" %d",elist[j].w);
    }
    printf(".\n");

    qsort(vlist,n,sizeof(vpair),comparedge);
    for (i = 1; i < n; ++i)
	if (vlist[i].v != vlist[i-1].v) break;
    if (i < n)
    {
	printf("f=[");
	for (i = 0; i < n; ++i)
	{
	   if (i > 0 && vlist[i].v != vlist[i-1].v)
		printf("\n |");
	   printf(" %d",vlist[i].w);
	}
	printf("]\n");
    }

    printf("$$\n");

    if (argc > 2) printf("%s\n",argv[2]);
    return 0;
}
