/* gutils.h - procedure declarations for gutil1.c and gutil2.c */

#ifdef __cplusplus
extern "C" {
#endif

extern void degstats(graph*,int,int,
                     unsigned long*,int*,int*,int*,int*,boolean*);
extern void degstats2(graph*,boolean,int,int,unsigned long*,int*l,
     int*,int*,int*,int*, int*, int*,int*,int*, boolean*);
extern void diamstats(graph*,int,int,int*,int*);
extern void find_dist(graph*,int,int,int,int*);
extern void find_dist2(graph*,int,int,int,int,int*);
extern int girth(graph*,int,int);
extern boolean isbiconnected1(graph*,int);
extern boolean isbiconnected(graph*,int,int);
extern boolean isbipartite(graph*,int,int);
extern int bipartiteside(graph*,int,int);
extern boolean twocolouring(graph*,int*,int,int);
extern boolean isconnected1(graph*,int);
extern boolean isconnected(graph*,int,int);
extern boolean issubconnected(graph*,set*,int,int); 
extern long maxcliques(graph*,int,int);

extern int loopcount(graph*,int,int);
extern long pathcount1(graph*,int,setword,setword);
extern long cyclecount1(graph*,int);
extern long cyclecount(graph*,int,int);
extern long indpathcount1(graph*,int,setword,setword);
extern long indcyclecount1(graph*,int);
extern long indcyclecount(graph*,int,int);
extern void commonnbrs(graph*,int*,int*,int*,int*,int,int);
extern void contract1(graph*,graph*,int,int,int);
extern int cstarcontent(graph*,int,int);
extern long numtriangles1(graph*,int);
extern long numtriangles(graph*,int,int);
extern long numdirtriangles(graph*,int,int);
extern void delete1(graph*,graph*,int,int);
extern int conncontent(graph*,int,int);
extern boolean stronglyconnected(graph*,int,int);

/* extern int diameter_sg(sparsegraph*,int*,int*); */

#ifdef __cplusplus
}
#endif
