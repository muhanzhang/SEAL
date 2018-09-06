// cc -O4 -o water2 -DWORDSIZE=32 -DMAXN=WORDSIZE nauty.c naugraph.c nautil.c gtools.c schreier.c naurng.c watercluster2.c

/*
Reads graphs in g6 code or multicode (optional) from stdin and directs them 

options: 

ix  means: the indegree of every vertex may be at most x.

oy  means: the outdegree of every vertex may be at most y.

  S  means: allow that for every pair of vertices x,y at most one of the edges x-->y 
     and y-->x may be present. By default both of them may be present in the same graph.


  T  means: Output directed graphs in T-code. This is a simple ASCII output format. Every line
     contains one graph. First the number of vertices, then the number of 
     directed edges and then the list of directed edges with the start first 
     and the end then. E.g.: 3 2 0 1 2 1 means 3 vertices, 2 directed edges:
     0-->1 and 2-->1

  B  means: Output the directed graphs in a binary code. Every item of the code is an unsigned
     char. The first unsigned char is the number nv of vertices. The vertices are numbered 1..nv
     Then the list of vertices x for which there is a directed edge 1->x follow. This list is
     ended by a 0. Then the list of outgoing neighbours of 2 follows -- again ended with a 0, etc.
     The code is complete with the 0 ending the list of outgoing neighbours of nv.

  C  means: Do really construct all the directed graphs in memory, but don't output them. This is not
     a big difference in case of restricted in- and outdegrees, because all that is done extra is that 
     edges are directed instead of just keeping track of in- and out-degrees. This option is intended only
     for testing purposes to test also routines that are normally not used when counting. Things that would 
     speed up the counting also in some cases of restricted in- and out-degrees -- like multiplying the 
     possibilities of assigning directions to edges that can be assigned directions independent 
     of each other (depending on the degrees of the endvertices and overlaps) -- are not included. 
     In case of not restrictive bounds on the in- and out-degree it not really constructing the graphs
     can be considerably faster. In cases of restricted in- and out-degrees the only difference is that
     the graph isn't modified...
     The fact that in case of no output the graph is not modified is mainly to save time for the one 
     case of waterclusters, where large numbers were determined. If large numbers (without output)
     for other cases shall be determined, one should think about adding the multiplication routines.

   m read multicode

This program uses different labelling routines -- all based on the ideas of 

G. Brinkmann, Generating water clusters and other directed graphs,
Journal of Mathematical Chemistry 46, 1112--1121 (2009)

October 10, 2011: corrected error caused by overflow of 32bit int used as hashvalue.

Sep, 2012: PROCESS feature added by BDM.
*/

/* PROCESS feature
 *
 * If PROCESS is defined, it must expand as the name of a procedure
 * with prototype like void PROCESS(FILE *f, graph *g, int n). This
 * procedure will be called for each output graph before it is
 * output (or before output is suppressed) with f being the output
 * file (possibly NULL).
 * It is an error if n > WORDSIZE.
 * If NOCONVERSE is also defined, the call is only made for one of
 * each digraph and its converse.
 */

/* SUMMARY feature
 *
 * If SUMMARY is defined, it must expand as the name of a procedure
 * with prototype  void SUMMARY(void).  It is called at the end after
 * the normal summary. 
 */

//#include<stdio.h>
#include "nauty.h"
#include<limits.h>
#include<string.h>

#ifdef PROCESS
extern void PROCESS(FILE*,graph*,int);
nauty_counter dg_nin,dg_nout;
#endif
#ifdef SUMMARY
extern void SUMMARY(void);
#endif

#define MAX_BOGEN ((MAXN*(MAXN-1))/2)
#define INFTY_UCHAR UCHAR_MAX
typedef int BOOG[2];

void nontrivlabels(), init_nauty_options();

BOOG edgelist[MAX_BOGEN+1]; /* de lijst van bogen */
BOOG edgelist_final[MAX_BOGEN+1]; /* de lijst van bogen die nadat er eerst nontriviale automorphismen waren die dan 
				  verdwenen nog gericht moeten worden.*/
BOOG *laatstepositie; /* The last position in one of these lists. Global in order not to have to copy it every time */

/* all the next variables must be evaluated immediately when orbits are constructed. They are reused in the next
   iteration */
unsigned char *operations=NULL;
int *root_op=NULL, size_root=0, blocklength, orbitblocklength[MAXN], size_operations=0, number_operations=0;
/* these arrays will be dynamically allocated and extended. Operations will be an array with "blocks" of length
   "blocklength=tobedirected[vertex_in_orbit]+3". They represent an operation the following way: the first entry is the central
   vertex then the list of vertices from which only edges come in -- ended by INFTY_UCHAR, then the list of vertices 
   to which only outgoing edges go -- ended by INFTY_UCHAR. Then the list of vertices to which incoming AND outgoing 
   edges go -- ended by INFTY_UCHAR. 

   root has as many entries as there are blocks in operations. root[i]=j means that when computing the orbits of
   operations and constructing a union-find-tree, block number j is the root of the tree containing vertex i. */

unsigned char *remember_operations[MAXN]; //the nonequivalent operations that are stored for every iteration
int remember_size[MAXN]; // remember_size[i] is the number of characters allocated for remember_operations[i]

#define COPYBOOG(a,b) { (a)[0]=(b)[0]; (a)[1]=(b)[1]; }


/* OPTIONS */
boolean dummybool;
int mingerichtdeg;
int remaining_doubles=0; /* hoeveel bogen kunnen ten hoogste nog in allebei 
			    richtingen gericht worden ? */
int watermaxedges; /* what is the theoretical maximum for the number of edges which can
		      directed in a way that doesn't give a conflict with the in/out-degree
		      bounds. */
int watermaxdeg; /* maxindeg+maxoutdeg */

int is_gericht[MAX_BOGEN][MAX_BOGEN]={{0}}; /*  is_gericht[i][j]==1 als {i,j} is al gericht -- anders 0 */
int virtual_gericht[MAX_BOGEN][MAX_BOGEN]={{0}}; /*  virtual_gericht[i][j]==1 als het al vastgelegd is dat {i,j} 
						    i->j gericht moet worden, 2 als hij j->i gericht moet worden
						    en 0 als het nog niet vastgelegd is. Vastgelegd betekent door
						    de beperkingen van in- en outgraad moet deze boog zo gericht 
						    worden. */
int positie[MAXN][MAXN]; /* the value of positie[i][j] is the position of edge {1,j} in edgelist in case of
			    nontrivial automorphisms of the underlying undirected graph. */

int virtual_indeg[MAXN], virtual_outdeg[MAXN]; 
                                               /* de graden als de nog niet gerichte maar al vastgelegde 
						  bogen meegerekend worden -- dat mag niet voor kanoniciteit 
						  gebruikt worden */
int nextstep_depth; 

#define SWITCHPAR_ORBSIZE 5 // when does serial orbit labelling switch to parallel -- and stay there
#define MAXPAR_ORBSIZE 9 // maximum 16, when changing, change following MAXPAROPS too
#define MAXPAROPS 19683 // set to minimum 3^MAXPAR_ORBSIZE
unsigned int parops[MAXPAROPS];
// in the case of no degree restrictions this could in fact be permanently filled in -- only depending
// on the edge orbit size. But the time for filling it in is so small that it is simply not worth it.
int number_parops; 
// an operation is defined as follows: if edge number i in kleinste_orbit must be directed kleinste_orbit[i][0]->kleinste_orbit[i][1],
// bits 2i and 2i+1 form the number 0 (so are both 0). For kleinste_orbit[i][0]<-kleinste_orbit[i][1] they form 2, for
// a double edge they form 1 ---- so 2-type is always the inverse
#define SETOP(op,edge,type) ((op) = ((op) & ~(3<<((edge)<<1))) | (type)<<((edge)<<1))
#define GETTYPE(op,edge) (((op)>>((edge)<<1))&3)

/* for the following macros global variables _x_ and _y_ are used */
int _x_, _y_;
#define BILDBOOG_UNDIR(startboog,bildboog,permnummer) \
  { _x_=generators[permnummer][(startboog)[0]]; _y_=generators[permnummer][(startboog)[1]]; \
if (_x_<_y_) { (bildboog)[0]=_x_; (bildboog)[1]=_y_; } else  { (bildboog)[0]=_y_; (bildboog)[1]=_x_; }}
#define BILDBOOG(startboog,bildboog,permnummer) \
  { (bildboog)[0]=generators[permnummer][(startboog)[0]]; (bildboog)[1]=generators[permnummer][(startboog)[1]]; }
#define POSBILD(startboog,permnummer) \
  (positie[generators[permnummer][(startboog)[0]]][generators[permnummer][(startboog)[1]]])

/* a macro that sets numbers a and b equivalent to a */
#define SETEQUIV(a,b) { _y_=(b); while ((_x_=number[_y_])!=_y_) { number[_y_]=(a); _y_=_x_; } number[_y_]=(a); } 
/* zorgt ervoor dat nummer[a] zo is dat nummer[nummer[a]]=nummer[a] */
#define UPDATENUMBER(a) { _y_=(a); while ((_x_=number[_y_])!=_y_) { number[_y_]=number[_x_]; _y_=number[_y_]; } number[a]=number[_y_]; } 

/* the quality of a directed edge */
#define QUALITY(a,b) (((saturated[a]+saturated[b])<<6)+colour[nextstep_depth][a])
// upper bound for the quality if the in- or out-deg of one of the entries is changed by one
#define QUALITY_P1(a,b) (((saturated[a]+saturated[b]+1)<<6)+colour[nextstep_depth][a])
// quality if the in- or out-deg of one of the entries are both changed by one
#define QUALITY_P2(a,b) (((saturated[a]+saturated[b]+2)<<6)+colour[nextstep_depth][a])

int _marks[MAX_BOGEN], _markvalue=INT_MAX;
#define RESETMARKS { int i; if (_markvalue<INT_MAX) _markvalue++; else { _markvalue=1; \
                     for (i=0;i<MAX_BOGEN;i++) _marks[i]=0; } }
#define ISMARKED(i) (_marks[i]==_markvalue)
#define UNMARKED(i) (_marks[i]!=_markvalue)
#define MARK(i) {_marks[i]=_markvalue;}

long long int onlevel[MAX_BOGEN]={0};

/* variabelen voor nauty: nvector sounds like a vector but is a typedef for int */
nvector lab[MAX_BOGEN][MAXN]={{0}}, ptn[MAX_BOGEN][MAXN]={{0}}, orbits[MAX_BOGEN], colour[MAX_BOGEN][MAXN]={{0}};
/* MAX_BOGEN is zeker te veel -- maar wat zou een goede bovengrens zijn ? */
/* in the arrays lab[i][] and ptn[i][] the partition is stored after having labelled i orbits 
   completely. lab[0][] and ptn[0][] are the partitions as computed for the original graph. */
int rememberorbits[MAXN][MAXN];

nvector bufferlab[MAXN], bufferptn[MAXN];
graph canong[MAXN];
/* make sure these are always evaluated at once, because they are reused whenever new
   edges are directed */
graph workg[MAXN], staticg[MAXN], canong[MAXN];
graph bit_orbit[MAXN];
int deg[MAXN]={0}, outdeg[MAXN]={0}, indeg[MAXN]={0};
int tobedirected[MAXN]={0}; /* how many edges have still to be directed? This can
   also be used to detect whether an edge that is in both directions is a double edge or 
   just not yet directed in the routines for nontrivial symmetry. It is a double edge if 
   and only if one of the endpoints is ready with ready defines as: */
#define READY(i) (tobedirected[i]==0)
#define NOT_READY(i) (tobedirected[i]>0)
int aantal_toppen, aantal_bogen, aantal_gerichte_bogen;
int max_doubles;
int double_free[MAXN]={0}; /* hoeveel dubbele bogen kan top [i] nog krijgen? */

long long int addnumber=0; /* How much must the counter be increased if a graph
			      is found. Just relevant in case not all are really constructed and coded.
			      Then this gives a speedup.*/ 

long long int aantal_grafen_met_triv_group=0LL;
long long int aantal_gerichte_grafen=0LL;

static DEFAULTOPTIONS(options_directed);
static DEFAULTOPTIONS(options_directed_canon);
static DEFAULTOPTIONS(options);
static DEFAULTOPTIONS(options_final);
static statsblk stats;
setword workspace[100*MAXN];

permutation generators[MAXN][MAXN];
int number_of_generators;
int group_up_to_date=0;

int indeg_free[MAXN]={0}, outdeg_free[MAXN]={0}, saturated[MAXN]={0};

graph all; // the set of all vertices

int orbitchoices[MAXN]={0};

#define CSIZE 32768
int _colourmarks[CSIZE], _cmark=INT_MAX;
#define RESETMARKS_COLOUR { int i; \
                            if (_cmark==INT_MAX)\
				     { for (i=0;i<CSIZE;i++) _colourmarks[i]=0; _cmark=1; }\
			    else _cmark++; }
#define MARK_COLOUR(a) {_colourmarks[a]=_cmark;}
#define ISMARKED_COLOUR(a) (_colourmarks[a]==_cmark)

void directorbit();
graph* readg();
void parallel_orbit_labelling();

//unsigned int waterclusteruse=0UL, water_v_use=0UL;

/****OPTIONS******/
int maxindeg=MAXN, maxoutdeg=MAXN, maxdirectdeg=MAXN, nodegbound=1;
int double_allowed=1, direct_output=0;

#define NEXTEL(a,b) (b)=FIRSTBIT((a)& BITMASK(b)) /*(b)++, b=FIRSTBIT((a)<<(b))+(b)*/
#define FORALLELEMENTS(x,y) for (y= FIRSTBIT(x); y<WORDSIZE; NEXTEL((x),(y)))
#define FORALLELEMENTS_BOUND(x,y,b) for (y=FIRSTBIT((x)& BITMASK(b)) ; y<WORDSIZE; NEXTEL((x),(y)))
/* Attention: bound is a lower bound and the first possible value is bound+1 ! */

#define WRITEUP() { aantal_gerichte_grafen+=addnumber; if (direct_output==1) {writeTcode(workg,aantal_toppen);\
if (addnumber==2) writeTcode_invers(workg,aantal_toppen);} else if (direct_output==3) \
{writeBcode(workg,aantal_toppen);\
if (addnumber==2) writeBcode_invers(workg,aantal_toppen);} }
#define WRITEUP_COUNT() { aantal_gerichte_grafen+=addnumber; }

#ifdef PROCESS
 void callprocess(graph *g, int aantal_toppen)
   { int i,j; graph gconv[MAXN];

     dg_nout = aantal_gerichte_grafen+1;
     PROCESS(stdout,g,aantal_toppen);
#ifndef NOCONVERSE
     if (addnumber==2)
       {
	  EMPTYSET(gconv,aantal_toppen);
          for (i=0; i<aantal_toppen;i++)
            { 
	       FORALLELEMENTS(g[i],j) ADDELEMENT(gconv+j,i);
	    }
            ++dg_nout;
          PROCESS(stdout,gconv,aantal_toppen);
       }
#endif
   }
#define MAYBEPROCESS callprocess(workg,aantal_toppen)
#else
#define MAYBEPROCESS
#endif

 void writeTcode(graph *g, int aantal_toppen)
   { int i,j;

   fprintf(stdout,"%d %d",aantal_toppen,aantal_bogen+max_doubles-remaining_doubles);
   for (i=0; i<aantal_toppen;i++) 
     { 
       FORALLELEMENTS(g[i],j) fprintf(stdout," %d %d",i,j);
     }
       fprintf(stdout,"\n");
   }

 void writeTcode_invers(graph *g, int aantal_toppen)
   { int i,j;

   fprintf(stdout,"%d %d",aantal_toppen,aantal_bogen+max_doubles-remaining_doubles);
   for (i=0; i<aantal_toppen;i++) 
     { 
       FORALLELEMENTS(g[i],j) fprintf(stdout," %d %d",j,i);
     }

       fprintf(stdout,"\n");
   }

 void writeBcode(graph *g, int aantal_toppen)
   { int i,j;
   putc(aantal_toppen,stdout);
   for (i=0; i<aantal_toppen;i++) 
     { 
       FORALLELEMENTS(g[i],j) putc(j+1,stdout);
       putc(0,stdout);
     }
   }

 void writeBcode_invers(graph *g, int aantal_toppen)
   { int i,j;
     int lijst[MAXN][MAXN], lengte[MAXN];
     putc(aantal_toppen,stdout);

     for (i=0; i<aantal_toppen;i++) lengte[i]=0;
     for (i=0; i<aantal_toppen;i++) 
       { 
	 FORALLELEMENTS(g[i],j) { lijst[j][lengte[j]]=i+1; lengte[j]++; }
       }
     for (i=0; i<aantal_toppen;i++)
       { for (j=0; j<lengte[i]; j++) putc(lijst[i][j],stdout);
         putc(0,stdout);
       }
   }

void writelab(nvector ptn[], nvector lab[])
{
  int i;

  for (i=0; i<aantal_toppen; i++) 
    { fprintf(stderr," %d ",lab[i]); if (ptn[i]==0) fprintf(stderr,"|| "); }
  fprintf(stderr,"\n");

  return;
}

void writeset(graph set)
{ int i;
  fprintf(stderr,"The set:\n");
  FORALLELEMENTS(set,i) fprintf(stderr,"%d ",i);
  fprintf(stderr,"\n");
}

void writeoperation(unsigned char op[])
{ int i,counter;

  fprintf(stderr,"Operation: ");
  for (i=counter=0; counter<3; i++)
    {
      if (op[i]==INFTY_UCHAR) { fprintf(stderr," ||"); counter++; }
      else fprintf(stderr," %d",op[i]);
    }
  fprintf(stderr,"\n");
}

void writelist(int list[])
{
  int i;

  for (i=0;list[i]>=0;i++) fprintf(stderr,"%d ",list[i]); fprintf(stderr,"\n");
}

 void writegraph(graph *g)
 // for graphs as they are used in the vertexorbit routines
   { int i,j;
     static int counter=0;
   fprintf(stderr,"---------------------------------------------------------\n");
   fprintf(stderr,"Graph number %d with %d vertices\n",++counter,aantal_toppen);
   for (i=0; i<aantal_toppen;i++) 
     { fprintf(stderr,"%d:",i);
       FORALLELEMENTS(g[i],j) 
	 { fprintf(stderr," %d", j);
	   if (READY(i) || READY(j)) 
	     { 
	       if (ISELEMENT(g+j,i)) fprintf(stderr,"[d]"); else fprintf(stderr,"[s]");
	     }
	       else fprintf(stderr,"[ng]");
	 }
       fprintf(stderr,"\n");
     }
   fprintf(stderr,"---------------------------------------------------------\n");
   }


/**********************************************************************************/

void sammle_permutationen(int count, permutation perm[], nvector orbits[],
                          int numorbits, int stabvertex, int n)
{
  memcpy(generators+number_of_generators,perm,sizeof(permutation)*n);

  number_of_generators++;
}


/**********************************************************************************/

void init_nauty_options()
/* initialises the nauty variables */
{
 /* tc_level = 0 is not the default in the most recent
    editions of nauty.  However, it is better for very small
    graphs so we set it here. */

options.getcanon=FALSE;
options.userautomproc = sammle_permutationen;
options.defaultptn=TRUE;
options.writeautoms=FALSE;
options.writemarkers=FALSE;
options.digraph=FALSE;
options.tc_level = 0;


options_directed.getcanon=FALSE;
options_directed.userautomproc = sammle_permutationen;
options_directed.defaultptn=FALSE;
options_directed.writeautoms=FALSE;
options_directed.writemarkers=FALSE;
options_directed.digraph=TRUE;
options_directed.tc_level = 0;

options_directed_canon.getcanon=TRUE;
options_directed_canon.userautomproc = sammle_permutationen;
options_directed_canon.defaultptn=FALSE;
options_directed_canon.writeautoms=FALSE;
options_directed_canon.writemarkers=FALSE;
options_directed_canon.digraph=TRUE;
options_directed_canon.tc_level = 0;


options_final.getcanon=TRUE;
options_final.defaultptn=TRUE;
options_final.writeautoms=FALSE;
options_final.writemarkers=FALSE;
options_final.digraph=TRUE;
options_final.tc_level = 0;
}

/****************************LESE_MULTICODE************************/

int lese_multicode(unsigned char**code, int *codelaenge, FILE *fil)

/* Liest den code und gibt EOF zurueck, wenn das Ende der Datei erreicht
   ist, 1 sonst. Der Speicher fuer den code wird alloziert, falls noetig,
   was entschieden wird anhand der lokalen Variablen maxknotenzahl */
{
static int maxknotenzahl= -1;
int codel, gepuffert=0;
int knotenzahl,a=0,b=0,nuller;
unsigned char ucharpuffer;

if ((knotenzahl=getc(fil))==EOF) return EOF;
if (knotenzahl==0) { fprintf(stderr,"Umschaltung auf short noch nicht implementiert.\n");
                     exit(0);
                    } 
nuller=0; codel=1;
if (knotenzahl=='>') /* koennte ein header sein -- oder 'ne 62, also ausreichend fuer
			     unsigned char */
      { gepuffert=1;
	a=getc(fil);
	if(a==0) nuller++; 
	b=getc(fil);
	if(b==0) nuller++; 
	/* jetzt wurden 3 Zeichen gelesen */
	if ((a=='>') && (b=='m')) /*garantiert header*/
	  { while ((ucharpuffer=getc(fil)) != '<') {}
	    /* noch zweimal: */ ucharpuffer=getc(fil); 
	    if (ucharpuffer!='<') { fprintf(stderr,"Problems with header -- single '<'\n"); exit(1); }
	    if ((knotenzahl=getc(fil))==EOF) return EOF;
	    /* kein graph drin */
	  }
	/* else kein header */
      }


if (knotenzahl > maxknotenzahl)
  { if (*code) free(*code);
    *code=(unsigned char *)malloc((knotenzahl*(knotenzahl-1)/2+knotenzahl)*sizeof(unsigned char));
    if (code==NULL) { fprintf(stderr,"Do not get memory for code\n"); exit(0); }
    maxknotenzahl=knotenzahl;
  }

(*code)[0]=knotenzahl; if (gepuffert) { codel=3; (*code)[1]=a; (*code)[2]=b; }

while (nuller<knotenzahl-1)
  { if (((*code)[codel]=getc(fil))==0) nuller++;
    codel++; }

*codelaenge=codel;
return 1;
}


void usage(char name[])
{

  fprintf(stderr,"usage: %s [ix] [oy] [m] [T] [C] [B] [S] .\n",name);
  fprintf(stderr,"The option ix restricts the maximum indegree to x.\n");
  fprintf(stderr,"The option oy restricts the maximum outdegree to y.\n");
  fprintf(stderr,"The default maximum in- and out-degrees are unlimited. \n");
  fprintf(stderr,"T means: Output directed graphs in T-code -- for details see header\n");
  fprintf(stderr,"B means: Output directed graphs in binary code -- for details see header\n");
  fprintf(stderr,"C means: Do really construct all the directed graphs in memory, but don't output them.\n");
  fprintf(stderr,"S means that for each edge only one direction must be chosen -- not both.\n");
  fprintf(stderr,"Default is that both are allowed\n");
  fprintf(stderr,"  -- so the edge a-b can become a->b AND b->a in the same output graph.\n");
  fprintf(stderr,"m means: read multicode instead of g6 code \n");
  exit(1);

}

/**********DECODE_TO_NAUTY****************************************************/

void decode_to_nauty(unsigned char *code, int codelaenge, graph *g, int degree[]) 


  /* Dekodiert multicode nach nauty bitcode */

  /* alle knotennamen muessen fuer nauty um 1 nach unten verschoben werden */

{

  int v1,v2;
  unsigned char *end;

 if (code[0]>32) { fprintf(stderr,"Not prepared for %d>32 vertices.\n",code[0]); exit(2); }
 aantal_toppen=code[0];
 aantal_bogen=codelaenge-aantal_toppen;
 for (v1=0; v1<aantal_toppen; v1++) { EMPTYSET1(g+v1,1); degree[v1]=0; }
 for (v1=0, end=code+codelaenge, code++; code<end; code++) 
   { if (*code==0) v1++;
     else { v2=(*code)-1;
            ADDELEMENT1(g+v1,v2); ADDELEMENT1(g+v2,v1); 
	    degree[v1]++; degree[v2]++;
          }
   }
 aantal_gerichte_bogen=0;
 //for (v1=0; v1<aantal_toppen; v1++) tobedirected[v1]=degree[v1];

 return;
}

/****************************************INIT_FOR_G6*****************************/

void init_for_g6(graph g[],int aantal_toppen, int degree[])
{
  int i;

  if (aantal_toppen>32) { fprintf(stderr,"Not prepared for %d>32 vertices.\n",aantal_toppen); exit(2); }
  aantal_bogen=0;

  for (i=0; i<aantal_toppen; i++) { degree[i]=POPCOUNT(g[i]); aantal_bogen+=degree[i]; }
  aantal_bogen = aantal_bogen>>1;
  aantal_gerichte_bogen=0;
  return;

}


/***************************FILL_EDGELIST**************************/

void fill_edgelist()
     /* writes the edges in g into the list in a lexicographic way. 
	Is used if no bounds for in- and out-degree are given
	or few edges are left. */

{
  int i,j,end;

  for (i=j=0; i<aantal_toppen; i++) /* j is de positie in edgelist */
    { indeg_free[i]=maxindeg-indeg[i]; outdeg_free[i]=maxoutdeg-outdeg[i];
      FORALLELEMENTS_BOUND(workg[i],end,i)
	{ 
	  edgelist[j][0]=i; edgelist[j][1]=end; j++;
	}
    }

}

/***************************FILL_EDGELIST_ORDER**************************/


void fill_edgelist_order()
  /* writes the edges in gg into the list in a way that for i in 1...numberedges we always have
     that the starting point of edge i has the minimum degree of all vertices in the graph you get
     when removing all edges i+1...aantal_bogen from gg.
     The global variable positie is not assigned any values.
 */
{
  int last_positie,i, buffer, olddeg, newdeg, beste, start, end;
  int list[MAXN][MAXN], listlen[MAXN]={0}; /* list[i] contains the vertices with degree i at that moment */
  int toppositie[MAXN], bufferdeg[MAXN]; /* toppositie[i] is de positie van top i in de lijst */
  int buren[MAXN]; /* de buren van de top waaraan gewerkt wordt */
  graph dummy[MAXN];

  if (nodegbound) { fill_edgelist(); return; }

  memcpy(dummy,workg,aantal_toppen*sizeof(graph));

    for (i=0; i<aantal_toppen; i++) 
      { buffer=bufferdeg[i]=deg[i]; list[buffer][listlen[buffer]]=i; 
	toppositie[i]=listlen[buffer]; (listlen[buffer])++; 
	indeg_free[i]=maxindeg-indeg[i]; outdeg_free[i]=maxoutdeg-outdeg[i];}


    last_positie=aantal_bogen-1;
    while (last_positie>=0)
      {
	for (beste=1;listlen[beste]==0; beste++);
	(listlen[beste])--; /* zo wordt hij ook verwijderd */
	start=list[beste][listlen[beste]];
	for (i=0; i<beste; i++) /* alle buren opslaan*/
	  { end=buren[i]=FIRSTBIT(dummy[start]); DELELEMENT(dummy+start,end); }

	for (i=0; i<beste; i++) /* alle buren */
	  { 
	    end=buren[i];
	    edgelist[last_positie][0]=start; edgelist[last_positie][1]=end; 
	    last_positie--;
	    DELELEMENT(dummy+end,start);
	    /* de buur verhuizen: */
	    olddeg=bufferdeg[end];
	    (bufferdeg[end])--;
	    newdeg=bufferdeg[end];
	    /* uit de oude lijst verwijderen */
	    if (listlen[olddeg]==1) listlen[olddeg]=0;
	    else { (listlen[olddeg])--;
	    buffer= list[olddeg][listlen[olddeg]];
	    list[olddeg][toppositie[end]]=buffer;
	    toppositie[buffer]=toppositie[end]; }
	    /* tot de nieuwe toevoegen */
	    if (newdeg)
	      { list[newdeg][listlen[newdeg]]=end;
	        toppositie[end]=listlen[newdeg];
	        (listlen[newdeg])++;
	      }
	  }

      }

}

/**********************************TRIVLABEL_ROUTINES****************************/

 void trivlabels_nowrite_nodouble(BOOG *positie)
   {
     int start, end;

     start=(*positie)[0]; end=(*positie)[1];

  if (positie==laatstepositie)
    {
     if (indeg_free[start] && outdeg_free[end]) /* bogen kan end->start gericht worden */
       {
	 WRITEUP_COUNT();
       }
    if (indeg_free[end] && outdeg_free[start]) /* bogen kan start->end gericht worden */
       { 
	 WRITEUP_COUNT();
       }

    return;
    }

  /* else -- niet op laatste positie */

     if (indeg_free[start] && outdeg_free[end]) /* bogen kan end->start gericht worden */
       { 
	 (indeg_free[start])--; (outdeg_free[end])--;
	 trivlabels_nowrite_nodouble(positie+1);
	 (indeg_free[start])++; (outdeg_free[end])++;
       }
    if (indeg_free[end] && outdeg_free[start]) /* bogen kan start->end gericht worden */
       {
	 (indeg_free[end])--; (outdeg_free[start])--;
	 trivlabels_nowrite_nodouble(positie+1);
	 (indeg_free[end])++; (outdeg_free[start])++;
       }
    return;
   }




 void trivlabels_nowrite(BOOG *positie)
   {
     int start, end, counter;

     if (!remaining_doubles) { trivlabels_nowrite_nodouble(positie); return; }

     start=(*positie)[0]; end=(*positie)[1];

  if (positie==laatstepositie)
    {
     if (indeg_free[start] && outdeg_free[end]) /* bogen kan end->start gericht worden */
       { counter=1;
	 WRITEUP_COUNT();
       }
     else counter=0;
    if (indeg_free[end] && outdeg_free[start]) /* bogen kan start->end gericht worden */
       { counter++;
	 WRITEUP_COUNT();
       }
    if (remaining_doubles && (counter==2)) WRITEUP_COUNT();

    return;
    }

  /* else -- niet op laatste positie */

     if (indeg_free[start] && outdeg_free[end]) /* bogen kan end->start gericht worden */
       { counter=1;
	 (indeg_free[start])--; (outdeg_free[end])--;
	 trivlabels_nowrite(positie+1);
	 (indeg_free[start])++; (outdeg_free[end])++;
       }
     else counter=0;
    if (indeg_free[end] && outdeg_free[start]) /* bogen kan start->end gericht worden */
       { counter++;
	 (indeg_free[end])--; (outdeg_free[start])--;
	 trivlabels_nowrite(positie+1);
	 (indeg_free[end])++; (outdeg_free[start])++;
       }
    if (remaining_doubles && (counter==2) && double_free[start] && double_free[end])
      {
	 (indeg_free[end])--; (outdeg_free[start])--;
	 (indeg_free[start])--; (outdeg_free[end])--;
	 double_free[start]--; double_free[end]--;
	 remaining_doubles--;
	 trivlabels_nowrite(positie+1);
	 remaining_doubles++;
	 double_free[start]++; double_free[end]++;
	 (indeg_free[end])++; (outdeg_free[start])++;
	 (indeg_free[start])++; (outdeg_free[end])++;
      }
    return;
   }


 void trivlabels(BOOG *positie)
   {
     int start, end, counter;

  start=(*positie)[0]; end=(*positie)[1];

  

  if (positie==laatstepositie)
    {
     if (indeg_free[start] && outdeg_free[end]) /* bogen kan end->start gericht worden */
       { counter=1;
	 //(indeg_free[start])--; (outdeg_free[end])--;
	 DELELEMENT(workg+start,end);
	 MAYBEPROCESS;
	 WRITEUP();
	 ADDELEMENT(workg+start,end);
	 //(indeg_free[start])++; (outdeg_free[end])++;
       }
     else counter=0;
    if (indeg_free[end] && outdeg_free[start]) /* bogen kan start->end gericht worden */
      {  counter++;
	//(indeg_free[end])--; (outdeg_free[start])--;
	 DELELEMENT(workg+end,start);
	 MAYBEPROCESS;
	 WRITEUP();
	 ADDELEMENT(workg+end,start);
	 //(indeg_free[end])++; (outdeg_free[start])++;
       }
    if (remaining_doubles && (counter==2)) { remaining_doubles--; MAYBEPROCESS; WRITEUP(); remaining_doubles++; }
    return;
    }

  /* else -- niet op laatste positie */

     if (indeg_free[start] && outdeg_free[end]) /* bogen kan end->start gericht worden */
       { counter=1;
	 (indeg_free[start])--; (outdeg_free[end])--;
	 DELELEMENT(workg+start,end);
	 trivlabels(positie+1);
	 ADDELEMENT(workg+start,end);
	 (indeg_free[start])++; (outdeg_free[end])++;
       }
     else counter=0;
    if (indeg_free[end] && outdeg_free[start]) /* bogen kan start->end gericht worden */
      {  counter++;
	 (indeg_free[end])--; (outdeg_free[start])--;
	 DELELEMENT(workg+end,start);
	 trivlabels(positie+1);
	 ADDELEMENT(workg+end,start);
	 (indeg_free[end])++; (outdeg_free[start])++;
       }
    if (remaining_doubles && (counter==2) && double_free[start] && double_free[end])
      {
	 (indeg_free[end])--; (outdeg_free[start])--;
	 (indeg_free[start])--; (outdeg_free[end])--;
	 double_free[start]--; double_free[end]--;
	 remaining_doubles--;
	 trivlabels(positie+1);
	 remaining_doubles++;
	 double_free[start]++; double_free[end]++;
	 (indeg_free[end])++; (outdeg_free[start])++;
	 (indeg_free[start])++; (outdeg_free[end])++;
      }

    return;
   }


 void trivlabels_init(BOOG *positie)
   {
     int remember;

     if (!direct_output) 
       { if (nodegbound)
	 { remember=addnumber;
	   if (remaining_doubles) { for (;positie<=laatstepositie;positie++) addnumber*=3;}
	   else { for (;positie<=laatstepositie;positie++) addnumber*=2;}
	   WRITEUP();
	   addnumber=remember;
	   return;
	   }
       /* else */
         if (!remaining_doubles) trivlabels_nowrite_nodouble(positie);
	   else trivlabels_nowrite(positie); 
	 return; 
       }
     /* else */
     trivlabels(positie); 
     return;
   }



/******************************DIRECT_ALL_TRIV********************************/

void direct_all_triv()
{
  int i, start, end;

  aantal_grafen_met_triv_group++; 
  if (nodegbound) fill_edgelist(); else fill_edgelist_order();
  for (i=0;i<aantal_toppen;i++) double_free[i]=maxdirectdeg-deg[i];
  laatstepositie=edgelist+aantal_bogen-1;
  if (maxoutdeg==maxindeg) /* then every valid graph for one direction is also valid with
			      the directions reversed */
    { 
      addnumber=2;
      start=edgelist[0][0]; end=edgelist[0][1];
      /* is_gericht[][] isn't used in the directing routine */
      (outdeg_free[start])--; (indeg_free[end])--;
      //virtual_outdeg[start]++; virtual_indeg[end]++;
      DELELEMENT(workg+end,start);
      aantal_gerichte_bogen=1;
      trivlabels_init(edgelist+1);
      aantal_gerichte_bogen=0;
      ADDELEMENT(workg+end,start);
      (outdeg_free[start])++; (indeg_free[end])++;
      
      if (remaining_doubles && double_free[start] && double_free[end])
	{
	  addnumber=1; 
	  (outdeg_free[start])--; (indeg_free[end])--;
	  (outdeg_free[end])--; (indeg_free[start])--;
	  double_free[start]--; double_free[end]--;
	  aantal_gerichte_bogen=1;
	  remaining_doubles--;
	  trivlabels_init(edgelist+1);
	  remaining_doubles++;
	  double_free[start]++; double_free[end]++;
	  (outdeg_free[start])++; (indeg_free[end])++;
	  (outdeg_free[end])++; (indeg_free[start])++;
	  aantal_gerichte_bogen=0;
	}
    }
  else
    {
      addnumber=1;
      aantal_gerichte_bogen=0;
      trivlabels_init(edgelist);
      aantal_gerichte_bogen=0;
    }
  return;

}

/*********************************SORT_DECREASING**********************************/

#define WISSEL(a,b) { buffer=(a); (a)=(b); (b)=buffer; }

void sort_decreasing(int list[], int number)
{
  // only very small numbers have to be sorted, so better simple

  int buffer,i,j;

  if (number==2)
    { if (list[0]<list[1]) WISSEL(list[0],list[1])
      return;			     
    }
  if (number==3)
    {
      if (list[0]<list[1])
	{ if (list[1]<list[2]) // l0<l1<l2
	    { WISSEL(list[0],list[2]) }
	  else // l1>l0 l1>l2
	    { if (list[0]<list[2]) // l0<l2<l1
		{ buffer=list[0]; list[0]=list[1]; list[1]=list[2]; list[2]=buffer; }
	      else // l2<l0<l1
		{ WISSEL(list[0],list[1]) }
	    }
	}
      else // l1<l0
	{ if (list[0]<list[2]) // l1<l0<l2
	    { buffer=list[0]; list[0]=list[2]; list[2]=list[1]; list[1]=buffer; }
	  else // l1<l0 l2<l0
	    { if (list[1]<list[2]) // l1<l2<l0
		{ WISSEL(list[2],list[1]) }
	      //else // l2<l1<l0 -- OK
	    }
	}
      return;
    }
  if (number>3) //should hardly ever happen and then "number" is still small -- bubblesort
    {
      for (i=number-1; i ; i--)
	for (j=0;j<i;j++)
	  if (list[j]<list[j+1]) WISSEL(list[j],list[j+1])
    }				   

  return;

}



/**********************************COMPUTE_IMAGE_OPERATION**************************/

int compute_image_operation(unsigned char image[], unsigned char original[], permutation generator[])
// returns 1 if image and original differ and 0 otherwise
{
  int in_image[MAXN], out_image[MAXN], double_image[MAXN];
  int run, num_in, num_out, num_double, change;

  change=0;

  image[0]=generator[original[0]];

  if (original[0]!=image[0]) change=1;

  for (run=1, num_in=0; original[run]!=INFTY_UCHAR; run++) 
    { in_image[num_in]=generator[original[run]]; num_in++; }

  for (run++, num_out=0; original[run]!=INFTY_UCHAR; run++) 
    { out_image[num_out]=generator[original[run]]; num_out++; }

  for (run++, num_double=0; (run<blocklength) && (original[run]!=INFTY_UCHAR); run++) 
    { double_image[num_double]=generator[original[run]]; num_double++; }

  if (num_in>1) sort_decreasing(in_image,num_in);
  if (num_out>1) sort_decreasing(out_image,num_out);
  if (num_double>1) sort_decreasing(double_image,num_double);


  for (run=1 ; num_in; run++) 
    { num_in--; image[run]=in_image[num_in]; if (original[run]!=image[run]) change=1;}
  image[run]=INFTY_UCHAR; run++;
  for ( ; num_out; run++) 
    { num_out--; image[run]=out_image[num_out]; if (original[run]!=image[run]) change=1;}
  image[run]=INFTY_UCHAR; run++;
  for ( ; num_double; run++) 
    { num_double--; image[run]=double_image[num_double]; if (original[run]!=image[run]) change=1;}
  image[run]=INFTY_UCHAR;
  //for ( ;run<blocklength;run++) image[run]=INFTY_UCHAR;

  return change;

}


/**********************************CONSTRUCT_OPERATIONS_ONE**************************/

void construct_operations_one(int center, int end, int do_inedge, int do_outedge, int do_double)
{
  static unsigned char *buffer;
  unsigned char *start;

  if (size_operations<=((number_operations+1)*blocklength)) 
    { buffer=operations;
      operations=malloc((size_t)(2*size_operations));
      if (operations==NULL) 
	{ fprintf(stderr,"Can't allocate %d bytes for operations -- exiting\n",2*size_operations); 
	  exit(1); }
      memcpy(operations,buffer,size_operations);
      free(buffer);
      size_operations *= 2;
    }

  if (do_inedge && (indeg[center] < maxindeg) && (outdeg[end]<maxoutdeg))
    { start=operations+(number_operations*blocklength);
      start[0]=center; start[1]=end; start[2]=start[3]=start[4]=INFTY_UCHAR;
      number_operations++;
    }

  if (do_outedge && (indeg[end] < maxindeg) && (outdeg[center]<maxoutdeg))
    { start=operations+(number_operations*blocklength);
      start[0]=center; start[2]=end; start[1]=start[3]=start[4]=INFTY_UCHAR;
      number_operations++;
    }

  if (do_double && double_free[center] && 
      (indeg[end] < maxindeg) && (indeg[center] < maxindeg) &&
      (outdeg[center]<maxoutdeg) && (outdeg[end]<maxoutdeg))
    { start=operations+(number_operations*blocklength);
      start[0]=center; start[3]=end; start[1]=start[2]=start[4]=INFTY_UCHAR;
      number_operations++;
    }

  return;

}



/**********************************CONSTRUCT_OPERATIONS_FINAL**************************/

void construct_operations_final(int list[], int decided[],unsigned char buffer_op[], 
				int positie, int numberin, int numberout)
{
  int i, center;
  static unsigned char *buffer;

  /* here about position "positie" of the operation is decided. If "positie==blocklength"
     we are done. 

     This function decides on double edges.

  */


  if (size_operations<=((number_operations+1)*blocklength)) 
    { buffer=operations;
      operations=malloc((size_t)(2*size_operations));
      if (operations==NULL) 
	{ fprintf(stderr,"Can't allocate %d bytes for operations -- exiting\n",2*size_operations); 
	  exit(1); }
      memcpy(operations,buffer,size_operations);
      free(buffer);
      size_operations *= 2;
    }

  center=buffer_op[0];
  // double_free[center] must be tested inside the loop to make sure
  // that it only leads to rejection if all edges have already been assigned
  
  if (double_free[center])// otherwise the rest was filled in in the previous step
	for (i=0;list[i]>=0;i++)
	  if (!decided[i])
	    { 
	      if ( double_free[center] && double_free[list[i]] &&
		   ((outdeg[center]+numberout < maxoutdeg) && (indeg[list[i]]<maxindeg)) &&
		   ((indeg[center]+numberin < maxindeg) && (outdeg[list[i]]<maxoutdeg)))
		{
		  buffer_op[positie]=list[i];
		  positie++; numberout++; numberin++;
		}
	      else return;
	    }

  buffer_op[positie]=INFTY_UCHAR; positie++;

  memcpy(operations+(number_operations*blocklength),buffer_op,positie);
  number_operations++;
  return;

}


/**********************************CONSTRUCT_OPERATIONS_OUT**************************/

void construct_operations_out(int list[],int liststart, int decided[],unsigned char buffer_op[], 
			      int positie, int numberin, int numberout, int lowerlimit_outdeg,
			      int mindouble)
{
  int i, center, buffer;
  static int outdeg_larger;

  /* here about position "positie" of the operation is decided. If "positie==blocklength"
     we are done. 

     This function decides on outgoing edges.

     liststart ensures that no earlier -- smaller -- entries are written after larger ones
     unless there has been an INFTY_UCHAR to separate them.
  */

  center=buffer_op[0];

  if (numberout==0) // just test once
    { buffer=outdeg[center] + tobedirected[center] -numberin;
      if (buffer<lowerlimit_outdeg) return;
      else
	{ if (buffer>lowerlimit_outdeg) outdeg_larger=1; else outdeg_larger=0; }
    }

  if (double_free[center]) // then there is still something to decide
    { if (outdeg[center]+numberout < maxoutdeg)
	for (i=liststart;list[i]>=0;i++)
	  if (!decided[i])
	    { 
	      if (indeg[list[i]]<maxindeg)
		{ decided[i]=1;
		  buffer_op[positie]=list[i];
		  construct_operations_out(list,i+1,decided,buffer_op,positie+1,
					     numberin, numberout+1,lowerlimit_outdeg,mindouble);
		  decided[i]=0;
		}
	    }
      buffer_op[positie]=INFTY_UCHAR;
      if (outdeg_larger || 
	  (indeg[center]+outdeg[center]+(tobedirected[center]<<1)-numberin-numberout-deg[center]>=mindouble))
	  // looks complicated but is just the final number of double edges
      construct_operations_final(list,decided,buffer_op,positie+1,
			   numberin, numberout);
      return;
    }
  //else // all the undecided entries have to be outgoing
  
    for (i=liststart;list[i]>=0;i++)
      if (!decided[i])
	{ 
	  if ((outdeg[center]+numberout < maxoutdeg) && (indeg[list[i]]<maxindeg))
	    { 
	      buffer_op[positie]=list[i];
	      positie++; numberout++;
	    }
	  else return;
	}
  buffer_op[positie]=INFTY_UCHAR;
  construct_operations_final(list,decided,buffer_op,positie+1,
			     numberin, numberout);
  return;
    
}





/**********************************CONSTRUCT_OPERATIONS_IN**************************/

void construct_operations_in(int list[],int liststart, int decided[],unsigned char buffer_op[], 
			     int positie, int numberin, int numberout, int lowerlimit_outdeg,
			     int mindouble)
{
  int i, center;

  /* here about position "positie" of the operation is decided. If "positie==blocklength"
     we are done. 

     This function decides on incoming edges.

     liststart ensures that no earlier -- smaller -- entries are written after larger ones
     unless there has been an INFTY_UCHAR to separate them.
  */


  center=buffer_op[0];

  if(indeg[center]+numberin < maxindeg)
    for (i=liststart;list[i]>=0;i++)
      if (!decided[i])
	{ 
	  if (outdeg[list[i]]<maxoutdeg)
	    { decided[i]=1;
	      buffer_op[positie]=list[i];
	      construct_operations_in(list,i+1,decided,buffer_op,positie+1,
				      numberin+1, numberout,lowerlimit_outdeg,mindouble);
	      decided[i]=0;
	    }
	}
  buffer_op[positie]=INFTY_UCHAR;
  construct_operations_out(list,0,decided,buffer_op,positie+1,
			   numberin, numberout,lowerlimit_outdeg,mindouble);
  return;
}


/**********************************CONSTRUCT_EXTENSIONS**************************/

void construct_extensions(int still_open[], int orbit[], graph touched, int first_in_orbit, graph sameorbit)
{
  int top, top2, j, end, list[MAXN], decided[MAXN], error, lowerlimit_outdeg, readylist[MAXN], *readyrun, dummy;
  int minout, do_double=0, i, mindouble;
  unsigned char buffer[MAXN+4];
  graph pre_free, tbd1, bufferg, freevertices; // tbd1 is the set of vertices with exactly one edge to be directed left
  // decided[i]=0 if and only if about list[i] is already decided //

  tbd1=(graph)0;
  pre_free= all & ~touched;
  readyrun=readylist;
  lowerlimit_outdeg=0;

  // Since either all or no vertices in the orbit have edges going to untouched vertices, one could split
  // this function to save some tests. But for the moment it should better stay like this


  // if a vertex v has at least two edges to be directed and one of them leads to a vertex in the same
  // orbit that will be ready afterwards, then directing the edges around v would not have a
  // reverse operation. So these operations must be avoided.

  // if only one edge is left to be directed and the endpoint is in the same orbit is must be directed double
  // or outgoing.

if (first_in_orbit)
  {
    for (number_operations=0; (top=(*still_open))>=0; still_open++)
      // number_operations is global, so it must not be passed to the next functions
      //if (NOT_READY(top))
      { 
	if (tobedirected[top]==1)
	  {
	    for (j= FIRSTBIT(workg[top]); READY(j) ; NEXTEL((workg[top]),(j))); 
	    if (ISELEMENT(&sameorbit,j) && (tobedirected[j]==1)) construct_operations_one(top, j, 0, 1, double_allowed);
	    else construct_operations_one(top, j, 1, 1, double_allowed);
	  }
	else
	  if (!(staticg[top] & tbd1))
	    { 
	      for (end=0, j= FIRSTBIT(workg[top]); (j<WORDSIZE); NEXTEL((workg[top]),(j))) 
		{ if (NOT_READY(j)) { list[end]=j; decided[end]=0; end++; }}
	      // the edge has already been directed if and only if the other top has already been finished
	      list[end]= -1;
	      buffer[0]= top;
	      construct_operations_in(list,0,decided,buffer,1,0,0,0,0);
	    }
      }
    return;
  }

// else: not first in orbit

 minout=INT_MAX;

 if (staticg[orbit[0]]& ~touched) // all operations will contain edges to free vertices
   {
     for (j=0; (top=orbit[j])>=0; j++)
       { if (NOT_READY(top)) 
	   { //ADDELEMENT(&sameorbit,top);
	     if (tobedirected[top]>=2) ADDELEMENT(&pre_free,top); // will stay free unless is chosen top
	     else ADDELEMENT(&tbd1,top);
	   }
	 else // that is: ready
	   { *readyrun=top; readyrun++; } // minout would not be used 
       }
   }
 else
   {
     do_double=double_allowed;
     for (j=0; (top=orbit[j])>=0; j++)
       { 
	 if (NOT_READY(top)) 
	   { //ADDELEMENT(&sameorbit,top);
	     if (tobedirected[top]>=2) ADDELEMENT(&pre_free,top); // will stay free unless is chosen top
	     else ADDELEMENT(&tbd1,top);
	   }
	 else // that is: ready
	   { bufferg=(workg[top]&sameorbit);
	     *readyrun=top; readyrun++; 
	     if ((bufferg) && (outdeg[top]<minout)) { minout=outdeg[top]; do_double=double_allowed; }
	     if (do_double && (outdeg[top]==minout))
	       { 
		 FORALLELEMENTS(bufferg,i)
		   { if (READY(i) && !ISELEMENT(workg+i,top)) // no double edge, so double edges inside the orbit would not be accepted
		       { do_double=0; i=WORDSIZE; }
		   }
	       }
	   }
       }
   }
 *readyrun= -1;

for (number_operations=0; (top=(*still_open))>=0; still_open++)
  // number_operations is global, so it must not be passed to the next functions
  //if (NOT_READY(top))
  { 
    if (tobedirected[top]==1) // already ready vertices can't have less...
      { 
	for (j= FIRSTBIT(workg[top]); READY(j) ; NEXTEL((workg[top]),(j))); 
	if (ISELEMENT(&sameorbit,j) && (tobedirected[j]==1)) 
	  { bufferg=workg[j]&sameorbit; DELELEMENT(&bufferg,top);
	    if (outdeg[top]<minout)
	      { if ((bufferg== (graph)0) || (outdeg[j]>outdeg[top])) // otherwise the end would be better 
		  {
		    if (outdeg[top]==minout-1) construct_operations_one(top, j, 0, 1, do_double);
		    else construct_operations_one(top, j, 0, 1, double_allowed);
		  }
		else 
		  if ((outdeg[j]==outdeg[top]) && double_allowed)
		  {
		    if (outdeg[top]==minout-1) construct_operations_one(top, j, 0, 0, do_double);
		    else construct_operations_one(top, j, 0, 0, double_allowed);
		  }
	      }
	  }
	else construct_operations_one(top, j, 1, 1, double_allowed);

      }
    else 
      { 
	if (!(staticg[top] & tbd1))
	  {
	    error=0;
	    freevertices = pre_free | (tbd1 & ~staticg[top]);
	    DELELEMENT1(&freevertices,top);
	    // now freevertices is the set of vertices that will still be free after  the addition of top
	    for (lowerlimit_outdeg=mindouble=0, readyrun=readylist; !error && (top2= *readyrun)>=0; readyrun++)
	      { bufferg= freevertices & staticg[top2];
		if ((dummy=POPCOUNT(bufferg)))// will be a candidate afterwards
		  { if (dummy < tobedirected[top]) error=1;
		    else
		      { if (dummy == tobedirected[top])
			  {
			    if (outdeg[top2]>lowerlimit_outdeg) 
			      { lowerlimit_outdeg=outdeg[top2]; mindouble=indeg[top2]+outdeg[top2]-deg[top2]; }
			    else 
			      if (double_allowed && (outdeg[top2]==lowerlimit_outdeg)) 
				{
				  if (indeg[top2]+outdeg[top2]-deg[top2]>mindouble) 
				    mindouble=indeg[top2]+outdeg[top2]-deg[top2];
				}
			  }
		      }
		  }
	      }
	    if (!error)
	      {
		for (end=0, j= FIRSTBIT(workg[top]); 
		     (j<WORDSIZE) ; NEXTEL((workg[top]),(j))) 
		  if (NOT_READY(j)) { list[end]=j; decided[end]=0; end++; }

		// the edge has already been directed if and only if the other top has already been finished
		list[end]= -1;
		buffer[0]= top;
		if (!error) construct_operations_in(list,0,decided,buffer,1,0,0,lowerlimit_outdeg,mindouble);
	      }
	  }
      }
  }
  return; 
}

int compare_op(unsigned char *op1, unsigned char *op2)
// returns 0 if operations are the same, something negative if op1
// is lexicographically smaller and something positive else
{
  int counter;

  if (*op1 != *op2) return (*op1 - *op2); // can never be INFTY_UCHAR
  op1++; op2++;
  counter=2;
  while (*op1 == *op2)
    { if (*op1 == INFTY_UCHAR) { if (!counter) return 0; counter--; }
      op1++; op2++;
    }
  return (*op1 - *op2);
}

/******************************SEARCH_OP*******************************/

int search_op(unsigned char *op)

/* searches the operation op in the list of operations and returns the index. */

{
  int min, max, center;

  min=0; max=number_operations-1; center= (min+max)/2;

  while (min<max)
    { 
      if (compare_op(op,operations+(blocklength*center))<=0)
	{ max=center; center= (min+max)/2; }
      else // op strictly larger 
	{ min=center+1; center= (min+max)/2; }
    }

  return min;
}

/******************************SEARCH_EDGE*******************************/

int search_edge(int start, int end, int edgelist[][2], int length)

/* searches the edge start->end in edgelist and returns the index. */

{
  int min, max, center;

  min=0; max=length-1; center= (min+max)/2;

  while (min<max)
    { 
      if ((start<edgelist[center][0]) || ((end<=edgelist[center][1]) && (start==edgelist[center][0]))) // start->end smaller or equal
	{ max=center; center= (min+max)/2; }
      else // start->end strictly larger 
	{ min=center+1; center= (min+max)/2; }
    }

  return min;
}


/******************************COMPUTE_ORBITS************************/

int compute_orbits()
/* compute the orbits of the group described by the global variable generators[][] on
   the operations described in the global variable operations.
   Returns the number of orbits.
*/
{
  int i,j,buffer,run, orbits, index;
  unsigned char image[MAXN+4];

  if (number_operations>size_root)
    { size_root=number_operations;
      free(root_op);
      root_op=malloc((size_t)size_root*sizeof(int));
    }
  if (root_op==NULL)  
    { fprintf(stderr,"Can't allocate %d items for root_op -- exiting.\n",size_root);
      exit(0); }

  for (i=0; i<number_operations; i++) root_op[i]=i;

  for (i=0, orbits=number_operations; i<number_operations; i++)
    for (j=0; j<number_of_generators; j++)
      { 
	if (compute_image_operation(image, operations+(i*blocklength), generators[j]))
	  {
	    index=search_op(image);
	    //unite:
	    while ((buffer=root_op[index])!=index) index=buffer; // no path compression here
	    run=i;
	    while ((buffer=root_op[run])!=run) { root_op[run]=index; run=buffer; }
	    if (run!=index) 
	      { orbits--; // the roots were different -- two trees are united
		root_op[run]=index; }
	  }
      }
  return orbits;
}

/******************************COMPUTE_EDGEORBITS****************************/

#define MAKEIMAGE(a,k) {imagex=generators[k][(a)[0]]; imagey=generators[k][(a)[1]];}

void compute_edgeorbits(int edgelist[][2],int orb[], int length)

{
  int i,j, index, buffer, run;
  int imagex, imagey;

  for (i=0; i<length; i++)
    for (j=0; j<number_of_generators; j++)
      { MAKEIMAGE(edgelist[i],j);
	index=search_edge(imagex, imagey, edgelist, length);
	while ((buffer=orb[index])!=index) index=buffer; // no path compression here
	run=i;
	while ((buffer=orb[run])!=run) { orb[run]=index; run=buffer; }
	orb[run]=index;
      }
}


/******************************CANONICAL****************************/

int canonical(unsigned char operation[], int vertexorbit[], int *newgroup, int orbitid, graph touched)

// checks whether the operation is canonical in vertexorbit. Returns 1 if yes, 0 otherwise.
// In newgroup it is stored whether a new group has been computed for this (*newgroup=1) or not
// (*newgroup=0)

// If no edges go to vertices outside the set of vertices in already "finished" orbits or to vertices in
// this orbit with still undecided edges, the canonical operation is directing one edge {x,y} x->y or x<->y.
// Edges completely in this orbit and obtained by operations center x end y are chosen as follows: 
// First the starting point x is chosen with minimum outdegree.
// For edges with starting points with the same minimal outdegree, edges are preferred which are single. 
// If the tested operation is double and single operations exist, the operation is rejected.
// Among the remaining edges those are chosen where y has minimum indegree.
// Among the remaining those are chosen where x is in the orbit of the one with smallest canonical number
// (after the operation). y is the neighbour with smallest indegree -- and among those with the same indegree
// the one with smallest canonical number (in the same given numbering of course).

// Otherwise:
// Define a vertex to be free if it is not in an orbit that was already chosen (we call that untouched)
// or it is in such an orbit but not all edges to it are already directed (that can only be in the last orbit).
// In this case the central vertex of the canonical operation is one that is chosen in steps:
// Among all with edges to free vertices it is chosen as one with a minimum number of these.
// Among those one with max outdegree is chosen.
// Among those one with maximum number of double edges is chosen.
// Finally a more complicated artificial criterion is tested -- not suitable for a look ahead during the construction
// of possible operations
// Among the remaining, one with smallest canonical label is chosen.

{ int i, j, k, min, top, buffer, candidatelist[MAXN], number_candidates, center;
  graph free, candidates, bufferset, dummy, sameorbit, bit_startlist;
  int edgelistcounter, edgelist[2*MAXN*MAXN][2], orb[2*MAXN*MAXN], canoncenter, canonend, endvertex=0;
  int finished[MAXN], finishedptn[MAXN], testoutdeg;
  int startlist[MAXN], *startrun; // list of possible starts in case of no candidates
  int numberends, ends[MAXN][MAXN], endin=0, *run;
  nvector *colour;
  long long int buffer2, k2;



  free= all & (~touched); // in fact all ready vertices are candidates if one has a non-empty
  // intersection with this as they are all in the same orbit -- so have the same number of edges
  // to this set
  sameorbit=bit_orbit[orbitid];

  for (run=vertexorbit; (*run)>=0; run++) if (NOT_READY(*run)) ADDELEMENT(&free,(*run));

  center=operation[0];

  candidates= (graph)0;
  number_candidates=0;
  min=INFTY_UCHAR;
  for (i=0; (top=vertexorbit[i])>=0; i++) 
    if (READY(top) && (bufferset=(staticg[top] & free)))
      {
	buffer=POPCOUNT(bufferset);
	if (buffer<min) { candidates= (graph)0; ADDELEMENT(&candidates,top); 
	                  candidatelist[0]=top; number_candidates=1; min=buffer; }
	else
	  if (buffer==min)
	    { ADDELEMENT(&candidates,top); candidatelist[number_candidates]=top; number_candidates++; }
      }
  if (number_candidates && !ISELEMENT(&candidates,center))
    { 
      return 0;
    }

  if (number_candidates==1) { *newgroup=0; return 1; }
  
  testoutdeg=outdeg[center];
  // First canonicity criterion: maximum outdegree
  if (number_candidates)
    {
      for (i=0; i<number_candidates; ) 
	{ k= outdeg[candidatelist[i]]-testoutdeg;
	  if (k>0) 
	    { return 0;}
	  else 
	    { if (k<0) { DELELEMENT(&candidates,candidatelist[i]); 
		         number_candidates--; 
			 candidatelist[i]=candidatelist[number_candidates];   }
	      else i++; // otherwise the new element has to be tested
	    }
	}
      
      if (number_candidates==1) { *newgroup=0; return 1; }
    // end number_candidates>0

      if (remaining_doubles != max_doubles) // there are double edges
	{ buffer=indeg[center]+outdeg[center]-deg[center]; // number of double edges starting at center
	  for (i=0; i<number_candidates; )
	    { j=candidatelist[i];
	      k= indeg[j]+outdeg[j]-deg[j]-buffer; 
	      if (k>0) 
		{ return 0;}
	      else 
		{ if (k<0) { DELELEMENT(&candidates,candidatelist[i]); 
		    number_candidates--; 
		    candidatelist[i]=candidatelist[number_candidates];   }
		  else i++; // otherwise the new element has to be tested
		}
	    }
	}
      if (number_candidates==1) { *newgroup=0; return 1; }
 
      //OK-- last try. A colour that can not really be used to avoid unnecessary operation in advance
      colour=rememberorbits[orbitid];
      dummy=workg[center]&(free | sameorbit);
      buffer2=1LL;
      FORALLELEMENTS(dummy,i) 
	{ buffer2 *= (tobedirected[i]<<12)+(indeg[i]<<9)+(outdeg[i]<<6)+(colour[i]<<3)+deg[i]+1;
	  buffer2 = buffer2%63241LL;
	}

      for (i=0; i<number_candidates; )
	{ 
	  dummy=workg[candidatelist[i]]&(free | sameorbit);
	  k2=1LL;
	  FORALLELEMENTS(dummy,j) 
	    { k2 *= (tobedirected[j]<<12)+(indeg[j]<<9)+(outdeg[j]<<6)+(colour[j]<<3)+deg[j]+1;
	      k2 = k2%63241LL;
	    }
	  k2 -= buffer2;

	  if (k2>0LL) 
	    { return 0;}
	  else 
	    { if (k2<0LL) { DELELEMENT(&candidates,candidatelist[i]); 
		number_candidates--; 
		candidatelist[i]=candidatelist[number_candidates];   }
	      else i++; // otherwise the new element has to be tested
	    }
	}
      
      if (number_candidates==1) { *newgroup=0; return 1; }
 


   } // end number candidates > 0


  // OK -- we have to work. First the canonical form must be computed:
  // Problem: the difference between double edges and edges that have not yet been directed cannot be 
  // detected in the datastructure graph. To this end we put "ready" vertices in an extra partition --
  // this way double edges cannot be mapped on undirected edges. This is only necessary if there are double edges.

  else // (number_candidates==0) // finished vertices in this orbit always have an internal edge!
    { startrun=startlist;
      bit_startlist=(graph)0;
      if (operation[2]==INFTY_UCHAR) // a double edge was added and all single edges starting at same outdeg are better
	{ endvertex= operation[3];
	  endin=indeg[endvertex];
	  for (i=0; (k=vertexorbit[i])>=0; i++)
	    if (READY(k) && (workg[k] & sameorbit))// has _outgoing_ edge into same orbit
	      // the end points are automatically ready and all in the same orbit -- otherwise there would have been a candidate
	      { 
		if (outdeg[k]<testoutdeg) return 0;
		if (outdeg[k]==testoutdeg) 
		  { numberends=0; 
		    bufferset= workg[k]&sameorbit;
		    FORALLELEMENTS(bufferset,j) 
		      { if (!ISELEMENT(workg+j,k)) return 0; // else: also double edge
			if (indeg[j]<endin) return 0;
			if (indeg[j]==endin) { ends[k][numberends]=j; numberends++; }
		      }
		    if (numberends)
		      { ends[k][numberends]= -1; 
			*startrun=k; startrun++; ADDELEMENT(&bit_startlist,k); 
		      }
		  }
	      }
	}
      else // het is dus een single edge
	{ endvertex= operation[2];
	  endin=indeg[endvertex];
	  for (i=0; (k=vertexorbit[i])>=0; i++)
	    if (READY(k) && (workg[k] & sameorbit))// has _outgoing_ edge into same orbit
	      // the end points are automatically ready and all in the same orbit -- otherwise there would have been a candidate
	      { 
		if (outdeg[k]<testoutdeg) return 0;
		if (outdeg[k]==testoutdeg) 
		  { numberends=0;
		    bufferset= workg[k]&sameorbit;
		    FORALLELEMENTS(bufferset,j) 
		      { if (!ISELEMENT(workg+j,k)) // also single
			  {
			    if (indeg[j]<endin) return 0;
			    if (indeg[j]==endin) { ends[k][numberends]=j; numberends++; }
			  }
		      }
		    if (numberends)
		      { ends[k][numberends]= -1;
			*startrun=k; startrun++; ADDELEMENT(&bit_startlist,k); 
		      }
		  }
	      }
	}
      *startrun= -1;
      dummy=workg[center] & sameorbit;
      if ((startlist[1]== -1) && (ends[center][1]== -1)) { *newgroup=0; return 1; } // maar 1 boog mogelijk
    } 


  if (double_allowed) // distinguish single and double edges for nauty()
    {
      for (i=j=k=0; i<aantal_toppen;i++)
	{ top=lab[orbitid][i];
	  if (READY(top)) { finished[k]=top; finishedptn[k]=1; k++; } 
	  else { bufferlab[j]=top; bufferptn[j]=1; j++; }
	  if (ptn[orbitid][i]==0)
	    { if (j) bufferptn[j-1]=0;
	      if (k) finishedptn[k-1]=0;
	    }
	}
      for (i=0;i<k;i++,j++) { bufferptn[j]=finishedptn[i]; bufferlab[j]=finished[i]; } 
    }
  else
    {
      memcpy(bufferlab,lab[orbitid],aantal_toppen*sizeof(nvector));
      memcpy(bufferptn,ptn[orbitid],aantal_toppen*sizeof(nvector));
    }
  number_of_generators=0;
  nauty(workg,bufferlab,bufferptn,NILSET,orbits,&options_directed_canon,&stats,workspace,100*MAXN,1,aantal_toppen,canong);
  *newgroup=1;

  //nautyuse[number_candidates]++;

  // Being a free vertex is invariant under automorphisms and the set of edges to free vertices depends only
  // on the central vertex. So the central vertex characterizes the whole operation and orbits of central
  // vertices uniquely correspond to orbits of operations.

  if (number_candidates>1)
    {      
      for (i=0; !ISELEMENT(&candidates,bufferlab[i]); i++); // zoekt candidaat met kleinste kanonische label
            
      if (orbits[center]==orbits[bufferlab[i]]) { return 1; }
      else { return 0; }
    }

// only remaining case: number_candidates=0 -- that is: only one edge between 2 vertices of vertexorbit was added

  // operation must be outgoing or double edge -- but this is guaranteed by the construction of operations

  //  for (i=0; !(READY(bufferlab[i]) && ISELEMENT(&sameorbit,bufferlab[i]) && (workg[bufferlab[i]] & sameorbit))
  for (i=0; !ISELEMENT(&bit_startlist,bufferlab[i]); i++);
  // Now we have the smallest vertex in vertexorbit with an edge to another vertex in vertexorbit

  canoncenter=bufferlab[i];

  if (orbits[center]!=orbits[canoncenter]) { return 0; }

  // Now look for the endvertex:

  min=INFTY_UCHAR;
  dummy= sameorbit & workg[canoncenter];

  // at this point it is known that none of the possible startvertices has an endvertex that
  // has smaller indegree than the endvertex of the operation tested (that is: endin)

  // now look for the smallest labelled neighbour:

  if (operation[2]==INFTY_UCHAR) // a double edge was added and any suitable neighbour is double
    for (i=0;!(ISELEMENT(&dummy,bufferlab[i]) && (indeg[bufferlab[i]]==endin)) && (i<aantal_toppen);i++);
  else // a single edge was added and only neighbours via single edges may be chosen
    for (i=0;ISELEMENT(workg+bufferlab[i],canoncenter) ||
	   !(ISELEMENT(&dummy,bufferlab[i]) && (indeg[bufferlab[i]]==endin)); i++);
  
  canonend=bufferlab[i];

  if (orbits[endvertex]!=orbits[canonend]) { return 0; }
    if ((canoncenter==center) && (canonend==endvertex)) { return 1; }

  // OK -- bad luck. We have to check whether the edges are really in the same orbit
  // First make a list of all edges that are candidates:

  edgelistcounter=0;
  for (startrun=startlist; (k= *startrun)>=0;startrun++)
    if (orbits[k]==orbits[canoncenter])
      { 
	for (i=0; (j=ends[k][i])>=0; i++)
	  { 
	    if (orbits[j]==orbits[canonend]) 
	      { edgelist[edgelistcounter][0]=k; edgelist[edgelistcounter][1]=j; 
		orb[edgelistcounter]=edgelistcounter; edgelistcounter++; }
	  }
      }

  compute_edgeorbits(edgelist,orb,edgelistcounter);

  i=search_edge(center,endvertex,edgelist,edgelistcounter);
  while(i!=orb[i]) i=orb[i];
  j=search_edge(canoncenter,canonend,edgelist,edgelistcounter);
  while(j!=orb[j]) j=orb[j];

  if (i==j) { return 1; } // same orbit as canonical edge
  else { return 0; }


}

int all_diff_colours(graph testset, int orbitid)
// returns 1 if all elements have some different vertex invariant and 0 otherwise
{
  nvector *colour;
  int i,buffer;

  RESETMARKS_COLOUR;
  colour=rememberorbits[orbitid];

  FORALLELEMENTS(testset,i) 
    { 
      buffer = (tobedirected[i]<<12)+(indeg[i]<<9)+(outdeg[i]<<6)+(colour[i]<<3)+deg[i];
      buffer &= 32767;
      if (ISMARKED_COLOUR(buffer)) return 0;
      MARK_COLOUR(buffer);
    }

  return 1;
}


void fill_edgelist_final(graph g[])
     /* writes the edges in gg into the list in a lexicographic way and initializes is_gericht 
        and positie*/

{
  int i,j,end;

  for (i=0, j=aantal_gerichte_bogen; i<aantal_toppen; i++) /* j is de positie in edgelist */
    while (g[i])
      { end=FIRSTBIT(g[i]); DELELEMENT(g+i,end);  DELELEMENT(g+end,i); 
	edgelist_final[j][0]=i; edgelist_final[j][1]=end;
	j++; 
      }
  return;
}



void fill_edgelist_order_final()
     /* works like fill_edgelist_order -- only that in this case some edges are already assigned a direction.
      They are not added to the list. The first entry is at position aantal_gerichte_bogen, so that the list
      can be used like the normal edgelist for directing in the trivial case.
      The global variable positie is not assigned any values.
*/

{
  int last_positie,i, buffer, olddeg, newdeg, beste, start, end;
  int list[MAXN][MAXN], listlen[MAXN]={0}; /* list[i] contains the vertices with degree i at that moment */
  int toppositie[MAXN], bufferdeg[MAXN]; /* toppositie[i] is de positie van top i in de lijst */
  int buren[MAXN]; /* de buren van de top waaraan gewerkt wordt */
  graph g[MAXN], notready;

  for (i=aantal_gerichte_bogen=0; i<aantal_toppen; i++) 
    { indeg_free[i]=maxindeg-indeg[i]; outdeg_free[i]=maxoutdeg-outdeg[i];
      aantal_gerichte_bogen += (deg[i]-tobedirected[i]); }
  aantal_gerichte_bogen = aantal_gerichte_bogen>>1; // every edge was counted twice


  memcpy(g,workg,aantal_toppen*sizeof(graph));
  EMPTYSET1(&notready,1);

  // make a list of vertices according to the whole degree. Vertices with large degree give stronger restrictions
    for (i=0; i<aantal_toppen; i++) 
      if (NOT_READY(i)) { buffer=bufferdeg[i]=deg[i]; list[buffer][listlen[buffer]]=i; 
	                  toppositie[i]=listlen[buffer]; (listlen[buffer])++; 
			  ADDELEMENT(&notready,i); }
    else EMPTYSET1(g+i,1);

    for (i=0; i<aantal_toppen; i++) g[i] &= notready; // now only edges that are still to be directed are in g[]

    if (nodegbound || (aantal_bogen-aantal_gerichte_bogen<7)) { fill_edgelist_final(g); return; }

    last_positie=aantal_bogen-1;

    while (last_positie>=aantal_gerichte_bogen)
      {
	for (beste=1;listlen[beste]==0; beste++); // zoek kleinste aanwezige nog toe te voegen graad
	(listlen[beste])--; /* zo wordt hij ook verwijderd */
	start=list[beste][listlen[beste]];
	for (i=0; g[start] != (graph)0; i++) /* alle buren opslaan*/
	  { end=buren[i]=FIRSTBIT(g[start]); DELELEMENT(g+start,end); }
	buren[i]= -1;
	for (i=0; (end=buren[i])>=0; i++) /* alle buren */
	  { 
	    edgelist_final[last_positie][0]=start; edgelist_final[last_positie][1]=end; 
	    last_positie--;
	    DELELEMENT(g+end,start);
	    /* de buur verhuizen: */
	    olddeg=bufferdeg[end];
	    (bufferdeg[end])--;
	    newdeg=bufferdeg[end];
	    /* uit de oude lijst verwijderen */
	    if (listlen[olddeg]==1) listlen[olddeg]=0;
	    else { (listlen[olddeg])--;
	    buffer= list[olddeg][listlen[olddeg]];
	    list[olddeg][toppositie[end]]=buffer;
	    toppositie[buffer]=toppositie[end]; }
	    /* tot de nieuwe toevoegen */
	    if (newdeg)
	      { list[newdeg][listlen[newdeg]]=end;
	        toppositie[end]=listlen[newdeg];
	        (listlen[newdeg])++;
	      }
	  }
      }
    return;
}

/******************************CHOOSEORBIT***************************/

void chooseorbit(graph *touched, int best_orbit[], int orbitid)
{ int i,j,k, num_in_orbit[MAXN], dummy, inorbit[MAXN][MAXN], best, bestroot=0;

  // the group must be up to date at this point!

     for (i=0;i<aantal_toppen;i++) { num_in_orbit[i]=0; }
      for (i=0;i<aantal_toppen;i++)
	{ dummy=orbits[i]; 
	  inorbit[dummy][num_in_orbit[dummy]]=i;
	  num_in_orbit[dummy]++;
	  if (READY(i)) ADDELEMENT(touched,i);
	}

      
      best=INT_MAX;
      for (i=k=0;i<aantal_toppen;i++) 
	{ 
	  for (j=0;j<num_in_orbit[i]-1; j++) { ptn[orbitid][k]=1; lab[orbitid][k]=inorbit[i][j]; k++; } 
	  if (num_in_orbit[i]) { ptn[orbitid][k]=0; lab[orbitid][k]=inorbit[i][j]; k++; }
	  if ((orbits[i]==i) && NOT_READY(i) && (num_in_orbit[i] + (tobedirected[i]<<1))<best)
	    { best=num_in_orbit[i] + (tobedirected[i]<<1); bestroot=i; }
	
	}

      memcpy(rememberorbits[orbitid],orbits,sizeof(nvector)*aantal_toppen);
      
      // here the "best" orbit is chosen by a combination of number of vertices and edges to be directed per vertex
	
      orbitchoices[num_in_orbit[bestroot]]++;

      bit_orbit[orbitid]= (graph)0;
      for (i=j=0;i<aantal_toppen;i++) 
	if (orbits[i]==bestroot) { best_orbit[j]=i; j++; ADDELEMENT(bit_orbit+orbitid,i); }
      (*touched) |= bit_orbit[orbitid];
      best_orbit[j]= -1;  // the sign that the end of the list has been reached

      return;
}

/******************************CHOOSE_TRIV_ORBIT***************************/

// This routine does not really choose an orbit based on group computations, but
// just fixes a one vertex orbit detected by having a unique colour.
// It also updates lab[], ptn[] etc. Orbitid must be at least 1

void choose_triv_orbit(graph *touched, int best_orbit[], int orbitid, int fixed_vertex)
{ int i,oud;

  oud=orbitid-1;

  // take the old partitition and only move fixed_vertex to its own partition
  for (i=0;lab[oud][i]!=fixed_vertex; i++) { ptn[orbitid][i]=ptn[oud][i]; lab[orbitid][i]=lab[oud][i]; }
  if ((ptn[oud][i]==0) && (i>0)) ptn[orbitid][i-1]=0;
  for ( i++; i<aantal_toppen; i++) { ptn[orbitid][i-1]=ptn[oud][i]; lab[orbitid][i-1]=lab[oud][i]; }
  lab[orbitid][aantal_toppen-1]=fixed_vertex; ptn[orbitid][aantal_toppen-1]=0;

  memcpy(rememberorbits[orbitid],rememberorbits[orbitid-1],sizeof(nvector)*aantal_toppen);

  bit_orbit[orbitid]= (graph)0;
  ADDELEMENT(bit_orbit+orbitid,fixed_vertex);
  ADDELEMENT(touched,fixed_vertex);
  best_orbit[0]=fixed_vertex; best_orbit[1]= -1;

  return;
}

void prepare_next_step(int group_uptodate, int orbitid, int iterationdepth, graph touched)
{
  int i, j, k, d, end;
  int local_todolist[MAXN+1];
  graph not_ready;
  int colour[MAXN],numbercolour[MAXN], element[MAXN], diffcolours;
  nvector *orbcolour; 


// all vertices were in orbits
      not_ready= all & ~touched; end=1;
      FORALLELEMENTS(not_ready,j)
	{ if (NOT_READY(j)) { end=0; } else DELELEMENT(&not_ready,j); }
      if (end) { MAYBEPROCESS; WRITEUP(); return; } // no edges left to direct
      //else 

       
      if (!group_uptodate)
	{ end=1;
	  orbcolour=rememberorbits[orbitid];
	  diffcolours=0; 
	  FORALLELEMENTS(not_ready,i)
	    { numbercolour[diffcolours]=1; element[diffcolours]=i;
	      colour[diffcolours]= (tobedirected[i]<<12)+(indeg[i]<<9)+(outdeg[i]<<6)+(orbcolour[i]<<3)+deg[i];
	      for (k=0, d=1; k<diffcolours; k++) if (colour[k]==colour[diffcolours]) { end=d=0; numbercolour[k]++; }
	      if (d) diffcolours++; // different from all earlier colours
	    }
	
	  if (end) // group acts definitely trivial on rest
	    { 
	      fill_edgelist_order_final();
	      laatstepositie=edgelist_final+aantal_bogen-1;
	      trivlabels_init(edgelist_final+aantal_gerichte_bogen);
	      return;
	    }
	  
	  d= -1;
	  for (i=0; i<diffcolours; i++)
	    if (numbercolour[i]==1)
	      { if (tobedirected[element[i]]==1)
		  { 
		    choose_triv_orbit(&touched, local_todolist, orbitid+1, element[i]);
		    directorbit(local_todolist,local_todolist,0,orbitid+1,iterationdepth,touched);
		    return;
		  }
		else
		  { if (all_diff_colours(workg[element[i]]&not_ready,orbitid)) d=element[i]; }
	      }
	  // best is just one to direct, but d is reserve

	  if (d>=0)
	    { 
	      choose_triv_orbit(&touched, local_todolist, orbitid+1, d);
	      directorbit(local_todolist,local_todolist,0,orbitid+1,iterationdepth,touched);
	      return;
	    }
	  
	  // nieuwe orbit berekenen:
	  //if (!group_uptodate)
	  //	{ 
	  number_of_generators=0;
	  // since the orbit is complete, we don't NEED to distinguish between ready vertices or not
	  // to mark double edges. Tests showed that it also doesn't help.
	  memcpy(bufferlab,lab[orbitid],aantal_toppen*sizeof(nvector));
	  memcpy(bufferptn,ptn[orbitid],aantal_toppen*sizeof(nvector));
	  number_of_generators=0;
	  nauty(workg,bufferlab,bufferptn,NILSET,orbits,&options_directed,&stats,workspace,100*MAXN,1,aantal_toppen,NULL);
	  //orbitchoose_nauty++;
	}    

      end=1; // can we end the groupcomputations and switch to trivgroup?
      if (stats.numorbits<aantal_toppen) FORALLELEMENTS(not_ready,j) if (orbits[j]!=j) end=0;

      if (end)
	{ 
	  fill_edgelist_order_final();
	  laatstepositie=edgelist_final+aantal_bogen-1;
	  trivlabels_init(edgelist_final+aantal_gerichte_bogen);
	  return;
	}
     
      // --ELSE--
      /* Now the colouring for level orbitid can be prepared */
      
      chooseorbit(&touched, local_todolist, orbitid+1);

      directorbit(local_todolist,local_todolist,1,orbitid+1,iterationdepth,touched); 
      // no vertexenvironment directed, so does not need to be increased. Otherwise MAXN would not be an upper bound
      return;
}


/******************************DIRECTORBIT***************************/

void directorbit(int todo_list[], int vertexorbit[], int group_uptodate, int orbitid, 
		 int iterationdepth, graph touched)

/* number_done is number of vertices in this orbit that is already directed. The automorphism group
   does of course not interchange completed vertices from the orbit and uncompleted ones... 

   orbitid is the number of the orbit to be directed. In the beginning you start directing the first orbit (id=1),
   when that is done the second orbit (id=2), etc...
*/

{

  int i, j, k, top, localblocklength, number_orbits, remembertobedirected, center, end, group_OK;
  unsigned char *run, *runoperation;
  int local_todolist[MAXN+1], local_left_to_do, first_in_orbit, doit;
  int finished[MAXN], finishedptn[MAXN];
  graph dummy;
  // colour stores the different colours, numbercolour[] the number of graphs with that colour and
  // element an example element with that colour


#define START(x) ((x)*localblocklength)

  for (i=local_left_to_do=0;todo_list[i]>=0; i++) 
     if (NOT_READY(todo_list[i])) 
       { local_todolist[local_left_to_do]=todo_list[i]; local_left_to_do++; }
  local_todolist[local_left_to_do]= -1;

  if (local_left_to_do==0) 
    { if (touched==all) { MAYBEPROCESS; WRITEUP(); return; }
      else prepare_next_step(group_uptodate, orbitid, iterationdepth, touched);
      return;
    }

  // else -- there are still vertices in the orbit
  if (vertexorbit[local_left_to_do]== -1) first_in_orbit=1; else first_in_orbit=0;
  if (first_in_orbit) // new orbit to start with
    orbitblocklength[orbitid]=tobedirected[vertexorbit[0]]+4;
  localblocklength=blocklength=orbitblocklength[orbitid];
  /* als een top al afgewerkt is kan de bloklengte voor de  verschillende toppen in wat vroeger
     een orbit was verschillen. orbitblocklength[orbitid] is altijd een bovengrens */

  construct_extensions(local_todolist,vertexorbit,touched, first_in_orbit, bit_orbit[orbitid]);

  if (number_operations==0) return;

  // in case of one element with tobedirected>1 the group is up to date

  if ((group_uptodate && (stats.numorbits==aantal_toppen))
      || ((vertexorbit[1]== -1) && (tobedirected[vertexorbit[0]]==1))) doit=0;
  else doit=1;
  // doit means: really check the group

  if (doit)
    {
      dummy= all & (~touched);
      for (i=0; (k=local_todolist[i])>=0; i++) ADDELEMENT(&dummy,k);
      if (all_diff_colours(dummy,orbitid)) doit=0;
      else
	{
	  if (!group_uptodate)
	    { /* we need only the group -- no canonical numbering */
	      number_of_generators=0;
	      if (double_allowed)
		{
		  for (i=j=k=0; i<aantal_toppen;i++)
		    { top=lab[orbitid][i];
		      if (READY(top)) { finished[k]=top; finishedptn[k]=1; k++; } 
		      else { bufferlab[j]=top; bufferptn[j]=1; j++; }
		      if (ptn[orbitid][i]==0)
			{ if (j) bufferptn[j-1]=0;
			  if (k) finishedptn[k-1]=0;
			}
		    }
		  for (i=0;i<k;i++,j++) { bufferptn[j]=finishedptn[i]; bufferlab[j]=finished[i]; } 
		}
	      else
		{
		  memcpy(bufferlab,lab[orbitid],aantal_toppen*sizeof(nvector));
		  memcpy(bufferptn,ptn[orbitid],aantal_toppen*sizeof(nvector));
		}
	      nauty(workg,bufferlab,bufferptn,NILSET,orbits,&options_directed,&stats,workspace,100*MAXN,1,aantal_toppen,NULL);
	      //orbitcomp_nauty++;
	    }
	  
	  doit=0;
	  for (i=0; (k=local_todolist[i])>=0; i++) if (orbits[k]!=k) doit=1;
	  if (!doit)
	    { dummy= all & (~touched);
	      FORALLELEMENTS(dummy,j) if (orbits[j]!=j) doit=1;
	    }
	}
    }

  if (doit) number_orbits=compute_orbits(); else number_orbits=number_operations;

  if (remember_size[iterationdepth]<(number_orbits*localblocklength))
    { free(remember_operations[iterationdepth]);
      remember_size[iterationdepth] = 2*number_orbits*localblocklength;
      remember_operations[iterationdepth]=malloc((size_t)remember_size[iterationdepth]);
      if (remember_operations[iterationdepth]==NULL)
	{ fprintf(stderr,"Can't allocate %d items to store orbits -- exiting.\n",remember_size[iterationdepth]); 
	  exit(3); }
    }

  if (doit)
    {
      for (i=j=0; i<number_operations; i++)
	{ if (root_op[i]==i)
	    { memcpy(remember_operations[iterationdepth]+START(j),operations+START(i),localblocklength);
	      j++; 
	    }
	}
    }
  else // just copy
    { number_orbits=number_operations;
      memcpy(remember_operations[iterationdepth],operations,localblocklength*number_operations);
    }

  // Now each possible nonequivalent extension will be applied exactly once:

  for(i=0;i<number_orbits; i++) // in the loop edges are directed and datatypes updated -- then the next iteration
    // is called and afterwards things are reset
    { 
      runoperation=remember_operations[iterationdepth]+START(i);
      center=runoperation[0];
      remembertobedirected=tobedirected[center];
      for (run=runoperation+1; *run!=INFTY_UCHAR; run++) // first incoming edges
	{ end= *run;
	  DELELEMENT(workg+center,end);
	  indeg[center]++; outdeg[end]++; tobedirected[end]--;
	}
      for (run++; *run!=INFTY_UCHAR; run++) // then outgoing edges
	{ end= *run;
	  DELELEMENT(workg+end,center);
	  indeg[end]++; outdeg[center]++; tobedirected[end]--;
	}
      for (run++; *run!=INFTY_UCHAR; run++) // then double edges
	{ end= *run; remaining_doubles--;
	  indeg[end]++; outdeg[center]++; tobedirected[end]--;
	  indeg[center]++; outdeg[end]++;
	}
      tobedirected[center]=0;
      // ready for next round

      group_OK=0;

      if ((first_in_orbit) || canonical(runoperation, vertexorbit, &group_OK, orbitid, touched))
	directorbit(local_todolist, vertexorbit, group_OK, orbitid, iterationdepth+1, touched);

      // now reset everything

     for (run=runoperation+1; *run!=INFTY_UCHAR; run++) // first incoming edges
	{ end= *run;
	  ADDELEMENT(workg+center,end);
	  indeg[center]--; outdeg[end]--; tobedirected[end]++;
	}
      for (run++; *run!=INFTY_UCHAR; run++) // then outgoing edges
	{ end= *run;
	  ADDELEMENT(workg+end,center);
	  indeg[end]--; outdeg[center]--; tobedirected[end]++;
	}
      for (run++; *run!=INFTY_UCHAR; run++) // then double edges
	{ end= *run; remaining_doubles++;
	  indeg[end]--; outdeg[center]--; tobedirected[end]++;
	  indeg[center]--; outdeg[end]--;
	}
      tobedirected[center]=remembertobedirected;
    }

  return;

}

/******************************DIRECT_ALL_NONTRIV***************************/

void direct_all_nontriv()
/* this functions starts to direct the edges in case of nontrivial automorphism.
   nauty() has just been called, so orbit etc are all up to date */

{
  int i,best_orbit[MAXN+1];
  int num_in_orbit[MAXN];
  graph touched;

  addnumber=1;
  memcpy(rememberorbits[0],orbits,sizeof(nvector)*aantal_toppen);

  all= (graph)0;
  for (i=0;i<aantal_toppen;i++) { num_in_orbit[i]=0; ADDELEMENT(&all,i); }

  touched= (graph)0;
  chooseorbit(&touched, best_orbit, 0);

  directorbit(best_orbit,best_orbit,1,0,0,touched);

  return;
}



void init_waterplugin(int toppen)
{ 
  
 watermaxdeg=maxindeg+maxoutdeg;
 //if (watermaxdeg<maxdeg) maxdeg=watermaxdeg; That was for the routine that was a plugin for geng
 /* grafen met een grotere graad kunnen niet gericht worden zodat ze nog aan de 
    voorwaarden voldoen */
 
 if (maxindeg<=maxoutdeg) { mingerichtdeg=maxindeg; watermaxedges=toppen*maxindeg; }
 else { mingerichtdeg=maxoutdeg; watermaxedges=toppen*maxoutdeg; }
 //if (watermaxedges<maxe) maxe=watermaxedges; That was for the routine that was a plugin for geng
 /* grafen met meer bogen kunnen niet gericht worden, zodat ze nog aan de 
    voorwaarden voldoen */

}


void writeop(unsigned int op, BOOG orbit[], int length)
{
  int i, start, end;
  fprintf(stderr,"\n");
  for (i=0;i<length;i++)
    { start=orbit[i][0]; end=orbit[i][1];
      if (GETTYPE(op,i)==0) fprintf(stderr,"(%d->%d) ",start, end);
      else if (GETTYPE(op,i)==1) fprintf(stderr,"(%d<->%d) ",start, end);
      else fprintf(stderr,"(%d<-%d) ",start, end);
    }
  fprintf(stderr,"\n");
  return;
}


 void writegraph_edgeorb(graph *g, int aantal_toppen,int aantal_bogen, int aantal_gerichte_bogen)
 // for graphs as they are used in the edgeorbit routines
   { int i,j;
   fprintf(stderr,"---------------------------------------------------------\n");
   fprintf(stderr,"Graph with %d vertices, %d edges, %d already directed\n",aantal_toppen,aantal_bogen, aantal_gerichte_bogen);
   for (i=0; i<aantal_toppen;i++) 
     { fprintf(stderr,"%d:",i);
       FORALLELEMENTS(g[i],j) 
	 { fprintf(stderr," %d", j);
	 if (aantal_gerichte_bogen)
	   {
	 if (!is_gericht[i][j]) fprintf(stderr," (ng)");
	 else { if (ISELEMENT(g+j,i)) fprintf(stderr," (D)");
	        else fprintf(stderr," (S)");
	      }
	   }
	 }
       fprintf(stderr,"\n");
     }
   fprintf(stderr,"---------------------------------------------------------\n");
   }




void fill_edgelist_edgeorb()
     /* writes the edges in gg into the list in a lexicographic way and initializes is_gericht 
        and positie*/

{
  graph g[MAXN];
  int i,j,end;


  memcpy(g,staticg,aantal_toppen*sizeof(graph));

  for (i=j=0; i<aantal_toppen; i++) /* j is de positie in edgelist */
    while (g[i])
      { end=FIRSTBIT(g[i]); DELELEMENT(g+i,end);  DELELEMENT(g+end,i); 
        edgelist[j][0]=i; edgelist[j][1]=end; 
	positie[i][end]=positie[end][i]=j;
	j++;
      }
}



void fill_edgelist_edgeorb_order()
  /* writes the edges in gg into the list in a way that for i in 1...numberedges we always have
     that the starting point of edge i has the minimum degree of all vertices in the graph you get
     when removing all edges i+1...aantal_bogen from gg.
     The global variable positie is not assigned any values.
 */
{
  int last_positie,i, buffer, olddeg, newdeg, beste, start, end;
  int list[MAXN][MAXN], listlen[MAXN]={0}; /* list[i] contains the vertices with degree i at that moment */
  int toppositie[MAXN], bufferdeg[MAXN]; /* toppositie[i] is de positie van top i in de lijst */
  int buren[MAXN]; /* de buren van de top waaraan gewerkt wordt */
  graph g[MAXN];

  if (nodegbound) { fill_edgelist_edgeorb(); return; }

  memcpy(g,staticg,aantal_toppen*sizeof(graph));

    for (i=0; i<aantal_toppen; i++) { buffer=bufferdeg[i]=deg[i]; list[buffer][listlen[buffer]]=i; 
                                      toppositie[i]=listlen[buffer]; (listlen[buffer])++; }


    last_positie=aantal_bogen-1;
    while (last_positie>=0)
      {
	for (beste=1;listlen[beste]==0; beste++);
	(listlen[beste])--; /* zo wordt hij ook verwijderd */
	start=list[beste][listlen[beste]];
	for (i=0; i<beste; i++) /* alle buren opslaan*/
	  { end=buren[i]=FIRSTBIT(g[start]); DELELEMENT(g+start,end); }

	for (i=0; i<beste; i++) /* alle buren */
	  { 
	    end=buren[i];
	    edgelist[last_positie][0]=start; edgelist[last_positie][1]=end; 
	    last_positie--;
	    DELELEMENT(g+end,start);
	    /* de buur verhuizen: */
	    olddeg=bufferdeg[end];
	    (bufferdeg[end])--;
	    newdeg=bufferdeg[end];
	    /* uit de oude lijst verwijderen */
	    if (listlen[olddeg]==1) listlen[olddeg]=0;
	    else { (listlen[olddeg])--;
	    buffer= list[olddeg][listlen[olddeg]];
	    list[olddeg][toppositie[end]]=buffer;
	    toppositie[buffer]=toppositie[end]; }
	    /* tot de nieuwe toevoegen */
	    if (newdeg)
	      { list[newdeg][listlen[newdeg]]=end;
	        toppositie[end]=listlen[newdeg];
	        (listlen[newdeg])++;
	      }
	  }

      }

}


void fill_edgelist_edgeorb_final()
     /* writes the edges in gg into the list in a lexicographic way and initializes is_gericht 
        and positie*/

{
  graph g[MAXN];
  int i,j,end;


  memcpy(g,staticg,aantal_toppen*sizeof(graph));

  for (i=0, j=aantal_gerichte_bogen; i<aantal_toppen; i++) /* j is de positie in edgelist */
    while (g[i])
      { end=FIRSTBIT(g[i]); DELELEMENT(g+i,end);  DELELEMENT(g+end,i); 
      if (!is_gericht[i][end])
	{ edgelist_final[j][0]=i; edgelist_final[j][1]=end;
	j++; }
      }
}



void fill_edgelist_edgeorb_order_final()
     /* works like fill_edgelist_edgeorb_order -- only that in this case some edges are already assigned a direction.
      They are not added to the list. The first entry is at position aantal_gerichte_bogen, so that the list
      can be used like the normal edgelist for directing in the trivial case.
      The global variable positie is not assigned any values.
*/

{
  int last_positie,i, buffer, olddeg, newdeg, beste, start, end;
  int list[MAXN][MAXN], listlen[MAXN]={0}; /* list[i] contains the vertices with degree i at that moment */
  int toppositie[MAXN], bufferdeg[MAXN]; /* toppositie[i] is de positie van top i in de lijst */
  int buren[MAXN]; /* de buren van de top waaraan gewerkt wordt */
  graph g[MAXN];

  if (nodegbound || (aantal_bogen-aantal_gerichte_bogen<7)) { fill_edgelist_edgeorb_final(); return; }

  memcpy(g,staticg,aantal_toppen*sizeof(graph));

    for (i=0; i<aantal_toppen; i++) { buffer=bufferdeg[i]=deg[i]; list[buffer][listlen[buffer]]=i; 
                                      toppositie[i]=listlen[buffer]; (listlen[buffer])++; }

    last_positie=aantal_bogen-1;
    while (last_positie>=aantal_gerichte_bogen)
      {
	for (beste=1;listlen[beste]==0; beste++);
	(listlen[beste])--; /* zo wordt hij ook verwijderd */
	start=list[beste][listlen[beste]];
	for (i=0; i<beste; i++) /* alle buren opslaan*/
	  { end=buren[i]=FIRSTBIT(g[start]); DELELEMENT(g+start,end); }

	for (i=0; i<beste; i++) /* alle buren */
	  if (!is_gericht[start][end=buren[i]])
	  { 
	    edgelist_final[last_positie][0]=start; edgelist_final[last_positie][1]=end; 
	    last_positie--;
	    DELELEMENT(g+end,start);
	    /* de buur verhuizen: */
	    olddeg=bufferdeg[end];
	    (bufferdeg[end])--;
	    newdeg=bufferdeg[end];
	    /* uit de oude lijst verwijderen */
	    if (listlen[olddeg]==1) listlen[olddeg]=0;
	    else { (listlen[olddeg])--;
	    buffer= list[olddeg][listlen[olddeg]];
	    list[olddeg][toppositie[end]]=buffer;
	    toppositie[buffer]=toppositie[end]; }
	    /* tot de nieuwe toevoegen */
	    if (newdeg)
	      { list[newdeg][listlen[newdeg]]=end;
	        toppositie[end]=listlen[newdeg];
	        (listlen[newdeg])++;
	      }
	  }

      }

}



void  mark_components(int graaf[][MAXN],int adj[],int aantal_toppen,int number[])
     /* Does a bfs on the vertices of graaf and makes that for every vertex i the
	value number[i] is the number of the smallest vertex in the component. 
     */
{
  int i,j,min, lijst[MAX_BOGEN], buffer, b;
  int *run, *end;

  RESETMARKS

  for (i=0;i<aantal_toppen;i++) 
    if (UNMARKED(i))
      { number[i]=i;
      if (adj[i])
	{ lijst[0]=min=i;
          MARK(i);
	  run=lijst; end=lijst+1;
	  while (run<end)
	    { buffer= *run;
	    for (j=0;j<adj[buffer];j++) 
	      if (UNMARKED(graaf[buffer][j]))
		{ b=graaf[buffer][j];
		MARK(b);
	        *end=b; end++;
		number[b]=min;
		}
	    run++;
	    }
	}
      }
}


#define ADDEDGE(g,ad,a,b) {(g)[a][(ad)[a]]=b; (ad)[a]++; g[b][(ad)[b]]=a; (ad)[b]++; }

void mark_orbitnumbers_edgelist(int number[], int *specialexists)
     /* Computes the orbits of not yet directed edges by assigning the same number to edges in the 
	same orbit. If an edge is on position i in the edgelist, number[i] is its orbitnumer.
	The orbitnumber is always the number of the smallest edge in the orbit -- this is used 
	in other routines !!
	If an edge is found that is stabilized by all elements of the group, but the endpoints
	are not, it's index is written 
	to specialexists. If no such edge exists, *specialexists = -1;
	In the first case the numbers in number[] are undefined. 
	It is assumed that generators and edgelist are up to date and especially that no undirected
	edges are mapped on directed ones by the given permutations.
     */

{ int i, j, good_special,pos2;
  BOOG boog; 
  int graaf[MAX_BOGEN][MAXN], adj[MAX_BOGEN];

  for (i=0; i<aantal_bogen; i++) { adj[i]=0; }

  for (i=0; i<aantal_bogen; i++)
    {COPYBOOG(boog,edgelist[i]);
    if (!is_gericht[boog[0]][boog[1]])
      { 
      good_special=1;
      for (j=0; j<number_of_generators; j++)
	{ 
	  pos2=POSBILD(boog,j);
	  if (i!=pos2) /* het beeld is verschillend */
	  { good_special=0;
	    ADDEDGE(graaf,adj,i,pos2);
	  }
	else { if (orbits[boog[0]]!=orbits[boog[1]]) good_special=0; }
	}
      if (good_special==1) { *specialexists=i; return; }
      }
    }

  mark_components(graaf,adj,aantal_bogen,number);

  *specialexists=-1;
}


void mark_orbitnumbers_edgelist_first(int number[], int *specialexists, int *completelyfixededge)
     /*
       The same as mark_orbitnumbers_edgelist -- only that *completelyfixededge is also filled in
       -- comments see get_orbit_first()
     */
{ int i, j, good_special,edge_fixed, pos2;
  BOOG boog; 
  int graaf[MAX_BOGEN][MAXN], adj[MAX_BOGEN];

  *completelyfixededge = -1;

  for (i=0; i<aantal_bogen; i++) adj[i]=0;

  for (i=0; i<aantal_bogen; i++)
    {COPYBOOG(boog,edgelist[i]);
    if (!is_gericht[boog[0]][boog[1]])
      { 
      good_special=edge_fixed=1;
      for (j=0; j<number_of_generators; j++)
	{ 
	pos2=POSBILD(boog,j);
	if (i!=pos2) /* het beeld is verschillend */
	  { good_special=edge_fixed=0; 
	    ADDEDGE(graaf,adj,i,pos2);
	  }
	else { if (orbits[boog[0]]!=orbits[boog[1]]) good_special=0; }
	}
      if (good_special) { *specialexists=i; if (*completelyfixededge >=0) return; }
      else { if (edge_fixed) *completelyfixededge=i; }
    }
    }

  mark_components(graaf,adj,aantal_bogen,number);

  *specialexists=-1;
}




void mark_orbitnumbers(int number[], BOOG list_of_dir_edges[], int listlength)
     /* Computes the orbits of the DIRECTED edges in list_of_dir_edges.
        All edges in list_of_dir_edges[] are interpreted as going from ...[][0] to ...[][1]*/

{ int i, j, pos2;
  BOOG boog; 
  int positie[2*MAX_BOGEN][2*MAX_BOGEN];
  int graaf[MAX_BOGEN][MAXN], adj[MAX_BOGEN];

  for (i=0; i<listlength; i++) 
    { positie[list_of_dir_edges[i][0]][list_of_dir_edges[i][1]]=i; adj[i]=0; }

  for (i=0; i<listlength; i++)
    { COPYBOOG(boog,list_of_dir_edges[i]);
      for (j=0; j<number_of_generators; j++)
	{ 
	pos2=POSBILD(boog,j);
	if (pos2!=i)  /* het beeld is verschillend */
	  { ADDEDGE(graaf,adj,i,pos2);
	  }
	}
    }
  mark_components(graaf,adj,listlength,number);
}


void mark_orbitnumbers_only_directed(int number[], BOOG list_of_edges[], int listlength)
     /* Computes the orbits of those edges in list_of_edges that are already directed --
	but interprets them as undirected. The subset of already directed edges in the list 
	must be closed under the automorphisms given in the global variable generators...*/

{ int i, j,a,b,pos2;
  BOOG boog; 
  int positie[2*MAX_BOGEN][2*MAX_BOGEN];
  int graaf[MAX_BOGEN][MAXN], adj[MAX_BOGEN];

  for (i=0; i<listlength; i++)
    { adj[i]=0;
      a=list_of_edges[i][0]; b=list_of_edges[i][1];
      if (is_gericht[a][b])
	{ positie[a][b]=positie[b][a]=i; number[i]=1; }
      else number[i]= 0;
    }

  for (i=0; i<listlength; i++)
    if (number[i])
      { 
      COPYBOOG(boog,list_of_edges[i]);
      for (j=0; j<number_of_generators; j++)
	{ 
	pos2=POSBILD(boog,j);
	if (pos2!=i)  /* het beeld is verschillend */
	  { ADDEDGE(graaf,adj,i,pos2);
	  }
	}
    }

  mark_components(graaf,adj,listlength,number);

}

void mark_orbitnumbers_only_candidates(int number[], BOOG list_of_edges[], int listlength, int candidate[])
     /* Computes the orbits of those edges in list_of_edges that are already directed --
	but interprets them as undirected. The subset of already directed edges in the list 
	must be closed under the automorphisms given in the global variable generators...*/

{ int i, j,a,b,pos2;
  BOOG boog; 
  int positie[2*MAX_BOGEN][2*MAX_BOGEN];
  int graaf[MAX_BOGEN][MAXN], adj[MAX_BOGEN];

  for (i=0; i<listlength; i++)
    { adj[i]=0;
      if (candidate[i])
	{ a=list_of_edges[i][0]; b=list_of_edges[i][1];
	  positie[a][b]=positie[b][a]=i; 
	}
    }

  for (i=0; i<listlength; i++)
    if (candidate[i])
      { 
      COPYBOOG(boog,list_of_edges[i]);
      for (j=0; j<number_of_generators; j++)
	{ 
	pos2=POSBILD(boog,j);
	if (pos2!=i)  /* het beeld is verschillend */
	  { ADDEDGE(graaf,adj,i,pos2); 
	  }
	}
    }
  mark_components(graaf,adj,listlength,number);

}





void get_orbit(BOOG kleinste_orbit[], int *orbitsize, int *biggest_orbit)
     /* writes an orbit of still undirected edges of minimal size into kleinste_orbit. It is assumed
	that generators and edgelist are up to date. 
	Orbits of size one with the endpoints in different orbits are not considered -- fixing
	this edge won't help at all, because it is already fixed by every automorphism.

*/
{
  /* the numbers of edges are quite small, so it is not necessary to use more difficult union-set-
     algorithms */
  int i,j, number[MAX_BOGEN],aantallen[MAX_BOGEN];
  int orbit_met_een,min,minorb,max;

  mark_orbitnumbers_edgelist(number, &orbit_met_een);

  if (orbit_met_een>=0) { COPYBOOG(kleinste_orbit[0],edgelist[orbit_met_een]); 
                          *orbitsize=1; return; }

  /* else */
  //for (i=0; i<aantal_bogen; i++) aantallen[i]=0;
  // Omdat number[i] altijd de kleinste in zijn orbit is, kan aantallen ook als volgt op 0 gezet worden:
  for (i=0; i<aantal_bogen; i++)
    { aantallen[i]=0;
    // if (!is_gericht[edgelist[i][0]][edgelist[i][1]]) al gerichte bogen hebben hun eigen number en
    // dus orbitsize aantallen 1 achteraf en worden dus toch al niet gekozen
    (aantallen[number[i]])++; 
    }


  for (i=max=0, min=INT_MAX, minorb=-1 ; i<aantal_bogen; i++)
    { 
      if ((aantallen[i]>1) && (aantallen[i]<min)) /* als er een goede met orbitsize 1 is, kom je hier niet
						     omdat die al als special is gekozen */
	{ min=aantallen[i]; minorb=i; }
      if ((aantallen[i])>max) max=aantallen[i];
    }
  *biggest_orbit=max;

  /* it is possible that no orbit was found -- that is: in spite of the fact that the group is NOT
     trivial, it does act trivially on the set of still undirected edges: */

  if (minorb== -1) { *orbitsize=0; return; }

  /* nu worden de bogen uit de gekozen orbit in de 
     lijst kleinste_orbit geschreven */


   for (i=j=0; i<aantal_bogen; i++) 
     { if (number[i]==minorb) 
       { COPYBOOG(kleinste_orbit[j],edgelist[i]); j++; }
     } 
   *orbitsize=j;
}



void get_orbit_first(BOOG kleinste_orbit[], int *orbitsize, int *fixedgeindex, int *biggest_orbit)
     /* works like get_orbit(0), but may only be called for the initial still undirected
	graph. *fixedgeindex is the index of an edge in the global edgelist where both
	endpoints are fixed under the group (if such an edge exists -- otherwise it is -1.
     */

{
  /* the numbers of edges are quite small, so it is not necessary to use more difficult union-set-
     algorithms */
  int i,j, number[MAX_BOGEN],aantallen[MAX_BOGEN];
  int orbit_met_een,min,minorb, max;

  *fixedgeindex= -1;

  mark_orbitnumbers_edgelist_first(number, &orbit_met_een, fixedgeindex);

  if (orbit_met_een>=0) 
    { COPYBOOG(kleinste_orbit[0],edgelist[orbit_met_een]); 
      *biggest_orbit=MAXPAR_ORBSIZE+1; *orbitsize=1; return; }
 
  /* else */
  for (i=0; i<aantal_bogen; i++) aantallen[i]=0;
  for (i=0; i<aantal_bogen; i++) (aantallen[number[i]])++;
  //if (!is_gericht[edgelist[i][0]][edgelist[i][1]]) (aantallen[number[i]])++;


  for (i=max=0, min=INT_MAX, minorb= -1 ; i<aantal_bogen; i++)
    { 
      if ((aantallen[i]) && (aantallen[i]<min) && 
	  ((aantallen[i]>1) || (orbits[edgelist[i][0]]==orbits[edgelist[i][1]]))) 
	{ min=aantallen[i]; minorb=i; }
      if ((aantallen[i])>max) max=aantallen[i];
    }

  *biggest_orbit=max;


  /* it is possible that no orbit was found -- that is: in spite of the fact that the group is NOT
     trivial, it does act trivially on the set of still undirected edges: */

  if (minorb == -1) { *orbitsize=0; return; }
  /* the group only permutes isolated vertices */

  /* nu worden de bogen uit de gekozen orbit in de 
     lijst kleinste_orbit geschreven */


   for (i=j=0; i<aantal_bogen; i++) 
     { if (number[i]==minorb) 
       { COPYBOOG(kleinste_orbit[j],edgelist[i]); j++; }
     } 
   *orbitsize=j;
}

void trynextstep()
/* This routine is called when a whole orbit has been labelled and new orbits have to be computed.
   If orbitsgiven != NULL it contains the actual information about the new orbits. */
{  int i, orbitsize, dummy, k,j,complete,start,end, biggest_orbit=0;
  BOOG kleinste_orbit[MAX_BOGEN+1];
  BOOG al_gericht[MAX_BOGEN+1];
  int inorbit[MAXN][MAXN], numberinorbit[MAXN];
  graph buffergraph[MAXN];

  if (aantal_gerichte_bogen==aantal_bogen) { MAYBEPROCESS; WRITEUP(); return; }

  if (!nodegbound)
    {
      for (i=0, complete=1; i<aantal_bogen; i++) 
	{ if (!is_gericht[edgelist[i][0]][edgelist[i][1]] && 
	      !virtual_gericht[edgelist[i][0]][edgelist[i][1]]) { complete=0; break; }
	}
 
      if (complete)  
	{ if (direct_output==0) { WRITEUP(); } 
	  else
	    { 
	      memcpy(buffergraph,workg,aantal_toppen*sizeof(graph));
              for (i=0; i<aantal_bogen; i++) 
		{ start=edgelist[i][0]; end=edgelist[i][1];
		  if (!is_gericht[start][end])
		  { if (virtual_gericht[start][end]==1) 
		      DELELEMENT(workg+end,start);
		    else  DELELEMENT(workg+start,end); 
		  }
		}
	      MAYBEPROCESS;
	      WRITEUP();
	      memcpy(workg,buffergraph,aantal_toppen*sizeof(graph));
	    }
	return; }
    }

  if (!group_up_to_date)
    { /* we need only the group -- no canonical numbering */
         number_of_generators=0;
	 memcpy(bufferlab,lab[nextstep_depth],aantal_toppen*sizeof(nvector));
	 memcpy(bufferptn,ptn[nextstep_depth],aantal_toppen*sizeof(nvector));
	 nauty(workg,bufferlab,bufferptn,NILSET,orbits,&options_directed,&stats,workspace,100*MAXN,1,aantal_toppen,NULL);
	 group_up_to_date=1;
  }

  nextstep_depth++;

  if (stats.numorbits==aantal_toppen) /* triviale groep */
    { 
      fill_edgelist_edgeorb_order_final();
      laatstepositie=edgelist_final+aantal_bogen-1;
      trivlabels_init(edgelist_final+aantal_gerichte_bogen);
      nextstep_depth--;
      return;
    }
  /* else: de groep is jammer genoeg niet triviaal */
  get_orbit(kleinste_orbit,&orbitsize,&biggest_orbit);


  if (orbitsize==0) /* group acts trivially on remaining edges */
    { 
      fill_edgelist_edgeorb_order_final();
      laatstepositie=edgelist_final+aantal_bogen-1;
      trivlabels_init(edgelist_final+aantal_gerichte_bogen);
      nextstep_depth--;
      return;
    }
  /* else */ 

  memcpy(colour[nextstep_depth],orbits,aantal_toppen*sizeof(nvector));
  for (i=0;i<aantal_toppen;i++) numberinorbit[i]=0;
  for (i=0;i<aantal_toppen;i++) { dummy=orbits[i]; 
                                  inorbit[dummy][numberinorbit[dummy]]=i;
				  (numberinorbit[dummy])++;
                                }
  for (i=k=0;i<aantal_toppen;i++) 
    { for (j=0;j<numberinorbit[i]-1; j++) { ptn[nextstep_depth][k]=1; lab[nextstep_depth][k]=inorbit[i][j]; k++; } 
      if (numberinorbit[i]) { ptn[nextstep_depth][k]=0; lab[nextstep_depth][k]=inorbit[i][j]; k++; }
    }



  /* maybe at least the whole orbit is forced... */
  if (!nodegbound)
    {
      for (i=0, complete=1; i<orbitsize; i++) 
	{ if (!virtual_gericht[kleinste_orbit[i][0]][kleinste_orbit[i][1]]) { complete=0; break; }
	}
      
      if (complete) 
	{ aantal_gerichte_bogen+=orbitsize;
	  for (i=0; i<orbitsize; i++) 
		{ start=kleinste_orbit[i][0]; end=kleinste_orbit[i][1];
		  saturated[start]++; saturated[end]++;
		  is_gericht[start][end]=is_gericht[end][start]=1;
		  if (virtual_gericht[start][end]==1) 
		      { DELELEMENT(workg+end,start);
			outdeg_free[start]--; indeg_free[end]--;
		      }
		    else  
		      { DELELEMENT(workg+start,end); 
			outdeg_free[end]--; indeg_free[start]--;
		      }
		  }
	  group_up_to_date=0;
	  // I guess the group is still correct -- but better write up a formal proof before
	  // assuming this in the program...

	  trynextstep();

	  aantal_gerichte_bogen-=orbitsize;
	  for (i=0; i<orbitsize; i++) 
		{ start=kleinste_orbit[i][0]; end=kleinste_orbit[i][1];
		  saturated[start]--; saturated[end]--;
		  is_gericht[start][end]=is_gericht[end][start]=0;
		  if (virtual_gericht[start][end]==1) 
		    { ADDELEMENT(workg+end,start);
		    outdeg_free[start]++; indeg_free[end]++;
		    }
		  else  
		    { ADDELEMENT(workg+start,end); 
		    outdeg_free[end]++; indeg_free[start]++;
		    }
		}
	  group_up_to_date=0;
	  nextstep_depth--;
	return; }
    }

  nontrivlabels(kleinste_orbit,0,orbitsize,al_gericht,biggest_orbit);
  nextstep_depth--;
  return;
}

void trynextstep_par() 
/* This routine is called when a whole orbit has been labelled and new orbits have to be computed.
   If orbitsgiven != NULL it contains the actual information about the new orbits. */
{  int i, orbitsize, dummy, k,j,biggest_orbit;
  BOOG kleinste_orbit[MAX_BOGEN+1];
  int inorbit[MAXN][MAXN], numberinorbit[MAXN];

  if (aantal_gerichte_bogen==aantal_bogen) { MAYBEPROCESS; WRITEUP(); return; }

  if (!group_up_to_date)
    { /* we need only the group -- no canonical numbering */
         number_of_generators=0;
	 memcpy(bufferlab,lab[nextstep_depth],aantal_toppen*sizeof(nvector));
	 memcpy(bufferptn,ptn[nextstep_depth],aantal_toppen*sizeof(nvector));
	 nauty(workg,bufferlab,bufferptn,NILSET,orbits,&options_directed,&stats,workspace,100*MAXN,1,aantal_toppen,NULL);
	 group_up_to_date=1;
  }

  nextstep_depth++;

  if (stats.numorbits==aantal_toppen) /* triviale groep */
    { 
      fill_edgelist_edgeorb_order_final();
      laatstepositie=edgelist_final+aantal_bogen-1;
      trivlabels_init(edgelist_final+aantal_gerichte_bogen);
      nextstep_depth--;
      return;
    }
  /* else: de groep is jammer genoeg niet triviaal */
  get_orbit(kleinste_orbit,&orbitsize,&biggest_orbit);


  if (orbitsize==0) /* group acts trivially on remaining edges */
    { 
      fill_edgelist_edgeorb_order_final();
      laatstepositie=edgelist_final+aantal_bogen-1;
      trivlabels_init(edgelist_final+aantal_gerichte_bogen);
      nextstep_depth--;
      return;
    }
  /* else */ 

  memcpy(colour[nextstep_depth],orbits,aantal_toppen*sizeof(nvector));
  for (i=0;i<aantal_toppen;i++) numberinorbit[i]=0;
  for (i=0;i<aantal_toppen;i++) { dummy=orbits[i]; 
                                  inorbit[dummy][numberinorbit[dummy]]=i;
				  (numberinorbit[dummy])++;
                                }
  for (i=k=0;i<aantal_toppen;i++) 
    { for (j=0;j<numberinorbit[i]-1; j++) { ptn[nextstep_depth][k]=1; lab[nextstep_depth][k]=inorbit[i][j]; k++; } 
      if (numberinorbit[i]) { ptn[nextstep_depth][k]=0; lab[nextstep_depth][k]=inorbit[i][j]; k++; }
    }

  parallel_orbit_labelling(kleinste_orbit, orbitsize);
  nextstep_depth--;
  return;
}



void compute_extensions(BOOG kleinste_orbit[], int orbitsize, BOOG extensionlist[], int *number_of_extensions)
     /* Computes the number of extensions by computing the orbits on the possibilities to direct the 
	not yet directed edges in kleinste_orbit. To be precise: First a list of directed versions of the
	not yet directed edges is made (2 directed for each undirected) and then the orbits are computed
	and one directed edge for every version is filled in. 
	Already at this stage indeg and outdeg are considered, so that directed edges that would violate
	these conditions aren't considered.
*/
{
BOOG bufferlist[2*MAX_BOGEN];
BOOG candidatelist[MAX_BOGEN];
int number[2*MAX_BOGEN], number_of_candidates;
int buffersize,n_ext;
int i,start,end, problem, sum, maxsum;


 /* the last edge should always be one with minimal QUALITY -- good edges early. See Definition of QUALITY
    in the header. */


 maxsum=INT_MAX;
 number_of_candidates=0;

 for (i=0; i<orbitsize; i++)
   { start=kleinste_orbit[i][0]; end=kleinste_orbit[i][1];
   if (!is_gericht[start][end]) { candidatelist[number_of_candidates][0]=start; 
                                  candidatelist[number_of_candidates][1]=end;
				  number_of_candidates++; }
   else /* one that is already directed */
     { 
     if (ISELEMENT(workg+start,end)) sum=QUALITY_P1(start,end); else sum=QUALITY_P1(end,start);
     /* this sum can be increased by one when the new edge is added */
       if (sum<maxsum) maxsum=sum;
     }
   }
 
 for (i=buffersize=0; i<number_of_candidates; i++)
     { start=candidatelist[i][0]; end=candidatelist[i][1];
       /* the tests are done in the loop to detect problems with unassignable edges also if they
	  do not fulfill the canonicity criteria */
       problem=1; /* maybe the edge cannot be directed at all !! */
       if (outdeg_free[start] && indeg_free[end]) 
       { sum=QUALITY_P2(start,end);
	 if (sum<=maxsum) {bufferlist[buffersize][0]=start; bufferlist[buffersize][1]=end; buffersize++;} 
         problem=0;}
     if (outdeg_free[end] && indeg_free[start]) 
       { sum=QUALITY_P2(end,start);
	 if (sum<=maxsum) {bufferlist[buffersize][0]=end; bufferlist[buffersize][1]=start; buffersize++;} 
         problem=0;}
     if (problem) { *number_of_extensions=0; return; }
     }

 if (buffersize>1)
  { 
    if (!group_up_to_date) 
      { number_of_generators=0;
        memcpy(bufferlab,lab[nextstep_depth],aantal_toppen*sizeof(nvector));
	memcpy(bufferptn,ptn[nextstep_depth],aantal_toppen*sizeof(nvector));
	nauty(workg,bufferlab,bufferptn,NILSET,orbits,&options_directed_canon,&stats,\
	      workspace, 100*MAXN,1,aantal_toppen,canong);
	group_up_to_date=1;
      }
     mark_orbitnumbers(number,bufferlist,buffersize);
     for (i=n_ext=0; i<buffersize; i++) 
     if (number[i]==i) { COPYBOOG(extensionlist[n_ext],bufferlist[i]);
                         n_ext++; }
     *number_of_extensions=n_ext;
   }
 else
   { if (buffersize)
     { COPYBOOG(extensionlist[0],bufferlist[0]);
       *number_of_extensions=1; }
   else *number_of_extensions=0;
   }
 return;

}


int getexpensivequality(graph x,graph y)
{
int j;
int qual=0;

 FORALLELEMENTS(x,j) qual+=(saturated[j]<<4)-indeg_free[j];
 qual = qual<<6;
 FORALLELEMENTS(y,j) qual+=(saturated[j]<<4)-indeg_free[j];

 return qual;
}

int is_canonical_edge(BOOG list[],int last_positie)
/* Checks whether the directed edge x->y is canonical -- that is: 
   choose one with minimal QUALITY 
   and amongst them one with biggest (deg[end]+maxoutdeg-outdeg_free[end])<<6+indeg_free[start] 
   the smallest among all lexicographic pairs canon_number(a),canon_number(b)
   with a->b or b-> a a directed edge in list[]. Then check whether
   a->b is in the same orbit as the edge of this smallest pair.
   The group must be up to date.
 */
{
  int i, minx=INT_MAX, miny=INT_MAX, a,b, which=INT_MAX;
  /* Some of the initializations are not necessary -- just to get rid of warnings
     if compiled with -Wall */
  int canonnumber[MAXN],x,y;
  int number[MAX_BOGEN];
  int candidate[MAX_BOGEN];
  int referencesum, sum, gotacandidate, endquality, lq, eq, expensivequality;

  x=list[last_positie][0]; y=list[last_positie][1];

  referencesum=QUALITY(x,y);
  endquality=((deg[y]+maxoutdeg-outdeg_free[y])<<6)-indeg_free[x];
  expensivequality= -1;
  candidate[last_positie]=1;

  gotacandidate=0;
  for (i=0; i<last_positie; i++)
    { a=list[i][0]; b=list[i][1];
      sum=QUALITY(a,b);
      if (sum<referencesum) return 0;
      if (sum==referencesum)
	{ 
	  
	  lq=((deg[b]+maxoutdeg-outdeg_free[b])<<6)-indeg_free[a]; 
	  if (lq>endquality) return 0;
	  if (lq==endquality) 
	    { /* OK -- nog een poging */
	      if (expensivequality== -1) 
		expensivequality=getexpensivequality(workg[x],workg[y]);
	      eq=getexpensivequality(workg[a],workg[b]);
	      if (eq<expensivequality) return 0;
	      else 
		if (eq==expensivequality)
		  { candidate[i]=1; gotacandidate=1; }
		else candidate[i]=0;
	    }
	  else candidate[i]=0;
	}
      else candidate[i]=0;
    }
   
  if (!gotacandidate) return 1; /* Only the edge itself */

  number_of_generators=0;
  memcpy(bufferlab,lab[nextstep_depth],aantal_toppen*sizeof(nvector));
  memcpy(bufferptn,ptn[nextstep_depth],aantal_toppen*sizeof(nvector));
  nauty(workg,bufferlab,bufferptn,NILSET,orbits,&options_directed_canon,&stats,\
	workspace, 100*MAXN,1,aantal_toppen,canong);
  group_up_to_date=1;
  
  for (i=0;i<aantal_toppen;i++) canonnumber[bufferlab[i]]=i;

  for (i=0; i<=last_positie; i++)
    if (candidate[i])
    {
      a=canonnumber[list[i][0]]; b=canonnumber[list[i][1]];
      if ((a<minx) || ((a==minx) && (b<miny))) { minx=a; miny=b; which=i; }
    }

  if (which==last_positie) return 1;
   
  mark_orbitnumbers_only_candidates(number, list, last_positie+1, candidate);

  if (number[which]==number[last_positie]) return 1; else return 0;

}


void force_edges(int top, int all_out, BOOG forcelist[], int *listlen, int *problem)
     /* directs forced edges, *problem is set to 1 if the partial
	direction cannot be completed according to the rules 

	all_out==1 betekent dat de bogen naar buiten gericht moeten worden ==0
	naar binnen.

	This routine assumes that if indeg==maxindeg there is still enough room
	to force all the rest outgoing (and reverse). So it must be tested in advance
	that maxindeg+maxoutdeg<=deg[i] for all vertices i!

*/

{ int j,i,end;
  int list[2*MAXN+2], richting[3*MAXN];

  end=0;

  if (virtual_indeg[top]+virtual_outdeg[top]==deg[top]) return;

  if (all_out) /* all the rest must be directed outwards */
    { 
      FORALLELEMENTS(workg[top],j)
	{
	  if (!is_gericht[top][j] && !virtual_gericht[top][j]) /* a new edge to be virtually directed */
	    { if (virtual_indeg[j]==maxindeg) { *problem=1; return; }
	    virtual_indeg[j]++; virtual_outdeg[top]++;
	    virtual_gericht[top][j]=1; virtual_gericht[j][top]=2;
	    forcelist[*listlen][0]=top; forcelist[*listlen][1]=j; (*listlen)++;
	    if (virtual_indeg[j]==maxindeg) { list[end]=j; richting[end]=1; end++; }
	    /* top kan maar 1 keer door virtual_indeg[top]==maxindeg toegevoegd worden */
	    }
	}
      if (*problem) return;
    }
  else
    /* all the rest must be directed inwards */
    { 
      FORALLELEMENTS(workg[top],j)
	{
	  if (!is_gericht[top][j] && !virtual_gericht[top][j]) /* a new edge to be virtually directed */
	    { if (virtual_outdeg[j]==maxoutdeg) { *problem=1; return; }
	    virtual_indeg[top]++; virtual_outdeg[j]++;
	    virtual_gericht[top][j]=2; virtual_gericht[j][top]=1;
	    forcelist[*listlen][1]=top; forcelist[*listlen][0]=j; (*listlen)++;
	    if (virtual_outdeg[j]==maxoutdeg) { list[end]=j; richting[end]=0; end++; }
	    /* top kan maar 1 keer door virtual_outdeg[top]==maxoutdeg toegevoegd worden */
	    }
	}
      
    }

  
  for (i=0; i<end; i++)
    { top=list[i];
    all_out=richting[i];
    if (virtual_indeg[top]+virtual_outdeg[top]<deg[top])
      {
	if (all_out) /* all the rest must be directed outwards */
	  { 
	    FORALLELEMENTS(workg[top],j)
	      {
		if (!is_gericht[top][j] && !virtual_gericht[top][j]) /* a new edge to be virtually directed */
		  { if (virtual_indeg[j]==maxindeg) { *problem=1; return; }
		  virtual_indeg[j]++; virtual_outdeg[top]++;
		  virtual_gericht[top][j]=1; virtual_gericht[j][top]=2;
		  forcelist[*listlen][0]=top; forcelist[*listlen][1]=j; (*listlen)++;
		  if (virtual_indeg[j]==maxindeg) { list[end]=j; richting[end]=1; end++; }
		  /* top kan maar 1 keer door virtual_indeg[top]==maxindeg toegevoegd worden */
		  }
	      }
	  if (*problem) return;
	}
	else
       /* all the rest must be directed inwards */
	{
	  FORALLELEMENTS(workg[top],j)
	    { 
	      if (!is_gericht[top][j] && !virtual_gericht[top][j]) /* a new edge to be virtually directed */
		{ if (virtual_outdeg[j]==maxoutdeg) { *problem=1; return; }
		virtual_indeg[top]++; virtual_outdeg[j]++;
		virtual_gericht[top][j]=2; virtual_gericht[j][top]=1;
		forcelist[*listlen][1]=top; forcelist[*listlen][0]=j; (*listlen)++;
		if (virtual_outdeg[j]==maxoutdeg) { list[end]=j; richting[end]=0; end++; }
		/* top kan maar 1 keer door virtual_outdeg[top]==maxoutdeg toegevoegd worden */
		}
	    }

	}
    }
    }
}

int allemaal_doubles_mogelijk(BOOG edgelist[],int orbitsize, int marker[], int *endlist)
     /* returns 1 if all edges in edgelist that are not yet directed can be made double edges,
      otherwise 0. marker [i] is 1 if edgelist[i] can be made double and is still to be directed
      and 0 otherwise.

      The fields fixing that they are directed are also filled in.

      maybe I should rewrite the program in a way to always have the list of undirected 
      edges available -- unfortunately that takes so little time that it is hardly worth it... */
{
  int i, start, end;

  *endlist= -1;
  for (i=0; i<orbitsize; i++)
  { start=edgelist[i][0]; end=edgelist[i][1];
  if (!is_gericht[start][end])
    {
      if (!nodegbound)
	{
	  if (!(double_free[start] && double_free[end])) { return 0; }
	  if (!((virtual_indeg[start]<maxindeg) && (virtual_indeg[end]<maxindeg) &&
		(virtual_outdeg[start]<maxoutdeg) && (virtual_outdeg[end]<maxoutdeg))) { return 0; }
	  if (virtual_gericht[start][end]) { return 0; } /* dat is altijd in maar een richting... */
	}
      /* else -- kann dubbel worden */
      marker[i]=1;
      *endlist=i; /* de laatste positie waar er iets veranderd werd */
      remaining_doubles--;
      aantal_gerichte_bogen++;
      is_gericht[start][end]=is_gericht[end][start]=1;
      /* a double edge can never be forced */
      virtual_indeg[start]++; virtual_indeg[end]++; 
      virtual_outdeg[start]++; virtual_outdeg[end]++; 
      saturated[start]+=2; saturated[end]+=2;
      /* note that this is only compatible with the symmetry group, because all
	 edges in the orbit are handled at the same time -- otherwise it would 
	 destroy symmetry -- not yet handled edges and edges that are decided to stay
	 double will be different for the quality criterion if they were identic before. */
      (indeg_free[end])--; (outdeg_free[start])--;
      (indeg_free[start])--; (outdeg_free[end])--;
      double_free[start]--; double_free[end]--;
    }
  else marker[i]=0;
  }

  return 1;

}

void reset_doubles(BOOG edgelist[], int marklist[], int marklistend)

{
  int i, start, end;

  for (i=0; i<=marklistend; i++)
    if (marklist[i])
      { start=edgelist[i][0]; end=edgelist[i][1];
        is_gericht[start][end]=is_gericht[end][start]=0;
	virtual_indeg[start]--; virtual_indeg[end]--; 
	virtual_outdeg[start]--; virtual_outdeg[end]--; 
	saturated[start]-=2; saturated[end]-=2;
      /* note that this is only compatible with the symmetry group, because all
	 edges in the orbit are handled at the same time -- otherwise it would 
	 destroy symmetry -- not yet handled edges and edges that are decided to stay
	 double will be different for the quality criterion if they were identic before. */
	remaining_doubles++;
	(indeg_free[end])++; (outdeg_free[start])++;
	(indeg_free[start])++; (outdeg_free[end])++;
	double_free[start]++; double_free[end]++;
	aantal_gerichte_bogen--;
    }

  return;

}

/*************************DO_EXTENSIONS_PAR********************************/

void do_extensions_par(BOOG orbit[],int positie,int numberin[],int numberout[], int numberdouble[],unsigned int op)
{
  int start, end;

  if (positie<0) { parops[number_parops]=op; number_parops++; return; }

  start=orbit[positie][0]; end=orbit[positie][1];

  if ((indeg_free[start]-numberin[start]>0) && (outdeg_free[end]-numberout[end]>0))
    // incoming possible
    { SETOP(op,positie,2); 
      numberin[start]++; numberout[end]++;
      do_extensions_par(orbit,positie-1,numberin,numberout,numberdouble,op);
      numberin[start]--; numberout[end]--;
    }

  if (double_allowed && (indeg_free[start]-numberin[start]>0) && (indeg_free[end]-numberin[end]>0) &&
      (outdeg_free[start]-numberout[start]>0) && (outdeg_free[end]-numberout[end]>0) &&
      (double_free[start]-numberdouble[start]>0) && (double_free[end]-numberdouble[end]>0))
    // double possible
    { SETOP(op,positie,1); 
      numberin[start]++; numberin[end]++; numberout[start]++; numberout[end]++;
      numberdouble[start]++; numberdouble[end]++;
      do_extensions_par(orbit,positie-1,numberin,numberout,numberdouble,op);
      numberin[start]--; numberin[end]--; numberout[start]--; numberout[end]--;
      numberdouble[start]--; numberdouble[end]--;
    }


  if ((indeg_free[end]-numberin[end]>0) && (outdeg_free[start]-numberout[start]>0))
    // outgoing possible
    { SETOP(op,positie,0); 
      numberin[end]++; numberout[start]++;
      do_extensions_par(orbit,positie-1,numberin,numberout,numberdouble,op);
      numberin[end]--; numberout[start]--;
    }

  return;


}

/*************************COMPUTE_PAR_EXTENSIONS********************************/

int compute_par_extensions(BOOG orbit[], int orbitsize)
// orbitsize must be at least 1
{
  int numberin[MAXN], numberout[MAXN], numberdouble[MAXN], i;
  unsigned int op;

  for (i=0; i<aantal_toppen; i++) numberin[i]=numberout[i]=numberdouble[i]=0;


  number_parops=0;
  op=0U;

  do_extensions_par(orbit,orbitsize-1,numberin,numberout,numberdouble,op);

  return number_parops;
}

unsigned int par_image(unsigned int orig, permutation aut[], BOOG orbit[], int orbitsize, int inv[][MAXN])
{
  int i, type, start, end, positie;
  unsigned int image;

  image=0U;

  for (i=0;i<orbitsize;i++)
    { type=GETTYPE(orig,i);
      start=aut[orbit[i][0]]; end=aut[orbit[i][1]];
      positie=inv[start][end];
      if (positie<0)
	{ positie= -1-positie; type=2-type; }
      SETOP(image,positie,type);
    }

  return image;

}

#define FIND_PAROP(op,pos) {int min,max; min=0; max=num_extensions-1; pos=(min+max)/2;\
    while (min!=max) { if (op<parops[pos]) min=pos+1; else max=pos; pos=(min+max)/2; }}


int compute_par_orbits(BOOG edgeorbit[],int orbitsize,int num_extensions, unsigned int nonequivextensions[])
{
  int root[MAXPAROPS], i, j, pos1, pos2, orbits, buf;
  unsigned int buffer;
  int inv[MAXN][MAXN]; 
  // inv[i][j] is the position of edge {i,j} in edgeorbit. If inv[i][j]<0 this means that in fact {j,i} is in
  // edgeorbit -- and that on position -(inv[i][j])-1 -- the one is necessary as -0 = 0

  for (i=0;i<orbitsize;i++) { inv[edgeorbit[i][0]][edgeorbit[i][1]]=i; inv[edgeorbit[i][1]][edgeorbit[i][0]]= -i-1; }

  for (i=0; i<num_extensions; i++) root[i]=i;

  for (i=0; i<num_extensions; i++)
    for (j=0; j<number_of_generators; j++)
      { 
	buffer=par_image(parops[i],generators[j],edgeorbit,orbitsize,inv);
	FIND_PAROP(buffer,pos1);
	while (root[pos1]!=pos1) pos1=root[pos1];
	pos2=i;
	while ((buf=root[pos2])!=pos2) { root[pos2]=pos1; pos2=buf; }
	root[pos2]=pos1;
      }

  for (orbits=i=0; i<num_extensions; i++)
    { if (root[i]==i) { nonequivextensions[orbits]=parops[i]; orbits++; } }

  return orbits;

}


void parallel_orbit_labelling(BOOG edge_orbit[], int orbitsize)
{
  int num_extensions, i,j,type, start, end;
  unsigned int nonequivextensions[MAXPAROPS], buf;

  
  if (!group_up_to_date) 
    { number_of_generators=0;
      memcpy(bufferlab,lab[nextstep_depth],aantal_toppen*sizeof(nvector));
      memcpy(bufferptn,ptn[nextstep_depth],aantal_toppen*sizeof(nvector));
      nauty(workg,bufferlab,bufferptn,NILSET,orbits,&options_directed_canon,&stats, \
	    workspace, 100*MAXN,1,aantal_toppen,canong);
      group_up_to_date=1;
      }

  num_extensions=compute_par_extensions(edge_orbit, orbitsize);

  num_extensions=compute_par_orbits(edge_orbit,orbitsize,num_extensions,nonequivextensions);

  


  aantal_gerichte_bogen += orbitsize;
  for (j=0;j<orbitsize;j++) 
    { start=edge_orbit[j][0]; end=edge_orbit[j][1];
      saturated[start]++; saturated[end]++; 
      is_gericht[start][end]=is_gericht[end][start]=1;
    }

  for (i=0;i<num_extensions;i++)
    { buf=nonequivextensions[i];
      for (j=0;j<orbitsize;j++)
	{ type=GETTYPE(buf,j);
	  start=edge_orbit[j][0]; end=edge_orbit[j][1];
	  
	  if (type==0) // outgoing
	    { outdeg_free[start]--; indeg_free[end]--;
	      DELELEMENT(workg+end,start);
	    }
	  else
	  if (type==2) // incoming
	    { outdeg_free[end]--; indeg_free[start]--;
	      DELELEMENT(workg+start,end);
	    }
	  else // (type==1)
	    { outdeg_free[end]--; outdeg_free[start]--; remaining_doubles--;
	      indeg_free[end]--; indeg_free[start]--;
	    }
	}

      group_up_to_date=0;
      trynextstep_par(); 
      // don't switch back to general routine as some fields are not filled in
      // properly -- e.g. virtual degrees and saturated
      
      for (j=0;j<orbitsize;j++)
	{ type=GETTYPE(buf,j);
	  start=edge_orbit[j][0]; end=edge_orbit[j][1];
	  
	  if (type==0) // outgoing
	    { outdeg_free[start]++; indeg_free[end]++;
	      ADDELEMENT(workg+end,start);
	    }
	  else
	  if (type==2) // incoming
	    { outdeg_free[end]++; indeg_free[start]++;
	      ADDELEMENT(workg+start,end);
	    }
	  else // (type==1)
	    { outdeg_free[end]++; outdeg_free[start]++; remaining_doubles++;
	      indeg_free[end]++; indeg_free[start]++;
	    }
	}
      group_up_to_date=0;
    }
  
  aantal_gerichte_bogen -= orbitsize;
  for (j=0;j<orbitsize;j++) 
    { start=edge_orbit[j][0]; end=edge_orbit[j][1];
      saturated[start]--; saturated[end]--; 
      is_gericht[start][end]=is_gericht[end][start]=0;
    }
  
  return;

}


void nontrivlabels(BOOG kleinste_orbit[], int done_in_orbit, int orbitsize, BOOG al_gericht[],
		   int maxorbit)
     /* geeft een richting aan de bogen in het geval van automorphismen */
{
  int i,j,start,end, number_of_extensions, problem, forcelistlen;
  BOOG extensionlist[2*MAX_BOGEN], forcelist[MAX_BOGEN];
  int marklist[MAX_BOGEN], marklistend;


  if ((done_in_orbit==0) && // just to make sure the other tests are just done once for each orbit
      (orbitsize>1) && (orbitsize<=SWITCHPAR_ORBSIZE) && (maxorbit<=MAXPAR_ORBSIZE))
    { parallel_orbit_labelling(kleinste_orbit, orbitsize); return; }

 compute_extensions(kleinste_orbit, orbitsize, extensionlist,&number_of_extensions);

 for (i=0; i<number_of_extensions; i++)
   { 
       start=extensionlist[i][0]; end=extensionlist[i][1];
       if (((virtual_outdeg[start]<maxoutdeg) && (virtual_indeg[end]<maxindeg)) || (virtual_gericht[start][end]==1))
        { problem=forcelistlen=0;
	  outdeg_free[start]--; indeg_free[end]--;
	  saturated[start]++; saturated[end]++;
	  is_gericht[start][end]=is_gericht[end][start]=1;
	  DELELEMENT(workg+end,start);
	  group_up_to_date=0;
	  aantal_gerichte_bogen++;
	  al_gericht[done_in_orbit][0]=start; al_gericht[done_in_orbit][1]=end;
	  if (!virtual_gericht[start][end]) 
	    { 
	      virtual_outdeg[start]++; virtual_indeg[end]++;
	      if (virtual_outdeg[start]==maxoutdeg) 
		force_edges(start, 0, forcelist, &forcelistlen, &problem);
	      if (!problem && (virtual_indeg[end]==maxindeg) ) 
		force_edges(end, 1, forcelist, &forcelistlen, &problem); }

	  if (!problem && ((done_in_orbit==0) || is_canonical_edge(al_gericht,done_in_orbit)))
	 /* if is_canonical_edge() is called, done_in_orbit>0, so nauty was just called and the 
	    group is up to date */
	    { 
	      if (done_in_orbit+1 < orbitsize)
		nontrivlabels(kleinste_orbit, done_in_orbit+1, orbitsize,al_gericht,maxorbit);
	    else
	      trynextstep();
	    }
	  outdeg_free[start]++; indeg_free[end]++;
	  saturated[start]--; saturated[end]--;
	  is_gericht[start][end]=is_gericht[end][start]=0;
	  ADDELEMENT(workg+end,start);
	  group_up_to_date=0;
	  aantal_gerichte_bogen--;
	  if (!virtual_gericht[start][end]) 
	    { virtual_outdeg[start]--; virtual_indeg[end]--;
	    for (j=0; j<forcelistlen;j++)
	      {
		start=forcelist[j][0]; end=forcelist[j][1];
		virtual_gericht[start][end]=virtual_gericht[end][start]=0;
		virtual_outdeg[start]--; virtual_indeg[end]--;
	      }
	    }
	}
   } /* einde for-loop */

 /* conventie: dubbel edges hebben de hoogste prioriteit om verwijdert te worden. Dat betekent 
    dat ze als laatste in een orbit die opgevuld moet worden toegevoegd worden */

 if ((remaining_doubles >= orbitsize-done_in_orbit) && (orbitsize>done_in_orbit))
   {
     if (allemaal_doubles_mogelijk(kleinste_orbit,orbitsize,marklist,&marklistend))
       { trynextstep(); }
     reset_doubles(kleinste_orbit,marklist,marklistend);
   }
}

int connected(graph g[], int aantal_toppen)
{
  int i, list[MAXN], length, j;
  graph reached, dummy;

  reached=(graph)0;
  ADDELEMENT(&reached,0);

  list[0]=0; list[1]= -1;
  length=1;

  for (i=0; i<length; i++)
    { dummy= g[list[i]] & ~reached;
      reached |= g[list[i]];
      FORALLELEMENTS(dummy,j)
	{ list[length]=j; length++; }
    }

  if (reached==ALLMASK(aantal_toppen))
    { return 1; }
  else 
    { return 0; }
}

int test_possible(graph globalg[],int globaldeg[], int n, int m, int min_direct_deg)
     /* Test of sommige deelgrafen aan |E|<= min_direct_deg*aantal_toppen voldoet.

	Geeft 1 terug als er misschien een manier bestaat om de bogen een richting 
	toe te kennen en 0 als zo'n manier zeker niet bestaat.

     */

{
  graph g[MAXN];
  int deg[MAXN];
  int i,j, removed, grens, top, m0;
  int *runp, *endp;
  int list[2*MAXN]; /* de lijst van toppen. In het begin staan ze er allemaal in. Als er een top al bekeken is
		       en dan word de graad min_direct_deg word hij er opnieuw toegevoegd -- maar dat kan
		       maar 1 keer gebeuren. 
		       Het probleem is te voorkomen dat een top 2 keer in het gedeelte van de lijst staat
		       waaraan nog gewerkt wordt. Dan kan het gebeuren dat de top 2 keer verwijderd word.

		    */

  grens=n*min_direct_deg;

  if (aantal_bogen>grens) return 0;

  memcpy(g,globalg,n*sizeof(graph));
  memcpy(deg,globaldeg,n*sizeof(int));

  for (i=0; i<n; i++) list[i]=i;

  removed=0;
  for (runp=list, endp=list+n; runp<endp; runp++)
    { top=*runp;
      if (deg[top]<=min_direct_deg)
      { 
	grens-=min_direct_deg;
	m-=deg[top];
	deg[top]=0;
	removed=1;
	if (m>grens) { return 0; }
        FORALLELEMENTS(g[top],j)
	  { 
	    DELELEMENT(g+j,top);
	    deg[j]--;
	    /* top<i betekent: de eerste keer dat hij in de lijst stond was al vroeger. En en tweede keer
	       kan de top er nog niet staan omdat de volgende voorwaarde deg[j]==min_direct_deg garandeerd 
	       dat de graad net voldoende gedaald is */
	    if ((top<(runp-list)) && (deg[j]==min_direct_deg)) { *endp=j; endp++; }
	  }
      }
    }


  if (connected(globalg, n) && (!removed)) return 1;

  /* Nu met BFS kijken of er een componente met te veel bogen is: */

  RESETMARKS

  for (i=0;i<n;i++) 
    if (UNMARKED(i) && deg[i])
      { m0=deg[i];
	list[0]=i;
	MARK(i);
	for (runp=list, endp=list+1; runp<endp; runp++)
	  { top=*runp;
	  FORALLELEMENTS(g[top],j) 
	    { if (UNMARKED(j))
	      { MARK(j);
	        m0+=deg[j];
		*endp=j;
		endp++;
	      }
	    }
	  }
	if (min_direct_deg*(endp-list)*2<m0) { return 0;}
      }


  return 1;
}



void waterclusters (graph g[], int n)
{ 
  int i, j, k, orbitsize, start, end, fixed_edge, dummy, maxgraphdeg, biggest_orbit;
  BOOG kleinste_orbit[MAX_BOGEN+1];  /* het is mogelijk dat er maar 1 orbit is */
  BOOG al_gericht[MAX_BOGEN+1];
  int inorbit[MAXN][MAXN], numberinorbit[MAXN];
  //static int counter=0;

  init_waterplugin(n);

  if (n<3) { fprintf(stderr,"Come on -- start with at least 3 vertices! \n"); exit(0); }

  memcpy(workg,g,n*sizeof(graph));

  aantal_toppen=n;

  for (i=0;i<aantal_toppen;i++) {indeg_free[i]=maxindeg; outdeg_free[i]=maxoutdeg;
                                 virtual_indeg[i]=virtual_outdeg[i]=0; }

  number_of_generators=0;
  nauty(staticg,bufferlab,bufferptn,NILSET,orbits,&options,&stats,workspace,100*MAXN,1,n,NULL);
  group_up_to_date=1;

  maxgraphdeg=0;
  for (i=aantal_bogen=0; i<n ; i++) { deg[i]=POPCOUNT(g[i]); if (deg[i]>maxgraphdeg) maxgraphdeg=deg[i];
                                      aantal_bogen+=deg[i]; }
  aantal_bogen /= 2;

   /* als maxindeg==maxoutdeg is gegarandeerd dat voor elke graaf die gegenereerd word ook een manier
     bestaat om richtingen toe te kennen zonder de voorwaarden te schenden */
  if (maxgraphdeg>2*mingerichtdeg) 
    { if (test_possible(g,deg,n,aantal_bogen,mingerichtdeg)==0) return; }

  //if ((maxgraphdeg<=maxindeg) && (maxgraphdeg<=maxoutdeg)) nodegbound=1; else nodegbound=0;

  if (aantal_bogen==0) { addnumber=1; MAYBEPROCESS; WRITEUP(); return; }

  if (double_allowed) 
    max_doubles=remaining_doubles=watermaxedges-aantal_bogen; else max_doubles=remaining_doubles=0;
  if (remaining_doubles)
    { /* deg[i] is altijd <= maxindeg+maxoutdeg */
      for (i=0;i<aantal_toppen;i++) double_free[i]=watermaxdeg-deg[i];
    }
  else 
    for (i=0;i<aantal_toppen;i++) double_free[i]=0;

  nextstep_depth=0;

  if (stats.numorbits!=n)
    { fill_edgelist_edgeorb();
      get_orbit_first(kleinste_orbit,&orbitsize,&fixed_edge,&biggest_orbit); 
    }

  /* Now is (orbitsize==0) if the group acts trivially on edges 
     -- only interesting in case of disconnected graphs */

  if ((stats.numorbits==n) || (orbitsize==0)) /* trviale groep */
    { aantal_grafen_met_triv_group++; 
      fill_edgelist_edgeorb_order();
      laatstepositie=edgelist+aantal_bogen-1;
      if (maxoutdeg==maxindeg) /* then every valid graph for one direction is also valid with
				  the directions reversed */
	{
	  addnumber=2;
	  start=edgelist[0][0]; end=edgelist[0][1];
	  /* is_gericht[][] isn't used in the directing routine */
	  (outdeg_free[start])--; (indeg_free[end])--;
	  //virtual_outdeg[start]++; virtual_indeg[end]++;
	  DELELEMENT(workg+end,start);
	  aantal_gerichte_bogen=1;
	  trivlabels_init(edgelist+1);
	  aantal_gerichte_bogen=0;
	  ADDELEMENT(workg+end,start);
	  (outdeg_free[start])++; (indeg_free[end])++;

	  if (remaining_doubles && double_free[start] && double_free[end])
	    {
	      addnumber=1; 
	      (outdeg_free[start])--; (indeg_free[end])--;
	      (outdeg_free[end])--; (indeg_free[start])--;
	      double_free[start]--; double_free[end]--;
	      aantal_gerichte_bogen=1;
	      remaining_doubles--;
	      trivlabels_init(edgelist+1);
	      remaining_doubles++;
	      double_free[start]++; double_free[end]++;
	      (outdeg_free[start])++; (indeg_free[end])++;
	      (outdeg_free[end])++; (indeg_free[start])++;
	      aantal_gerichte_bogen=0;
	    }
	}
      else
	{
	  addnumber=1;
	  aantal_gerichte_bogen=0;
	  trivlabels_init(edgelist);
	  aantal_gerichte_bogen=0;
	}
      return;
    }
  /* else: de groep is jammer genoeg niet triviaal */


  /* first: write good beginning colours for nauty */
  memcpy(colour[0],orbits,aantal_toppen*sizeof(nvector));
  for (i=0;i<aantal_toppen;i++) numberinorbit[i]=0;
  for (i=0;i<aantal_toppen;i++) { dummy=orbits[i]; 
                                  inorbit[dummy][numberinorbit[dummy]]=i;
				  (numberinorbit[dummy])++;
                                }
  for (i=k=0;i<aantal_toppen;i++) 
    { for (j=0;j<numberinorbit[i]-1; j++) { ptn[0][k]=1; lab[0][k]=inorbit[i][j]; k++; } 
      if (numberinorbit[i]) { ptn[0][k]=0; lab[0][k]=inorbit[i][j]; k++; }
    }

  /* else */
  if ((fixed_edge<0)|| (maxindeg!=maxoutdeg)) 
    { addnumber=1; nontrivlabels(kleinste_orbit,0,orbitsize,al_gericht, biggest_orbit); return; }
  /* else */
  /* one edge can already be directed -- but because it is fixed anyway, the group is still correct. */

  addnumber=2;
  start=edgelist[fixed_edge][0]; end=edgelist[fixed_edge][1];
  (outdeg_free[start])--; (indeg_free[end])--;
  virtual_outdeg[start]++; virtual_indeg[end]++;
  saturated[start]=saturated[end]=1;
  is_gericht[start][end]=is_gericht[end][start]=1;
  DELELEMENT(workg+end,start);
  aantal_gerichte_bogen=1;
  nontrivlabels(kleinste_orbit,0,orbitsize,al_gericht, biggest_orbit);
  aantal_gerichte_bogen=0;
  ADDELEMENT(workg+end,start);
  is_gericht[start][end]=is_gericht[end][start]=0;
  virtual_outdeg[start]--; virtual_indeg[end]--;
  (outdeg_free[start])++; (indeg_free[end])++;
  saturated[start]=saturated[end]=0; 

  if (remaining_doubles && double_free[start] && double_free[end])
    { 
      addnumber=1;
      (outdeg_free[start])--; (indeg_free[end])--;
      (outdeg_free[end])--; (indeg_free[start])--;
      double_free[start]--; double_free[end]--;
      virtual_outdeg[start]++; virtual_outdeg[end]++;
      virtual_indeg[start]++; virtual_indeg[end]++;
      saturated[start]=saturated[end]=2;
      is_gericht[start][end]=is_gericht[end][start]=1;
      aantal_gerichte_bogen=1;
      remaining_doubles--;
      nontrivlabels(kleinste_orbit,0,orbitsize,al_gericht, biggest_orbit);
      remaining_doubles++;
      aantal_gerichte_bogen=0;
      is_gericht[start][end]=is_gericht[end][start]=0;
      double_free[start]++; double_free[end]++;
      virtual_outdeg[start]--; virtual_outdeg[end]--;
      virtual_indeg[start]--; virtual_indeg[end]--;
      (outdeg_free[start])++; (indeg_free[end])++;
      (outdeg_free[end])++; (indeg_free[start])++;
      saturated[start]=saturated[end]=0; 
    }


  return;
}








/******************************DIRECT_EDGES********************************/

void direct_edges(void) /* graph workg[] and int aantal_toppen, aantal_bogen are global */

{ int i, maxedges, minrestriction, maxdeg, free_vertices, maxgraphdeg, regular;

  regular=1;
  for (i=maxgraphdeg=0; i<aantal_toppen; i++)
    { if (deg[i]>maxgraphdeg) maxgraphdeg=deg[i];
      if (deg[i]!=deg[0]) regular=0;
    }

  if ((maxgraphdeg<=maxindeg) && (maxgraphdeg<=maxoutdeg)) nodegbound=1; else nodegbound=0;

  if (maxindeg<=maxoutdeg) { maxedges=aantal_toppen*maxindeg; minrestriction=maxindeg;}
  else { maxedges=aantal_toppen*maxoutdeg; minrestriction=maxoutdeg; }

  if (maxedges<aantal_bogen) return;
 /* grafen met meer bogen kunnen niet gericht worden, zodat ze nog aan de 
    voorwaarden voldoen */

  maxdeg=maxindeg+maxoutdeg;
  if (!nodegbound)
    {
      for (i=free_vertices=0;i<aantal_toppen;i++) 
	{ if (deg[i]>maxdeg) return;
	  if (deg[i]<maxdeg) free_vertices++;
	}
    }
  else free_vertices=aantal_toppen;

  group_up_to_date=0;

  // when to use what is just some heuristic -- a more elaborate one could be helpful
    if ((!regular) &&
	(nodegbound || (maxindeg<=2 && maxoutdeg<=2) || (free_vertices>((2*aantal_toppen)/3))))
    { //waterclusteruse++; 
      waterclusters (staticg, aantal_toppen); return; }

  //water_v_use++;

  memcpy(workg,staticg,sizeof(graph)*aantal_toppen);

  for (i=0; i<aantal_toppen; i++) 
    { indeg[i]=outdeg[i]=0;
      if (double_allowed) 
	{ double_free[i]=maxdirectdeg-deg[i];
	  if (minrestriction<double_free[i]) double_free[i]=minrestriction;
	}
      else double_free[i]=0;
    }

  aantal_gerichte_bogen=0;

  if (double_allowed) 
    max_doubles=remaining_doubles=maxedges-aantal_bogen; else max_doubles=remaining_doubles=0;


  number_of_generators=0;
  nauty(workg,bufferlab,bufferptn,NILSET,orbits,&options,&stats,workspace,100*MAXN,1,aantal_toppen,NULL);

  if (stats.numorbits==aantal_toppen) 
    { 
      direct_all_triv();
    }
  else 
    { 
      for (i=0; i<aantal_toppen; i++) { tobedirected[i]=deg[i]; }

      direct_all_nontriv(); }


  return;

}

void init_allocated_fields()
{
  int i;

  size_root=1000; root_op=malloc((size_t)size_root*sizeof(int)); 
  if (root_op==NULL) 
    { fprintf(stderr,"Can't allocate %d items for root_op in the beginning -- exiting.\n",size_root);
      exit(0); }


  operations=malloc((size_t)4096);
  if (operations==NULL) 
    { fprintf(stderr,"Can't allocate initial memory for operations -- exiting\n"); 
      exit(1); }
  size_operations=4096;

  for (i=0; i<MAXN; i++)
    {
      remember_operations[i]=malloc((size_t)4096);
      if (remember_operations[i]==NULL) 
	{ fprintf(stderr,"Can't allocate initial memory for operations -- exiting\n"); 
	  exit(1); }
      remember_size[i]=4096;
    }
  return;
}



/******************************MAIN********************************/

int main(int argc, char *argv[])

// reads from stdin 
{ 

  int i, m, zaehlen=0;
  unsigned char *code=NULL;
  int codelaenge;
  int multicode=0, g6code=1;
  long long int last=0LL;

  if (sizeof(long long int)<8) 
    { 
      fprintf(stderr,"This may cause problems with the hashing function for large degree -- exit().\n");
      exit(1); 
    }

  for (i=1; i<argc; i++)
    {
      if (argv[i][0]=='i') maxindeg=atoi(argv[i]+1);
      else  if (argv[i][0]=='o') maxoutdeg=atoi(argv[i]+1);
      else  if (argv[i][0]=='T') direct_output=1; 
	else  if (argv[i][0]=='C') direct_output=2;
	  else  if (argv[i][0]=='B') direct_output=3;
	    else  if (argv[i][0]=='S') double_allowed=0;
	      else  if (argv[i][0]=='m') { g6code=0; multicode=1; }
      else usage(argv[0]);
    }

#ifdef PROCESS
  if (direct_output==0) direct_output=2;
#endif

  init_allocated_fields();

  //if (maxindeg==MAXN && maxoutdeg==MAXN) nodegbound=1; else nodegbound=0;
  // wordt in direct_edges individueel vastgelegd
  maxdirectdeg=maxindeg+maxoutdeg;

  init_nauty_options();

  while((g6code && (readg(stdin,staticg,1,&m,&aantal_toppen) != NULL)) ||
	(multicode && (lese_multicode(&code, &codelaenge, stdin) != EOF)))
    { 
      zaehlen++; 
#ifdef PROCESS
      dg_nin = zaehlen;
#endif

      if (aantal_toppen>=MAXN)
	{ fprintf(stderr,"At most %d vertices possible -- exiting.\n",MAXN-1);
	  exit(1); }

      if (multicode) decode_to_nauty(code,codelaenge,staticg,deg);
      else /* g6code */ init_for_g6(staticg,aantal_toppen,deg);

      // direct_edges(staticg, aantal_toppen);
      direct_edges();   // BDM
      last=aantal_gerichte_grafen;



  }

  fprintf(stderr,"Number of directed graphs: %lld \n",aantal_gerichte_grafen);
#ifdef SUMMARY
  SUMMARY();
#endif

  //fprintf(stderr,"Used edge routines %u times.\n",waterclusteruse);
  //fprintf(stderr,"Used vertex routines %u times.\n",water_v_use);

  return 0;
}
