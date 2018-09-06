/******************************************************************************
 *                                                                            *
 * This is the main file for traces() version 2.1, which is included into     *
 *   nauty() version 2.6.                                                     *
 *                                                                            *
 *   nauty is Copyright (1984-2016) Brendan McKay.  All rights reserved.      *
 *   Traces is Copyright Adolfo Piperno, 2008-2016.  All rights reserved.     *
 *   See the file COPYRIGHT for the details of the software license.          *
 *                                                                            *
 *   CHANGE HISTORY                                                           *
 *       28-Dec-12 : final changes for version 2.0                            *
 *       20-Jan-13 : add code for ^C catching in Traces                       *
 *       29-Mar-13 : bug correction in automorphism mode                      *
 *       02-Apr-13 : add preprocessing                                        *
 *       21-May-13 : bug correction (coloured lists)                          *
 *       29-Jun-13 : bug correction (coloured lists and cycles)               *
 *       07-Dec-13 : bug correction in automorphism mode (wrong group size    *
 *                   due to randomness in Schreier-Sims orbit computation)    *
 *                   bug correction (discrete initial partition)              *
 *       15-Feb-14 : CPUDEFS removed (already declared in gtools.h)           *
 *       01-Sep-15 : add weighted edges (not active)                          *
 *       28-Jan-16 : version ready for nauty and Traces v.2.6 distribution    *
 *       12-Jul-16 : bug correction (reaching degree 2 vertices)              *
 *****************************************************************************/

#include "traces.h"

#ifdef NAUTY_IN_MAGMA
#include "cleanup.e"
#endif

#define SORT_OF_SORT 2
#define SORT_NAME sort2ints
#define SORT_TYPE1 int
#define SORT_TYPE2 int
#include "sorttemplates.c"

typedef struct weightwhere {
    int weight;
    int *ref;
} weightwhere;

#define SORT_OF_SORT 2
#define SORT_NAME sortweights
#undef SORT_TYPE2
#define SORT_TYPE1 int
#define SORT_TYPE2 weightwhere
#include "sorttemplates.c"

#define NAUTY_ABORTED (-11)
#define NAUTY_KILLED (-12)

typedef struct Candidate {
    boolean sortedlab;
    int *invlab;
    int *lab;
    int code;
    int do_it;
    int indnum;
    int name;
    int vertex;
    struct Candidate *next;
    struct searchtrie *stnode;
    unsigned int firstsingcode;
    unsigned int pathsingcode;
    unsigned int singcode;
} Candidate;

typedef struct Partition {
    int *cls;
    int *inv;
    int active;
    int cells;
    int code;
} Partition;

typedef struct trielist {
    struct searchtrie *triearray;
    struct trielist *prev;
    struct trielist *next;
} trielist;

typedef struct TracesVars {
    char digstring[25];
    double autchk;
    double expaths;
    double schreier1;
    double schreier2;
    double schreier3;
    int augmented_cells;
    int build_autom;
    int *currorbit;
    int *orbits;
    int answ;
    int brkstpcount;
    int compstage;
    int cand_level;
    int canlist;
    int digits;
    int expathlength;
    int firstpathlength;
    int fromlevel;
    int group_level;
    int indivend;
    int indivstart;
    int indiv_vtx;
    int lastcell;
    int lastlev;
    int lev_of_lastauto;
    int levelfromCS0;
    int linelgth;
    int mark;
    int treemark;
    int autmark;
    int markcell1;
    int markcell2;
    int maxdeg;
    int maxtreelevel;
    int maxspineorblevel;
    int mindeg;
    int name;
    struct searchtrie *gotonode;
    struct searchtrie *newgotonode;
    struct searchtrie *newst_stage1;
    int newindex;
    int nextlevel;
    int nfix;
    int finalnumcells;
    int permInd;
    int preprocessed;
    int samepref;
    int specialgens;
    int stackmark;
    int steps;
    int strategy;
    trielist *strielist;
    int strienext;
    int tcell_sz;
    int tcell;
    int tcellevel;
    int tcellexpath_sz;
    int tcellexpath;
    int tolevel_tl;
    int tolevel;
    int treedepth;
    int trienext;
    int triepos;
    TracesOptions *options;
    TracesStats *stats;
    unsigned int singlongcode;
    sparsegraph *graph;
    sparsegraph *cangraph;
    sparsegraph *input_graph;
    int conta0;
    int conta1;
    int conta2;
    int conta3;
    int conta4;
    int conta5;
    int conta6;
    int conta7;
    int contatc;
} TracesVars;

typedef struct TracesInfo {
    boolean autofound;
    boolean deg_one;
    boolean first_matching;
    boolean regular;
    boolean exitfromref;
    boolean identitygroup;
    boolean minimalinorbits;
    boolean thegraphisparse;
    boolean thegrouphaschanged;
    boolean thereisnextlevel;
    boolean useTempOrbits1;
    boolean useTempOrbits2;
} TracesInfo;

typedef struct TracesSpine {
    boolean thetracexists;
    Candidate *listend;
    Candidate *liststart;
    int ccend;
    int ccstart;
    int listcounter;
    int stpend;
    int stpstart;
    int tgtcell;
    int tgtend;
    int tgtfrom;
    int tgtpos;
    int tgtsize;
    int trcend;
    int trcstart;
    int singend;
    int singstart;
    int updates;
    unsigned long keptcounter;
    unsigned long levelcounter;
    Partition *part;
    unsigned int singcode;
} TracesSpine;

typedef struct trie {
    int value;
    struct trie *first_child;
    struct trie *next_sibling;
} trie;

typedef struct searchtrie {
    int index;
    int name;
    int vtx;
    int level;
    struct searchtrie *father;
    struct searchtrie *first_child;
    struct searchtrie *last_child;
    struct searchtrie *next_sibling;
    struct searchtrie *goes_to;
} searchtrie;

typedef struct pair {
    int arg;
    int val;
} pair;

typedef struct grph_strct {
    int *e;
    int *w;
    int d;
    boolean one;
} grph_strct;

typedef struct ExpPathInfo {
    int code;
    int cell;
    int info;
} ExpPathInfo;

static boolean traces_degree_refine(sparsegraph*, Candidate*, int, Partition*,
                                    struct TracesVars*, struct TracesInfo*, int, int*);
static int  traces_vertexclass_refine (int, int*, int*, Candidate*, Partition*, int*);
static int  traces_refine(Candidate*, int, Partition*,
                          struct TracesVars*, struct TracesInfo*, int, boolean);
static void traces_refine_notrace(Candidate*, int, Partition*,
                                  struct TracesVars*, struct TracesInfo*);
static void traces_refine_maketrie(Candidate*, int, Partition*,
                                   struct TracesVars*, struct TracesInfo*);
static int  traces_refine_comptrie(Candidate*, int, Partition*,
                                   struct TracesVars*, struct TracesInfo*);
static int  traces_refine_sametrace(Candidate*, int, Partition*,
                                    struct TracesVars*, struct TracesInfo*);
static int  traces_refine_refine(sparsegraph*, Candidate*, int, Partition*,
                                 struct TracesVars*, struct TracesInfo*);
static int  refine_tr_refine(Candidate*, int, Partition*,
                             struct TracesVars*, struct TracesInfo*);
static int given_gens(sparsegraph*, permnode*, int*, boolean);
static void quickSort(int*, int);
static struct Partition* NewPartition(int);
static struct Candidate* NewCandidate(int, Candidate**, int);
static void NewPartSpine(int, int);
static int FreeList(Candidate*, int);
static int FixBase(int*, struct TracesVars*, Candidate*, int, int);
static boolean FixedBase(int*, struct TracesVars*, Candidate*, int, int);
static void factorial(double*, int*, int);
static void factorial2(double*, int*, int);
static int CheckForAutomorphisms(Candidate*, Candidate*, struct TracesVars*, struct TracesInfo*, int, int, Partition*);
static int CheckForSingAutomorphisms(Candidate*, Partition*, Candidate*, struct TracesVars*, struct TracesInfo*, int, int);
static int CheckForMatching(Candidate*, Candidate*, Partition*, struct TracesVars*, struct TracesInfo*, int, int);
static void Individualize(Partition*, Candidate*, int, int, int, int);
static boolean TreeFyTwo(int, Candidate*, Candidate*, Partition*, int, struct TracesVars*, struct TracesInfo*);
static void ExperimentalStep(Partition*, Candidate*, struct TracesVars*, struct TracesInfo*, int, int);
static boolean TargetCell(Candidate*, Partition*, int, struct TracesVars*, int);
static boolean TargetCellFirstPath(Candidate*, Partition*, struct TracesVars*);
static int TargetCellExpPath(Candidate*, Partition*, struct TracesVars*);
static boolean TargetCellSmall(Candidate*, Partition*, int, struct TracesVars*, int);
static boolean TargetCellFirstPathSmall(Candidate*, Partition*, struct TracesVars*);
static int TargetCellExpPathSmall(Candidate*, Partition*, struct TracesVars*);
static boolean SelectNextLevel(int, struct TracesVars*, struct TracesInfo*);
static void CopyCand(Candidate*, Candidate*, int, int*, int*);
static struct trie* trie_new(int, struct TracesVars*);
static struct trie* trie_make(trie*, int, int, struct TracesVars*);
static struct trie* trie_comp(trie*, int);
static void trie_dump(trie*);
static void trie_class(trie*, int*);
static void RemoveFromLevel(int, int, int, boolean);
static int CompStage0(Partition*, Partition*, Candidate*, Candidate*, int, int, struct TracesVars*, struct TracesInfo*);
static int CompStage1(Partition*, Partition*, Candidate*, Candidate*, int, int, struct TracesVars*, struct TracesInfo*);
static int CompStage2(Partition*, Partition*, Candidate*, Candidate*, int, int, struct TracesVars*, struct TracesInfo*);
static void grouporderplus(sparsegraph*, Candidate*, Partition*, permnode**, double*, int*, int, struct TracesVars*, struct TracesInfo*);
static boolean Prefix(Candidate*, Candidate*, int);
static boolean findperm(permnode*, int*, int);
static int spinelementorbsize(int*, int*, int, int);
static trielist* searchtrie_new(int, struct TracesVars*);
static searchtrie* searchtrie_make(Candidate*, Candidate*, int, struct TracesVars*);
static boolean lookup(searchtrie*);
static int* findcurrorbits(schreier*, int);
static int Preprocess(sparsegraph*, permnode**, Candidate*, int, Partition*, struct TracesVars*);
static int Preprocess_refine(sparsegraph*, permnode**, Candidate*, int, Partition*, struct TracesVars*);
static void MakeTree(int, int, sparsegraph*, int, struct TracesVars*, boolean);
static void MakeCanTree(int, sparsegraph*, int, Candidate*, Partition*, struct TracesVars*);
static int max(int, int);
static int min(int, int);
static void orbjoin_sp_perm(int*, int*, int*, int, int*);
static void orbjoin_sp_pair(int*, int*, int, int, int, int*);
static boolean isautom_sg_pair(graph*, int*, boolean, int, int, struct TracesVars*);
static void SetAutom(int, int, struct TracesVars*);
static void ResetAutom(int, int, struct TracesVars*);
static void PrintVect(int*, int, int, int);
static void putgraphplus_sg(FILE*, sparsegraph*, int);
static boolean VerifyId(int *p, int n);
static void PrintPartition(int*, int*, int, int, int);
static void Place(int, Candidate*, Partition*);
static int NonSingDeg(int, Candidate*, Partition*);
static int NonSingDegPlus1(Candidate*, Partition*, int, TracesVars*);
static void NonSingDegPlus2(Candidate*, Partition*, int, TracesVars*);
static void Edge_Delete(int, int, Candidate*, TracesVars*);
static boolean VerifyPart(Partition*, int, int);
static int VerifyPerm(int*, int,int);
static boolean VerifyCand(Candidate*, int, int);
static int FirstNeighbour(int, Candidate*, Partition*, int*, int, int*, int);
static int NextNeighbour(int, Candidate*, Partition*, int*, int, int*, int);
static sparsegraph* copy_sg_structure(sparsegraph*, sparsegraph*);
static void WeightCodes (int);
static void PrintWeightedGraph1(sparsegraph*, int, char[30]);
static void PrintWeightedGraph2(int n, char msg[30]);
static void MakeDiscrete(Partition*, int);
static void Complete(sparsegraph*, Candidate*, Partition*, int, TracesVars*, double*, int*,permnode**, int);
static void Allocate_Traces_Structures(int);
static void Allocate_refine_Structures(int);
static void Initialize_Traces_Variables(TracesVars*, TracesOptions*, TracesStats*, int*, sparsegraph*, sparsegraph*, int);
static void Initialize_Traces_Statistics (TracesStats*, int);
static void Initialize_Traces_Time_Variables (TracesVars*);
static int trie_classify(int, TracesVars*);
static int Check_degree_one(sparsegraph*, Candidate*, Partition*, int);
static void sort_Split_Array(int*, int);
static const unsigned int fuzz1[] = {037541, 061532, 005257, 026416};
static const unsigned int fuzz2[] = {006532, 070236, 035523, 062437};
static int Select_from_CStack(int*, int);
static void PrintBlissGraph(int);
static void CodeClassify(int, int, int);

#define FUZZ1(x) ((x) ^ fuzz1[(x)&3])
#define FUZZ2(x) ((x) ^ fuzz2[(x)&3])

#define MASHCOMM(l, i) ((l) + FUZZ1(i))
#define MASHNONCOMM(l, i) (FUZZ2(l) + (i))
#define MASH(l, i) ((((l) ^ 065435) + (i)) & 077777)
#define MASH1(l, i) ((l + (i*i)) & 077777)
#define CLEANUP(l) ((int)((l) % 0x7FFF))
#define SS(n, sing, plur)  (n), ((n) == 1?(sing):(plur))

#define SETMARK(Arr, Mrk) if (Mrk > (NAUTY_INFINITY-2)) { memset(Arr, 0, n*sizeof(int)); Mrk = 0; } Mrk++;

#define COPYNODE(W, V) { \
memcpy(W->lab, V->lab, n*sizeof(int)); \
memcpy(W->invlab, V->invlab, n*sizeof(int)); \
W->code = V->code; \
W->singcode = V->singcode; \
W->do_it = V->do_it; }

#define NEXTLINE fprintf(outfile, "\n");

#define PRINTCHAR(c) fprintf(outfile, "%s", c);

#define PRINTCAND(V, Lev) PRINTCHAR(" ") for (tmp=1; tmp<=Lev; tmp++) {fprintf(outfile, tv->digstring, V->lab[Spine[tmp].tgtpos]+labelorg);}

#define PRINTCANDF(V, Lev) { NEXTLINE for (tmp=1; tmp<=Lev; tmp++) {fprintf(outfile, "F%di", V->lab[Spine[tmp].tgtpos]+labelorg);} NEXTLINE }

#define PRINTCANDBIG(V, Lev) { PRINTCHAR(" ") \
for (tmp=1; tmp<=5; tmp++) {fprintf(outfile, tv->digstring, V->lab[Spine[tmp].tgtpos]+labelorg);} \
fprintf(outfile, "... "); \
for (tmp=Lev-4; tmp<=Lev; tmp++) {fprintf(outfile, tv->digstring, V->lab[Spine[tmp].tgtpos]+labelorg);} }

#define LINE(K, c) { PRINTCHAR(c) for (tmp=1; tmp<=K; tmp++) {fprintf(outfile, c);} }

#define TRACE_CHECK(Tr, Ind, Arg, End) { TracePos = Tr+Ind; \
if (newtrace) { \
*TracePos = Arg; \
} \
else { \
if (Ind < *End) { \
if (*TracePos != Arg) { \
if (*TracePos > Arg) { \
return FALSE; \
} \
else { \
*TracePos = Arg; \
newtrace = TRUE; \
} \
} \
} \
else { \
*TracePos = Arg; \
newtrace = TRUE; \
} \
} \
Ind++; }

#define SAMETRACE_CHECK(Tr, Ind, Arg, End) { TracePos = Tr+Ind; \
if (Ind < *End) { \
if (*TracePos != Arg) { \
return FALSE; \
} \
} \
else { \
return FALSE; \
} \
Ind++; }

#define NEWPARTSPINE(Lev) { if (Lev > 3) { \
Spine[Lev].part = malloc(sizeof(*(Spine[Lev].part))); \
if (Spine[Lev].part == NULL) { \
fprintf(ERRFILE, "\nError, memory not allocated.\n"); \
exit(1); \
} \
Spine[Lev].part->cls = Spine[Lev-3].part->cls; \
Spine[Lev].part->inv = Spine[Lev-3].part->inv; \
Spine[Lev-3].part->cls = Spine[Lev-3].part->inv = NULL; \
Spine[Lev].part->code = -1; \
Spine[Lev].part->cells = 0; \
} \
else { \
Spine[Lev].part = NewPartition(n); \
} }

#define FIND_SPLIT_CELLS SplInd = 0; \
for (j = 0; j < HitClsInd; j++) { \
ind1 = HitCls[j]; \
ElmHitCll[ind1] -= ind1; \
if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) { \
SplCls[SplInd++] = ind1; \
} \
}

#define FREEPART(Part) { if (Part) { \
if (Part->cls) free(Part->cls); \
if (Part->inv) free(Part->inv); \
free(Part); } \
}

#define FREECAND(Cand) { if (Cand) { \
if (Cand->lab) free(Cand->lab); \
if (Cand->invlab) free(Cand->invlab); \
free(Cand); \
} }

#define COPYPART(P, Q) { memcpy(P->cls, Q->cls, n*sizeof(int)); \
memcpy(P->inv, Q->inv, n*sizeof(int)); \
P->cells = Q->cells; \
P->code = Q->code; } \

#define ADDTONEXTLEVEL { if (SpineTL->listend) { \
(SpineTL->listend)->next = NewCandidate(n, &GarbList, TRUE); \
if ((tv->compstage < 2) && (SpineTL->listcounter <= (NAUTY_INFINITY-2))) SpineTL->listcounter++; \
SpineTL->listend = (SpineTL->listend)->next; \
CopyCand(SpineTL->listend, NextCand, n, NULL, NULL); \
} \
else { \
SpineTL->liststart = NewCandidate(n, &GarbList, TRUE); \
if (tv->compstage < 2) SpineTL->listcounter = 1; \
SpineTL->listend = SpineTL->liststart; \
CopyCand(SpineTL->liststart, NextCand, n, NULL, NULL); \
} }

#define ORBITSIZES { memset(OrbSize, 0, n*sizeof(int)); \
for (i=0; i<n; i++) { \
OrbSize[tv->orbits[i]]++; \
} }

#define CURRORBITSIZES { memset(CurrOrbSize, 0, n*sizeof(int)); \
for (i=SpineTL->tgtcell; i<SpineTL->tgtend; i++) { \
CurrOrbSize[tv->currorbit[CurrCand->lab[i]]]++; \
} }

#define EXITFROMSTAGE0REFINE { PRINT_LINE_PLUS(tv->tolevel) \
if (tv->options->verbosity >= 2) fprintf(outfile, "-=="); \
CurrCand->indnum--; \
RemoveFromLevel(tv->tolevel, tv->treedepth, tv->strategy, FALSE); \
tv->compstage = 1; \
TempOrbits = NULL; \
trieroot = trie_new(n, tv); \
trieref = trieroot; \
tv->nextlevel = tv->maxtreelevel = tv->fromlevel; \
ti->thereisnextlevel = TRUE; \
ti->exitfromref = TRUE; \
return 0; }

#define EXITFROMSTAGE0EXPATH2 { PRINT_LINE_PLUS(tv->tolevel) \
if (tv->options->verbosity >= 2) fprintf(outfile, "=-="); \
tv->compstage = 1; \
TempOrbits = NULL; \
trieroot = trie_new(n, tv); \
trieref = trieroot; \
tv->nextlevel = tv->maxtreelevel = tv->tolevel; \
ti->thereisnextlevel = TRUE; \
ti->exitfromref = FALSE; \
return 0; }

#define EXITFROMSTAGE0EXPATH1 { PRINT_RETURN PRINT_LINE_PLUS(tv->tolevel) \
if (tv->options->verbosity >= 2) fprintf(outfile, "==-"); \
if (SpineTL->liststart) { \
AuxCand = SpineTL->liststart; \
SpineTL->liststart = NewCandidate(n, &GarbList, TRUE); \
CopyCand(SpineTL->liststart, NextCand, n, NULL, NULL); \
SpineTL->liststart->next = AuxCand; \
} \
else { \
SpineTL->liststart = NewCandidate(n, &GarbList, TRUE); \
SpineTL->listend = SpineTL->liststart; \
SpineTL->liststart->next = NULL; \
CopyCand(SpineTL->liststart, NextCand, n, NULL, NULL); \
} \
tv->compstage = 1; \
TempOrbits = NULL; \
trieroot = trie_new(n, tv); \
trieref = trieroot; \
tv->nextlevel = tv->maxtreelevel = tv->tolevel; \
ti->thereisnextlevel = TRUE; \
ti->exitfromref = FALSE; \
return 0; }

#define UPDATE_LINELGTH { if (tv->options->verbosity >= 2) { \
if (tv->tolevel < 12) { \
tv->linelgth = (tv->digits+1)*tv->tolevel+16; \
} \
else { \
tv->linelgth = (tv->digits+1)*10+20; \
} \
} }

#define PRINT_LINE { if ((tv->options->verbosity >= 1) && (tv->strategy == 0)) { \
if (tv->options->verbosity >= 2) { LINE(tv->linelgth, "-"); \
NEXTLINE} \
} \
}

#define PRINT_LINE_PLUS(Lev) { if ((tv->options->verbosity >= 1) && (tv->strategy == 0)) { \
if (tv->options->verbosity >= 2) { LINE(16+tv->digits+1, "-");} \
fprintf(outfile, " Level %d:  %d cell%s; %d singleton%s; target cell: %d; %d orbit%s; %lu node%s (%lu kept); %d update%s;", \
Lev, SS(Spine[Lev].part->cells, "", "s"), SS((Spine[Lev].part->cells == n) ? n : Spine[Lev].singend, "", "s"), Spine[Lev].tgtsize, SS(tv->stats->numorbits, "", "s"), \
SS(Spine[Lev].levelcounter, "", "s"), Spine[Lev].keptcounter, SS(Spine[Lev].updates, "", "s")); \
if (Lev <= tv->group_level) fprintf(outfile, " group_level: %d", tv->group_level); \
NEXTLINE \
} \
}

#define PRINT_CANDIDATE(Cand, Lev) { \
for (tmp = Cand->name, cu = 0; tmp > 0; tmp /= 10, ++cu) {} \
for (tmp = Lev, cu1 = 0; tmp > 0; tmp /= 10, ++cu1) {} \
cu = 14-cu-cu1; \
LINE(cu, "-") \
fprintf(outfile, " %d, %d) ", Lev % 10000, Cand->name % 10000000); \
if (Lev < 12) { \
PRINTCAND(Cand, Lev) \
} \
else { \
PRINTCANDBIG(Cand, Lev) \
} \
PRINTCHAR("| ")	\
fprintf(outfile, "{%x, %x} ", Cand->code, Cand->singcode); \
}

#define PRINT_CANDIDATEPLUS(PrevCand, Cand, Lev) { \
for (tmp = Cand->name, cu = 0; tmp > 0; tmp /= 10, ++cu) {} \
for (tmp = Lev, cu1 = 0; tmp > 0; tmp /= 10, ++cu1) {} \
cu = 14-cu-cu1; \
LINE(cu, "-") \
fprintf(outfile, " %d, %d, %d) ", Lev % 10000, Cand->name % 10000000, PrevCand->name); \
if (Lev < 12) { \
PRINTCAND(Cand, Lev) \
} \
else { \
PRINTCANDBIG(Cand, Lev) \
} \
PRINTCHAR("| ")	\
fprintf(outfile, "{%x, %x} ", Cand->code, Cand->singcode); \
}

#define PRINT_EXPPATHSTEP(Cand, Boo) { \
if (tv->options->verbosity >= 2) { \
if ((tv->tolevel_tl-tv->tolevel < 6) || !has_nexttcell) { \
fprintf(outfile, "%d ", tv->indiv_vtx+labelorg); \
if (Boo) { \
if (tv->options->verbosity >= 2) fprintf(outfile, "{%d:%x} ", tv->tcellexpath, Cand->code); \
} \
else fprintf(outfile, "{interr.(%d)} ", NextPart->cells); \
if ((!has_nexttcell) && (tv->compstage == 0)) { \
fprintf(outfile, "(%d) ", tv->tolevel_tl); \
} \
} \
else { \
if (tv->tolevel_tl-tv->tolevel == 6) { \
fprintf(outfile, "... "); \
} \
} \
} \
}

#define PRINT_RETURN { if (tv->options->verbosity >= 2) { \
fprintf(outfile, "\n"); \
} }

#define PRINT_FROM_VERB(Verb,Lev) { if (tv->options->verbosity >= Verb) { \
fprintf(outfile, "FROM: "); \
if (Lev < 12) { \
PRINTCAND(CurrCand, tv->fromlevel) \
} \
else { \
PRINTCANDBIG(CurrCand, tv->fromlevel) \
} \
fprintf(outfile, " do_it: %d, indnum: %d, stnode->index: %d ", CurrCand->do_it, CurrCand->indnum, CurrCand->stnode->index); \
PRINT_RETURN \
} }

#define PRINT_NOTMIN_VERB(Verb) { if (tv->options->verbosity >= Verb) { \
fprintf(outfile, " is NOT minimal in orbits (1, %d) [%d]; ", gom_level, CurrCand->lab[Spine[gom_level+1].tgtpos]+labelorg); \
fprintf(outfile, "at lev %d, orb[%d] = %d.\n", gom_level+1, CurrCand->lab[Spine[gom_level+1].tgtpos]+labelorg, tv->currorbit[CurrCand->lab[Spine[gom_level+1].tgtpos]]+labelorg); } }

#define PRINT_SKIPPED_VERB(Verb) { if (tv->options->verbosity >= Verb) \
fprintf(outfile, " skipped (0) (orbit[%d] = %d)\n", \
NextCand->vertex+labelorg, tv->currorbit[NextCand->vertex]+labelorg); }

#define PRINT_REFINE_VERB(Verb,where) { if (tv->options->verbosity >= Verb) \
fprintf(outfile, " REFINE(%c) (orbit[%d] = %d)\n",where, NextCand->vertex+labelorg, tv->currorbit[NextCand->vertex]+labelorg); }

#define PRINT_INDIV_VERB(Verb,Lev) { if (tv->options->verbosity >= Verb) { \
if (Lev < 12) { \
PRINTCAND(CurrCand, tv->fromlevel) \
} \
else { \
PRINTCANDBIG(CurrCand, tv->fromlevel) \
} \
fprintf(outfile, "| "); \
fprintf(outfile, tv->digstring, NextCand->vertex+labelorg); \
} }

#define PRINT_INDEX(V,Verb,where) if (tv->options->verbosity >= Verb) \
fprintf(outfile,"Set index @ %d: Name %d, index %d\n",where,V->name,V->index);

#define SPECIALGENERATORS { if (tv->options->generators) addpermutation(ring, AUTPERM, n); \
tv->stats->numgenerators++; \
tv->specialgens++; \
if (tv->options->writeautoms) { \
fprintf(outfile, "Gen #%d: ", tv->stats->numgenerators); \
writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n); \
} \
if (tv->options->userautomproc) { \
(*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n); \
} }

#define UPDATEMIN(A, B) if (B < A) A = B;

#define PAIRORBJOIN(A, V) { if (A != V) { \
PrmPairs[tv->permInd].arg = A; \
PrmPairs[tv->permInd].val = V; \
tv->permInd++; \
orbjoin_sp_pair(tv->orbits, OrbList, n, A, V, &tv->stats->numorbits); \
MakeTree(A, V, sg, n, tv, FALSE); \
} }

#define SETPAIRS(A, V) { if (A != V) { \
PrmPairs[tv->permInd].arg = A; \
PrmPairs[tv->permInd].val = V; \
tv->permInd++; \
} }

#define SETPAIRSAUT(A, V) { if ((A != V) && (AUTPERM[A] != V)) { \
AUTPERM[A] = V; \
PrmPairs[tv->permInd].arg = A; \
PrmPairs[tv->permInd].val = V; \
tv->permInd++; \
} }

#define SETPAIRSAUTANDTREE(arg, val) { if (tv->build_autom) { SETPAIRSAUT(arg, val) } \
if (arg != val) orbjoin_sp_pair(tv->orbits, OrbList, n, arg, val, &tv->stats->numorbits); \
MakeTree(arg, val, sg_orig, n, tv, FALSE); }

#define PRINTF2(A, B) if (tv->options->verbosity > 3) printf(A, B)
#define PRINTF2_2(A, B, C) if (tv->options->verbosity > 3) printf(A, B, C)
#define PRINTF2_3(A, B, C, D) if (tv->options->verbosity > 3) printf(A, B, C, D)
#define PRINTF2_4(A, B, C, D, E) if (tv->options->verbosity > 3) printf(A, B, C, D, E)

/* data decls. for CPUTIME */
#ifdef  CPUDEFS
CPUDEFS
#endif

#if !MAXN
DYNALLSTAT(int, AUTPERM, AUTPERM_sz);
DYNALLSTAT(int, BreakSteps, BreakSteps_sz);
DYNALLSTAT(int, CStack, CStack_sz);
DYNALLSTAT(int, CurrOrbSize, CurrOrbSize_sz);
DYNALLSTAT(int, CurrRefCells, CurrRefCells_sz);
DYNALLSTAT(boolean, Diff, Diff_sz);
DYNALLSTAT(boolean, Factorials, Factorials_sz);
DYNALLSTAT(int, fix, fix_sz);
DYNALLSTAT(int, IDENTITY_PERM, IDENTITY_PERM_sz);
DYNALLSTAT(int, Markers, Markers_sz);
DYNALLSTAT(int, TreeMarkers, TreeMarkers_sz);
DYNALLSTAT(int, AutMarkers, AutMarkers_sz);
DYNALLSTAT(int, MarkHitVtx, MarkHitVtx_sz);
DYNALLSTAT(int, MultRefCells, MultRefCells_sz);
DYNALLSTAT(int, NghCounts, NghCounts_sz);
DYNALLSTAT(int, OrbSize, OrbSize_sz);
DYNALLSTAT(int, OrbList, OrbList_sz);
DYNALLSTAT(pair, PrmPairs, PrmPairs_sz);
DYNALLSTAT(int, TempOrbList, TempOrbList_sz);
DYNALLSTAT(int, RefCells, RefCells_sz);
DYNALLSTAT(searchtrie*, RefPath, RefPath_sz);
DYNALLSTAT(int, Singletons, Singletons_sz);
DYNALLSTAT(int, SplCls, SplCls_sz);
DYNALLSTAT(int, SplCnt, SplCnt_sz);
DYNALLSTAT(int, SplPos, SplPos_sz);
DYNALLSTAT(int, StackMarkers, StackMarkers_sz);
DYNALLSTAT(int, TheTrace, TheTrace_sz);
DYNALLSTAT(int, TheTraceCC, TheTraceCC_sz);
DYNALLSTAT(int, TheTraceSplNum, TheTraceSplNum_sz);
DYNALLSTAT(int, TheTraceSteps, TheTraceSteps_sz);
DYNALLSTAT(int, TEMPLAB, TEMPLAB_sz);
DYNALLSTAT(int, TEMPINVLAB, TEMPINVLAB_sz);
DYNALLSTAT(int, WeightsSeq, WeightsSeq_sz);
DYNALLSTAT(int, WorkArray, WorkArray_sz);
DYNALLSTAT(int, WorkArray0, WorkArray0_sz);
DYNALLSTAT(int, WorkArray1, WorkArray1_sz);
DYNALLSTAT(int, WorkArray2, WorkArray2_sz);
DYNALLSTAT(int, WorkArray3, WorkArray3_sz);
DYNALLSTAT(int, WorkArray4, WorkArray4_sz);
DYNALLSTAT(int, WorkArray5, WorkArray5_sz);
DYNALLSTAT(int, WorkArray6, WorkArray6_sz);
DYNALLSTAT(int, WorkArray7, WorkArray7_sz);
DYNALLSTAT(int, Neighbs1, Neighbs1_sz);
DYNALLSTAT(int, Neighbs2, Neighbs2_sz);
DYNALLSTAT(int, TreeStack, TreeStack_sz);
DYNALLSTAT(TracesSpine, Spine, Spine_sz);
DYNALLSTAT(trie*, TrieArray, TrieArray_sz);
DYNALLSTAT(grph_strct, TheGraph, TheGraph_sz);
DYNALLSTAT(ExpPathInfo, EPCodes, EPCodes_sz);
#else
static TLS_ATTR int CStack[MAXN];
static TLS_ATTR int AUTPERM[MAXN];
static TLS_ATTR int BreakSteps[MAXN];
static TLS_ATTR int CurrOrbSize[MAXN];
static TLS_ATTR int CurrRefCells[MAXN];
static TLS_ATTR boolean Diff[MAXN];
static TLS_ATTR boolean Factorials[MAXN];
static TLS_ATTR int fix[MAXN];
static TLS_ATTR int IDENTITY_PERM[MAXN];
static TLS_ATTR int Markers[MAXN];
static TLS_ATTR int TreeMarkers[MAXN];
static TLS_ATTR int AutMarkers[MAXN];
static TLS_ATTR int MarkHitVtx[MAXN];
static TLS_ATTR int MultRefCells[MAXN];
static TLS_ATTR int NghCounts[MAXN];
static TLS_ATTR int OrbSize[MAXN];
static TLS_ATTR int OrbList[MAXN];
static TLS_ATTR pair PrmPairs[MAXN];
static TLS_ATTR int TempOrbList[MAXN];
static TLS_ATTR int RefCells[MAXN];
static TLS_ATTR searchtrie* RefPath[MAXN];
static TLS_ATTR int Singletons[MAXN];
static TLS_ATTR int SplCls[MAXN];
static TLS_ATTR int SplCnt[MAXN];
static TLS_ATTR int SplPos[MAXN];
static TLS_ATTR int StackMarkers[MAXN];
static TLS_ATTR int TheTrace[MAXN];
static TLS_ATTR int TheTraceCC[MAXN];
static TLS_ATTR int TheTraceSplNum[MAXN];
static TLS_ATTR int TheTraceSteps[MAXN];
static TLS_ATTR int TEMPLAB[MAXN];
static TLS_ATTR int TEMPINVLAB[MAXN];
static TLS_ATTR int WeightsSeq[MAXN];
static TLS_ATTR int WorkArray[MAXN];
static TLS_ATTR int WorkArray0[MAXN];
static TLS_ATTR int WorkArray1[MAXN];
static TLS_ATTR int WorkArray2[MAXN];
static TLS_ATTR int WorkArray3[MAXN];
static TLS_ATTR int WorkArray4[MAXN];
static TLS_ATTR int WorkArray5[MAXN];
static TLS_ATTR int WorkArray6[MAXN];
static TLS_ATTR int WorkArray7[MAXN];
static TLS_ATTR int TreeStack[MAXN];
static TLS_ATTR TracesSpine Spine[MAXN];
static TLS_ATTR trie* TrieArray[MAXN];
static TLS_ATTR grph_strct TheGraph[MAXN];
static TLS_ATTR int Neighbs1[MAXN];
static TLS_ATTR int Neighbs2[MAXN];
static TLS_ATTR ExpPathInfo EPCodes[MAXN];
#endif

#define PERMSTACK WorkArray1
#define CYCLES WorkArray1
#define HitCls WorkArray2
#define CYCOLR WorkArray2
#define HitVtx WorkArray3
#define CYLGTH WorkArray3
#define CYMULT WorkArray4
#define HitCount WorkArray5
#define ElmHitCll WorkArray5
#define CYCPOS WorkArray5
#define CYCHIT TempOrbList
#define LGHATTR RefCells
#define CYCREP MultRefCells
#define TempOrbSize TEMPLAB
#define AutomCount TEMPINVLAB
#define CanonIndices MarkHitVtx
#define NSFCells NghCounts
#define TreeNodes AutMarkers
#define CellMarkers1 WorkArray6
#define CellMarkers2 WorkArray7
#define SingNonSing Singletons

static TLS_ATTR FILE *outfile;

/* Brendan's SCHREIER */
static TLS_ATTR schreier  *gpB;				/* This will point to the Schreier structure */
static TLS_ATTR permnode  *gensB;			/* This will point to the stored generators */

static TLS_ATTR Candidate *GarbList, *SpOrd, *SpCyc, *SpSwp;
static TLS_ATTR Partition *SpPart1, *SpPart2;
static TLS_ATTR TracesSpine *SpineTL, *SpineFL, *SpineTL_tl;
static TLS_ATTR trie *trieroot, *trieref;
static TLS_ATTR int *TempOrbits = NULL;
static TLS_ATTR sparsegraph redgraph;


void
Traces(sparsegraph *g_arg, int *lab, int *ptn,
       int *orbits_arg, TracesOptions *options_arg, TracesStats *stats_arg,
       sparsegraph *canong_arg) {
    int i, j;
    int tmp;
    int deg, vtx1, vtx2, *ngh1, *ngh2, *wgh1, *wgh2, ord;
    size_t j1;
    
    trielist *STStart, *STAux;
    searchtrie *TrieNode;
    int retval;
    Partition *CurrPart, *NextPart;
    Candidate *CurrCand, *NextCand, *BestCand, *AuxCand;
    
    const int n = g_arg->nv;
    const int m = SETWORDSNEEDED(n);
    
    if (g_arg->nv > (NAUTY_INFINITY-2))
    {
        fprintf(ERRFILE, "Traces: need n <= %d, but n=%d\n\n",
                NAUTY_INFINITY-2, g_arg->nv);
        return;
    }
    
    Allocate_Traces_Structures(n);
    
    struct TracesVars *tv = malloc(sizeof(struct TracesVars));
    if (tv == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
    struct TracesInfo *ti = malloc(sizeof(struct TracesInfo));
    if (ti == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
    
    trieroot = NULL;
    NextCand = GarbList = NULL;
    DYNFREE(g_arg->w,g_arg->wlen);   /* to be removed in presence of weightd edges */
    
    Initialize_Traces_Variables(tv, options_arg, stats_arg, orbits_arg, g_arg, canong_arg, n);
    
    outfile = (tv->options->outfile == NULL ? stdout : tv->options->outfile);
    
    SpOrd = SpCyc = SpSwp = NULL;
    SpPart1 = SpPart2 = NULL;
    
    if (tv->options->verbosity >= 2) {
        for (i = n, tv->digits = 0; i > 0; i /= 10, ++tv->digits) {}
        sprintf(tv->digstring, "%s%dd ", "%", tv->digits);
    }
    
    /* Initialize group and statistics */
    Initialize_Traces_Statistics(stats_arg,n);
    
    if (tv->options->verbosity >= 2) {
        Initialize_Traces_Time_Variables(tv);
    }
    
    /* Initialize lab and ptn when in the unit partition case */
    if (tv->options->defaultptn) {
        for (i = 0; i < n; i++) {
            IDENTITY_PERM[i] = i;
            ptn[i] = NAUTY_INFINITY;
        }
        ptn[n-1] = 0;
        memcpy(lab, IDENTITY_PERM, n*sizeof(int));
    } else {
        for (i = 0; i < n; i++) {
            IDENTITY_PERM[i] = i;
        }
    }
    
    memcpy(orbits_arg, IDENTITY_PERM, n*sizeof(int));
    
    if (tv->options->generators) {
        tv->stats->numorbits = given_gens(g_arg, *tv->options->generators,
                                          orbits_arg, tv->options->digraph);
        newgroup(&gpB, NULL, n);
        gensB = *tv->options->generators;
        expandschreier(gpB, &gensB, n);
        ti->thegrouphaschanged = TRUE;
        ti->identitygroup = FALSE;
        memcpy(OrbList, gpB->orbits, n*sizeof(int));
    }
    else {
        newgroup(&gpB, &gensB, n);
        memcpy(OrbList, IDENTITY_PERM, n*sizeof(int));
        tv->stats->numorbits = n;
        ti->thegrouphaschanged = FALSE;
        ti->identitygroup = TRUE;
    }
    
    copy_sg_structure(&redgraph, g_arg);
    
    tv->graph = &redgraph;
    if (g_arg->w) memcpy(tv->graph->w, g_arg->w, tv->graph->wlen*sizeof(int));
    memcpy(tv->graph->e, g_arg->e, tv->graph->elen*sizeof(int));
    
    for (i=0; i<n; i++) {
        EPCodes[i].info = 0;
        TheGraph[i].d = g_arg->d[i];
        if (TheGraph[i].d > tv->maxdeg) {
            tv->maxdeg = TheGraph[i].d;
        }
        if (TheGraph[i].d < tv->mindeg) {
            tv->mindeg = TheGraph[i].d;
        }
        TheGraph[i].e = tv->graph->e + g_arg->v[i];
        if (g_arg->w)
            TheGraph[i].w = tv->graph->w + g_arg->v[i];
        else
            TheGraph[i].w = NULL;
        TheGraph[i].one = FALSE;
    }
    
    ord = 0;
    
    /*----------- WEIGHTS --------------*/
    if (tv->options->weighted) {
        WeightCodes(n);
        ord = trie_classify(n,tv);
    }
    /*----------------------------------*/
    
    if ((tv->maxdeg == tv->mindeg) && (ord == 0)) ti->regular = TRUE; else ti->regular = FALSE;
    
    tv->currorbit = gpB->orbits;
    
    memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
    tv->permInd = 0;
    
    memset(fix, 0, n*sizeof(int));
    memset(TheTraceCC, 0, n*sizeof(int));
    memset(Factorials, 0, n*sizeof(int));
    /* ran_init(1234);  any long int as an argument */
    
    /* The graph is sparse? */
    if (g_arg->nde < n || g_arg->nde / n < n / (g_arg->nde / n)) {
        ti->thegraphisparse = TRUE;
    }
    else {
        ti->thegraphisparse = FALSE;
    }
    
    tv->preprocessed = 0;
    ti->deg_one = FALSE;
    ti->first_matching = FALSE;
    retval = 0;
    
    /* Initialize candidate, partition, cells, orbits */
    Spine[0].part = NewPartition(n);
    CurrPart = Spine[0].part;
    memset(CurrPart->inv, 0, n*sizeof(int));
    
    NextPart = NewPartition(n);
    CurrCand = NewCandidate(n, &GarbList, TRUE);
    
    CurrCand->singcode = 0;
    TempOrbits = NULL;
    STStart = NULL;
    
    if (ti->regular) {
        if (tv->options->defaultptn) {
            memcpy(CurrCand->lab, IDENTITY_PERM, n*sizeof(int));
            memcpy(CurrCand->invlab, IDENTITY_PERM, n*sizeof(int));
            CurrPart->cells = 1;
            CurrPart->cls[0] = n;
            TheTrace[0] = 0;
        }
        else {
            memcpy(CurrCand->lab, lab, n*sizeof(int));
            CurrPart->cells = 0;
            j = 0;
            for (i = 0; i < n; i++) {
                if (j) CurrPart->inv[i] = j;
                CurrCand->invlab[CurrCand->lab[i]] = i;
                if (!ptn[i]) {
                    CurrPart->cls[j] = i-j+1;
                    if (CurrPart->cls[j] == 1) {
                        CurrCand->singcode = MASHCOMM(CurrCand->singcode, CurrCand->lab[j]);
                    }
                    TheTrace[CurrPart->cells++] = j;
                    j = i+1;
                }
            }
        }
    } else {
        if (tv->options->weighted) {
            CurrPart->cells = traces_vertexclass_refine(n, lab, ptn,
                                                        CurrCand, CurrPart, WeightsSeq);
        }
        else {
            CurrPart->cells = traces_vertexclass_refine (n, lab, ptn,
                                                         CurrCand, CurrPart, g_arg->d);
        }
    }
    
    memset(NghCounts,0,n*sizeof(int));
    if (tv->options->verbosity == 7) PrintPartition(CurrCand->lab,CurrPart->cls,n,labelorg,1323);
    
    /* Check for deg 1 vertices */
    ti->deg_one = Check_degree_one(g_arg, CurrCand, CurrPart, n);
    
#if !MAXN
    DYNALLOC1(int, Neighbs1, Neighbs1_sz, tv->maxdeg, "Traces");
    DYNALLOC1(int, Neighbs2, Neighbs2_sz, tv->maxdeg, "Traces");
#endif
    
    if (ti->deg_one) {
        tv->preprocessed = Preprocess(g_arg, &gensB, CurrCand, n, CurrPart, tv);
    }
    
    if (tv->preprocessed) {
        memset(Diff,0,n*sizeof(boolean));
        for (i=0; i<n; i++) {
            if ((TheGraph[i].d > 1) && (tv->input_graph->d[i] != TheGraph[i].d))
                Diff[i] = TRUE;
        }
    }
    
    /* Initialization of Spine structure */
    SpineFL = Spine;
    SpineFL->tgtcell = SpineFL->tgtpos = 0;
    SpineFL->tgtend = n;
    SpineFL->tgtfrom = -1;
    SpineFL->trcstart = 0;
    SpineFL->trcend = CurrPart->active = CurrPart->cells;
    SpineFL->ccstart = SpineFL->ccend = 0;
    SpineFL->stpstart = SpineFL->stpend = 0;
    SpineFL->singstart = SpineFL->singend = 0;
    SpineFL->thetracexists = FALSE;
    SpineFL->liststart = SpineFL->listend = CurrCand;
    SpineFL->levelcounter = 1;
    SpineFL->keptcounter = 1;
    SpineFL->updates = 1;
    
    /* Further initializations */
    tv->maxtreelevel = 0;
    tv->tolevel = 0;
    tv->tcell = 0;
    UPDATE_LINELGTH
    
    /* First refinement */
    if (tv->preprocessed < 2)
        traces_refine(CurrCand, n, CurrPart, tv, ti, 0, FALSE);
    
    CurrCand->name = 0;
    
    if (CurrPart->cells == n) {
        tv->stats->canupdates++;
        
        /* CANONICAL FORM ? */
        if (tv->options->getcanon) {
            memcpy(lab, CurrCand->lab, n*sizeof(int));
            if (canong_arg) memcpy(CurrPart->inv, CurrCand->invlab, n*sizeof(int));
            if (tv->options->usercanonproc != NULL)
            {
                (*tv->options->usercanonproc)((graph*)g_arg, lab, (graph*)canong_arg, tv->stats->canupdates, CurrCand->code, m, n);
            }
        }
    }
    else {
        STStart = searchtrie_new(n, tv);
        CurrCand->stnode = tv->strielist->triearray;
        NextCand = NewCandidate(n, &GarbList, TRUE);
        SpineFL->listcounter = 1;
        tv->tcellevel = 0;
        ti->exitfromref = FALSE;
        Spine[1].levelcounter = 0;
        Spine[1].updates = 0;
        Spine[1].tgtfrom = 0;
        
        memset(WorkArray, 0, n*sizeof(int));
        
        do {
            tv->fromlevel = tv->tolevel;
            SpineFL = Spine+tv->fromlevel;
            
            if (CurrCand) {
                switch (tv->compstage) {
                    case 0: retval = CompStage0(CurrPart, NextPart, CurrCand, NextCand, m, n, tv, ti);
                        break;
                    case 1:
                        if (!TempOrbits) {
                            memcpy(WorkArray1, tv->currorbit, n*sizeof(int));
                            TempOrbits = WorkArray1;
                        }
                        memset(TempOrbSize, 0, n*sizeof(int));
                        for (i=SpineTL->tgtcell; i<SpineTL->tgtend; i++) {
                            TempOrbSize[TempOrbits[CurrCand->lab[i]]]++;
                        }
                        retval = CompStage1(CurrPart, NextPart, CurrCand, NextCand, m, n, tv, ti);
                        break;
                    case 2:
                        retval = CompStage2(CurrPart, NextPart, CurrCand, NextCand, m, n, tv, ti);
                        break;
                    default:
                        break;
                }
                if (retval == NAUTY_ABORTED)
                    tv->stats->errstatus = NAUABORTED;
                else if (retval == NAUTY_KILLED) {
                    tv->stats->errstatus = NAUKILLED;
                    return;
                }
            }
            
            /* NEXT CANDIDATE */
            if (ti->thereisnextlevel) {
                if (tv->nextlevel != tv->fromlevel) {
                    UPDATE_LINELGTH
                }
                tv->tolevel = tv->nextlevel;
                CurrCand = Spine[tv->nextlevel].liststart;
                CurrPart = Spine[tv->nextlevel].part;
            }
        }
        while (ti->thereisnextlevel);
        
        if (!retval) {
            if (tv->compstage) {
                memset(CurrOrbSize, 0, n*sizeof(int));
                for (i=0; i<n; i++) {
                    CurrOrbSize[TempOrbits[i]]++;
                }
            }
            
            if (!tv->options->getcanon) {
                if (tv->compstage) {
                    tv->maxtreelevel++;
                    TrieNode = Spine[tv->maxtreelevel].liststart->stnode;
                    if (TrieNode->father) {
                        TrieNode = TrieNode->father;
                    }
                    if (!ti->exitfromref) {
                        TrieNode->index = 0;
                        for (i=1; i<AutomCount[0]; i++) {
                            TrieNode->index += CurrOrbSize[AutomCount[i]];
                            PRINT_INDEX(TrieNode,4,30)
                        }
                    }
                }
            }
            
            if (tv->maxtreelevel) PRINT_LINE_PLUS(tv->maxtreelevel);
            
            AuxCand = Spine[tv->maxtreelevel].liststart;
            while (!AuxCand->do_it) {
                AuxCand = AuxCand->next;
            }
            
            /* CANONICAL FORM ? */
            if (tv->options->getcanon) {
                BestCand = AuxCand;
                AuxCand = AuxCand->next;
                while (AuxCand) {
                    if (AuxCand->do_it) {
                        if (comparelab_tr(g_arg, BestCand->lab, BestCand->invlab, AuxCand->lab, AuxCand->invlab, Spine[tv->maxtreelevel].part->cls, Spine[tv->maxtreelevel].part->inv) == 1) {
                            BestCand = AuxCand;
                            tv->stats->canupdates++;
                            if (tv->options->usercanonproc != NULL)
                            {
                                (*tv->options->usercanonproc)((graph*)g_arg, lab, (graph*)canong_arg, tv->stats->canupdates, CurrCand->code, m, n);
                            }
                        }
                    }
                    AuxCand = AuxCand->next;
                }
                
                grouporderplus(g_arg, BestCand, Spine[tv->maxtreelevel].part, &gensB, &(tv->stats->grpsize1), &(tv->stats->grpsize2), n, tv, ti);
                
                if (tv->options->verbosity >= 2) {
                    LINE(32, "—")
                    NEXTLINE
                    fprintf(outfile, "Canonical:");
                    PRINTCAND(BestCand, tv->maxtreelevel)
                    PRINT_RETURN
                }
                memcpy(lab, BestCand->lab, n*sizeof(int));
                if (canong_arg) memcpy(CurrPart->inv, BestCand->invlab, n*sizeof(int));
            }
            else {
                grouporderplus(g_arg, AuxCand, Spine[tv->maxtreelevel].part, &gensB, &(tv->stats->grpsize1), &(tv->stats->grpsize2), n, tv, ti);
            }
            
            if (tv->options->verbosity >= 2) {
                if (tv->linelgth < 40) {
                    tv->linelgth = 40;
                }
                LINE(32, "—")
                NEXTLINE
            }
            
        }
    }
    tv->stats->treedepth = tv->treedepth;
    if (Spine[tv->treedepth].part->code == -1) {
        tv->stats->treedepth--;
    }
    
    if (tv->options->verbosity >= 2) {
        fprintf(outfile, "group time: %.2f, %.2f, %.2f; total: %.2f; exp_paths time: %.2f; aut_check time: %.2f\n%lu refinement%s interrupted by trace comparison (%s); special cells: %d\n------", tv->schreier1, tv->schreier2, tv->schreier3, tv->schreier1+tv->schreier2+tv->schreier3,
                tv->expaths, tv->autchk, SS(tv->stats->interrupted, "", "s"), (tv->options->strategy == 0 ? "breadth-first" : "depth-first"), tv->specialgens);
        PRINT_RETURN
    }
    if (tv->options->verbosity >= 3) fprintf(outfile, "CPYCAND(0): %d, ID<-TMPORB(1): %d, LAB(2): %d, PART(3): %d, TEMP(4)->: %d, TEMP(5)<-: %d, CHKFORAUT(6): %d, ISAUT(7): %d, ContaTC: %d\n", tv->conta0, tv->conta1, tv->conta2, tv->conta3, tv->conta4, tv->conta5, tv->conta6, tv->conta7, tv->contatc);
    
    if (tv->options->getcanon && canong_arg) {
        canong_arg->nv  = g_arg->nv;
        canong_arg->nde  = g_arg->nde;
        SG_ALLOC(*canong_arg, g_arg->nv, g_arg->nde, "traces canong");
        updatecan_tr(g_arg, canong_arg, lab, CurrPart->inv, 0);
    }
    
    if (tv->options->generators) {
        deleteunmarked(&gensB);
        *tv->options->generators = gensB;
        freeschreier(&gpB, NULL);
    }
    else {
        freeschreier(&gpB, &gensB);
        schreier_freedyn();
    }
    
    while (STStart) {
        STAux = STStart;
        free(STAux->triearray);
        STStart = STStart->next;
        free(STAux);
    }
    
    tv->canlist = 0;
    for (i=0; i<=tv->treedepth; i++) {
        if (Spine[i].liststart) {
            tv->canlist += FreeList(Spine[i].liststart, TRUE);
            Spine[i].liststart = Spine[i].listend = NULL;
        }
    }
    
    if (GarbList) {
        tv->stats->peaknodes = FreeList(GarbList, FALSE);
    }
    
    FREECAND(NextCand)
    FREECAND(SpOrd)
    FREECAND(SpCyc)
    FREECAND(SpSwp)
    FREEPART(NextPart)
    FREEPART(SpPart1)
    FREEPART(SpPart2)
    
    if (!tv->options->getcanon && trieroot) {
        for (i=0; i<=tv->triepos; i++) {
            free(TrieArray[i]);
        }
    }
    
    for (i=0; i <= tv->treedepth; i++) {
        FREEPART(Spine[i].part)
    }
    
    CurrCand = GarbList = NULL;
    tv->stats->peaknodes += tv->canlist;
    
    if (tv->graph != g_arg) {
        SG_FREE(redgraph);
    }
    free(tv);
    free(ti);
    traces_freedyn();
    
    return;
}

int traces_vertexclass_refine (int n, int *lab, int *ptn, Candidate *Cand, Partition *Part, int *RefArray) {
    
    int i, j, k, aux, cells, end;
    
    memcpy(Cand->lab, lab, n*sizeof(int));
    
    cells = 0;
    j = 0;
    
    for (i = 0; i < n; i++) {
        WorkArray1[i] = RefArray[Cand->lab[i]];
        if (!ptn[i]) {
            TheTrace[cells++] = j;
            
            sort2ints(WorkArray1+j,Cand->lab+j,i-j+1);
            
            aux = WorkArray1[j];
            Part->cls[j] = 1;
            Part->inv[j] = j;
            Cand->invlab[Cand->lab[j]] = j;
            
            if (i == j) {
                Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[j]);
                j++;
            } else {
                end = i+1;
                for (k=j+1; k<end; k++) {
                    if (WorkArray1[k] == aux) {
                        (Part->cls[j])++;
                        Part->inv[k] = j;
                        Cand->invlab[Cand->lab[k]] = k;
                    } else {
                        if (Part->cls[j] == 1) {
                            Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[j]);
                        }
                        
                        TheTrace[cells++] = k;
                        j = k;
                        aux = WorkArray1[j];
                        Part->cls[j] = 1;
                        Part->inv[j] = j;
                        Cand->invlab[Cand->lab[j]] = j;
                    }
                }
                j = end;
            }
        }
    }
    return cells;
}

int traces_refine(Candidate *Cand,
                  int n,
                  Partition *Part,
                  struct TracesVars* tv,
                  struct TracesInfo *ti,
                  int num_indv,
                  boolean make_code) {
    
    int i, j, k, jk, sc, ind0, ind1, ind2, ind3, ind4, tlp1, labi;
    int value, iend, newcell;
    int HitClsInd, SplInd, SplCntInd, CStackInd, TraceInd, TraceCCInd, TraceStepsInd, SingInd;
    int j1int, iend1int;
    unsigned int longcode;
    int newtrace = FALSE;
    int Sparse = TRUE;
    int *lab, *cls, *InvLab, *TracePos, *SplitCell, *LabCell, *TraceEnd, Traceccend, *Tracestpend;
    int BigCell, BigCellPos, BigCellSize;
    boolean TraceCell = FALSE;
    int *nghb;
    int conta;
    const int variation = 0;
    int currentweight, weightstart, weightend, currentcell, currentsize;
    
    HitClsInd = 0;
    if (tv->stackmark > (NAUTY_INFINITY-2)) {
        memset(StackMarkers, 0, n*sizeof(int));
        tv->stackmark = 0;
    }
    tv->stackmark++;
    tv->augmented_cells = Part->cells;
    
    SpineTL = Spine+tv->tolevel;
    TraceEnd = &(SpineTL->trcend);
    Traceccend = SpineTL->ccend;
    Tracestpend = &(SpineTL->stpend);
    TraceCCInd = SpineTL->ccstart;
    TraceStepsInd = SpineTL->stpstart;
    
    SingInd = SpineTL->singstart + num_indv;
    
    lab = Cand->lab;
    InvLab = Cand->invlab;
    cls = Part->cls;
    
    UPDATEMIN(Part->active, n-1);
    memcpy(CStack+1, TheTrace+SpineTL->trcstart, (Part->active)*sizeof(int));
    CStackInd = Part->active;
    for (i = 1; i <= CStackInd; i++) {
        StackMarkers[CStack[i]] = tv->stackmark;
    }
    
    longcode = Part->cells;
    TraceInd = SpineTL->trcstart+Part->active;
    if (!SpineTL->thetracexists) {
        newtrace = TRUE;
    }
    conta=0;
    
    while (CStackInd > 0) {
        
        weightend = 0;
        
        if (tv->mark > (NAUTY_INFINITY-2)) {
            memset(Markers, 0, n*sizeof(int));
            memset(MarkHitVtx, 0, n*sizeof(int));
            tv->mark = 0;
        }
        tv->mark++;
        
        if (Part->cells == n) break;
        
        k = Select_from_CStack(cls, CStackInd);
        
        currentcell = CStack[k];
        currentsize = currentcell+cls[currentcell];
        CStack[k] = CStack[CStackInd--];
        StackMarkers[currentcell] = 0;
        
        labi = lab[currentcell];
        iend1int = TheGraph[labi].d;
        nghb =  TheGraph[labi].e;
        
        do {
            ind0 = currentcell;
            ind2 = currentsize;
            weightstart = weightend;
            
            if (tv->options->weighted) {
                currentweight = (TheGraph[labi].w)[weightstart];
                while ((iend1int > weightend) && ((TheGraph[labi].w)[weightend] == currentweight)) {
                    weightend++;
                }
            } else {
                weightend = TheGraph[labi].d;
            }
            
            if (!newtrace) {
                TraceCell = ((ind0 == TheTraceCC[TraceCCInd]) && (TraceCCInd < Traceccend));
            }
            
            /* Analysis of occurrences of neighbors of the current cell */
            /* The list of cells with neighbors in the current cell is  built */
            if (cls[ind0] == 1) {			/* SINGLETON CURRENT CELL CASE */
                
                /* NEIGHCOUNT_SING_MULT */
                HitClsInd = 0;
                for (j1int = weightstart; j1int < weightend; ++j1int) {
                    k = nghb[j1int];
                    value = Part->inv[InvLab[k]];
                    if (cls[value] > 1) {
                        if (Markers[value] != tv->mark) {
                            HitCls[HitClsInd++] = value;
                            Markers[value] = tv->mark;
                            ElmHitCll[value] = value;
                        }
                        HitVtx[ElmHitCll[value]++] = k;
                    } else {
                        switch (variation) {
                            case 1:
                                longcode = MASHCOMM(longcode, value);
                                break;
                            default:
                                break;
                        }
                    }
                }
                /* end NEIGHCOUNT_SING_MULT */
                
                tv->mark++;
                FIND_SPLIT_CELLS;
                
                /* SINGLETON CC CASE */
                if (SplInd) {
                    if (newtrace) {
                        TheTraceCC[TraceCCInd] = ind0;
                    }
                    if (!TraceCell) {
                        TheTraceCC[TraceCCInd] = ind0;
                        newtrace = TRUE;
                    }
                    
                    TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
                    
                    /* Sorting the cells to be split */
                    sort_Split_Array(SplCls, SplInd);
                    
                    for (j = 0; j < SplInd; j++) {
                        ind1 = SplCls[j];
                        i = ind1+cls[ind1]-ElmHitCll[ind1];
                        TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
                    }
                    
                    for (j = 0; j < SplInd; j++) {
                        /* REARRANGE_CELLS */
                        ind1 = SplCls[j];
                        cls[ind1] = cls[ind1]-ElmHitCll[ind1];
                        newcell = ind1+cls[ind1];
                        cls[newcell] = ElmHitCll[ind1];
                        Part->cells++;
                        if (StackMarkers[ind1] != tv->stackmark) {
                            if (cls[newcell] < cls[ind1]) {
                                CStack[++CStackInd] = newcell;
                                StackMarkers[newcell] = tv->stackmark;
                            }
                            else {
                                CStack[++CStackInd] = ind1;
                                StackMarkers[ind1] = tv->stackmark;
                            }
                        }
                        else {
                            CStack[++CStackInd] = newcell;
                            StackMarkers[newcell] = tv->stackmark;
                        }
                        SplitCell = HitVtx+ind1;
                        ind3 = cls[newcell];
                        LabCell = lab+newcell;
                        for (jk = 0; jk < ind3; jk++) {
                            k = SplitCell[jk];
                            i = LabCell[jk];
                            Part->inv[newcell+jk] = newcell;
                            lab[InvLab[k]] = i;
                            InvLab[i] = InvLab[k];
                            LabCell[jk] = k;
                            InvLab[k] = newcell+jk;
                        }
                        /* END REARRANGE_CELLS */
                        
                        if (cls[ind1] == 1) {
                            Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[ind1]);
                            if (newtrace) Singletons[SingInd] = ind1;
                            SingInd++;
                        }
                        if (cls[newcell] == 1) {
                            Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[newcell]);
                            if (newtrace) Singletons[SingInd] = newcell;
                            SingInd++;
                        }
                    }
                }
                else {
                    if ((!newtrace) && TraceCell) {
                        if (!tv->options->weighted) return 0;
                    }
                }
            }
            else {
                if (ti->thegraphisparse) {
                    
                    /* NEIGHCOUNT_SPARSE_MULT */
                    if (cls[ind0] != n) {
                        HitClsInd = 0;
                        for (i = ind0; i < ind2; i++) {
                            labi = lab[i];
                            nghb = TheGraph[labi].e;
                            for (j = weightstart; j < weightend; j++) {
                                k = nghb[j];
                                if (MarkHitVtx[k] == tv->mark) {
                                    NghCounts[k]++;
                                }
                                else {
                                    value = Part->inv[InvLab[k]];
                                    if (cls[value] > 1) {
                                        MarkHitVtx[k] = tv->mark;
                                        NghCounts[k] = 1;
                                        if (Markers[value] != tv->mark) {
                                            HitCls[HitClsInd++] = value;
                                            Markers[value] = tv->mark;
                                            HitVtx[value] = k;
                                            ElmHitCll[value] = 1;
                                        }
                                        else {
                                            HitVtx[value+ElmHitCll[value]++] = k;
                                        }
                                    }
                                    else {
                                        switch (variation) {
                                            case 1:
                                                longcode = MASHCOMM(longcode, value);
                                                break;
                                            default:
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /* End NEIGHCOUNT_SPARSE_MULT */
                    
                    tv->mark++;
                    
                    SplInd = 0;
                    SplCls[0] = n;
                    for (j = 0; j < HitClsInd; j++) {
                        ind1 = HitCls[j];
                        if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
                            SplCls[SplInd++] = ind1;
                        }
                        else {
                            ind2 = ind1+cls[ind1];
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    break;
                                }
                            }
                        }
                    }
                    
                    /* SPARSE CASE */
                    if (SplInd) {
                        if (newtrace) {
                            TheTraceCC[TraceCCInd] = ind0;
                        }
                        if (!TraceCell) {
                            TheTraceCC[TraceCCInd] = ind0;
                            newtrace = TRUE;
                        }
                        TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+n, &Traceccend)
                        
                        /* Sorting the cells to be split */
                        sort_Split_Array(SplCls, SplInd);
                        
                        for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
                            ind0 = SplCls[sc];
                            ind1 = ind0 + cls[ind0];
                            SplCntInd = 0;
                            if (ElmHitCll[ind0] < cls[ind0]) {
                                SplCnt[SplCntInd++] = 0;
                                SplPos[0] = cls[ind0] - ElmHitCll[ind0];
                            }
                            
                            /* According to the numbers of neighbors of C into the current cell */
                            /* compute how many vertices in C will be placed into the same new cell */
                            iend = ind0 + ElmHitCll[ind0];
                            for (i = ind0; i < iend; i++) {
                                value = NghCounts[HitVtx[i]];
                                if (Markers[value] != tv->mark) {
                                    Markers[value] = tv->mark;
                                    SplCnt[SplCntInd++] = value;
                                    SplPos[value] = 1;
                                }
                                else {
                                    SplPos[value]++;
                                }
                            }
                            tv->mark++;
                            
                            if (SplCntInd) {
                                TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+n, Tracestpend)
                            }
                            
                            /* Sort the values deriving from the previous step */
                            sort_Split_Array(SplCnt,SplCntInd);
                            Part->cells += SplCntInd-1;
                            
                            /* Split the cell C and update the information for sizes of new cells */
                            /* Put the new cells into the stack */
                            i = ind0;
                            if (StackMarkers[i] != tv->stackmark) {
                                BigCellSize = 0;
                            }
                            for (k = 0; k < SplCntInd; k++) {
                                value = SplPos[SplCnt[k]];
                                cls[i] = value;
                                if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                    BigCell = i;
                                    BigCellPos = CStackInd;
                                    BigCellSize = cls[i];
                                }
                                SplPos[SplCnt[k]] = i;
                                i += value;
                                if (i < ind1) {
                                    CStack[++CStackInd] = i;
                                    StackMarkers[i] = tv->stackmark;
                                    TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
                                }
                            }
                            
                            if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                CStack[BigCellPos] = ind0;
                                StackMarkers[BigCell] = 0;
                                StackMarkers[ind0] = tv->stackmark;
                            }
                            /* Permute elements of the cell C */
                            iend = ind0 + ElmHitCll[ind0];
                            for (i = ind0; i < iend; i++) {
                                value = HitVtx[i];
                                j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
                                k = InvLab[value];				/* where HitVtx[i] is in lab */
                                lab[k] = lab[j];
                                lab[j] = value;
                                InvLab[value] = j;
                                InvLab[lab[k]] = k;
                                NghCounts[value] = 0;
                            }
                            
                            /* Reconstruct the cell C and update the inverse partition */
                            newcell = ind1 - ElmHitCll[ind0];
                            i = newcell;
                            ind2 = newcell+cls[newcell]-1;
                            do {
                                Part->inv[i] = newcell;
                                if (i == ind2) {
                                    newcell = i+1;
                                    if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                }
                            }
                            while (++i < ind1);
                            
                            for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
                                if (cls[i] == 1) {
                                    Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                    if (newtrace) Singletons[SingInd] = i;
                                    SingInd++;
                                }
                            }
                            
                        }
                    }
                    else {
                        if ((!newtrace) && TraceCell) {
                            return 0;
                        }
                    }
                    
                }
                else {
                    if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
                        Sparse = FALSE;
                    }
                    else {
                        Sparse = TRUE;
                    }
                    if (Sparse) {
                        /* Counting occurrences of neighbors of the current cell */
                        /* The list of cells with neighbors in the current cell is also built */
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                            HitCls[0] = 0;
                            HitClsInd = 1;
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_SPARSE_MULT */
                            HitClsInd = 0;
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb =  TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                    value = Part->inv[InvLab[k]];
                                    if (Markers[value] != tv->mark) {
                                        if (cls[value] > 1) HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                    }
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_SPARSE_MULT */
                            
                        }
                        
                        tv->mark++;
                        
                        
                        SplInd = 0;
                        for (j = 0; j < HitClsInd; j++) {
                            ind1 = HitCls[j];
                            ind2 = ind1+cls[ind1];
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    break;
                                }
                            }
                        }
                        
                        /* DENSE-SPARSE CASE */
                        if (SplInd) {
                            if (newtrace) {
                                TheTraceCC[TraceCCInd] = ind0;
                            }
                            if (!TraceCell) {
                                TheTraceCC[TraceCCInd] = ind0;
                                newtrace = TRUE;
                            }
                            TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+2*n, &Traceccend)
                            
                            /* Sorting the cells to be split */
                            sort_Split_Array(SplCls, SplInd);
                            
                            for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
                                ind0 = SplCls[j];
                                ind1 = ind0+cls[ind0];
                                SplCntInd = 0;
                                
                                /* According to the numbers of neighbors of C into the current cell */
                                /* compute how many vertices in C will be placed into the same new cell */
                                for (i = ind0; i < ind1; i++) {
                                    value = NghCounts[lab[i]];
                                    if (Markers[value] != tv->mark) {
                                        Markers[value] = tv->mark;
                                        SplCnt[SplCntInd++] = value;
                                        SplPos[value] = 1;
                                    }
                                    else {
                                        SplPos[value]++;
                                    }
                                }
                                tv->mark++;
                                
                                if (SplCntInd) {
                                    TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+2*n, Tracestpend)
                                }
                                
                                /* Sort the values deriving from the previous step */
                                sort_Split_Array(SplCnt, SplCntInd);
                                
                                Part->cells += SplCntInd-1;
                                
                                /* Split the cell C and update the information for sizes of new cells */
                                /* Put the new cells into the stack */
                                i = ind0;
                                if (StackMarkers[i] != tv->stackmark) {
                                    BigCellSize = 0;
                                }
                                for (k = 0; k < SplCntInd; k++) {
                                    value = SplPos[SplCnt[k]];
                                    cls[i] = value;
                                    
                                    if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                        BigCell = i;
                                        BigCellPos = CStackInd;
                                        BigCellSize = cls[i];
                                    }
                                    SplPos[SplCnt[k]] = i;
                                    i += value;
                                    if (i < ind1) {
                                        CStack[++CStackInd] = i;
                                        StackMarkers[i] = tv->stackmark;
                                        TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
                                    }
                                }
                                
                                if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                    CStack[BigCellPos] = ind0;
                                    StackMarkers[BigCell] = 0;
                                    StackMarkers[ind0] = tv->stackmark;
                                }
                                
                                /* Permute elements of the cell C */
                                i = ind0;
                                do {
                                    SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                                }
                                while(++i < ind1);
                                
                                /* Reconstruct the cell C and update the inverse partition */
                                newcell = ind0;
                                i = ind0;
                                ind2 = newcell+cls[newcell]-1;
                                do {
                                    lab[i] = SplCnt[i];
                                    InvLab[lab[i]] = i;
                                    
                                    Part->inv[i] = newcell;
                                    if (i == ind2) {
                                        newcell = i+1;
                                        if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                    }
                                }
                                while (++i < ind1);
                                
                                for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
                                    if (cls[i] == 1) {
                                        Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                        if (newtrace) Singletons[SingInd] = i;
                                        SingInd++;
                                    }
                                }
                                
                            }
                        }
                        else {
                            if ((!newtrace) && TraceCell) {
                                return 0;
                            }
                        }
                        
                    }
                    else {
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_DENSE_MULT */
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb = TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_DENSE_MULT */
                            
                        }
                        SplInd = 0;
                        ind4 = 0;
                        while (ind4 < n) {	/* For each cell C with size(C) > 1 */
                            ind1 = ind4+cls[ind4];
                            if (cls[ind4] > 1) {
                                
                                
                                /* Determine whether C must be split */
                                SplCntInd = 0;
                                value = NghCounts[lab[ind4]];
                                for (i = ind4+1; i < ind1; i++) {
                                    if (NghCounts[lab[i]] != value)
                                    {
                                        SplInd++;
                                        Markers[value] = tv->mark;
                                        SplCnt[SplCntInd++] = value;
                                        SplPos[value] = i-ind4;
                                        do {
                                            value = NghCounts[lab[i++]];
                                            if (Markers[value] != tv->mark) {
                                                Markers[value] = tv->mark;
                                                SplCnt[SplCntInd++] = value;
                                                SplPos[value] = 1;
                                            }
                                            else {
                                                SplPos[value]++;
                                            }
                                        }
                                        while(i != ind1);
                                        break;
                                    }
                                }
                                tv->mark++;
                                
                                if (SplInd && !TraceCell) newtrace = TRUE;
                                
                                if (SplCntInd) {
                                    TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+3*n, Tracestpend)
                                    
                                    /* Sort the values deriving from the previous step */
                                    sort_Split_Array(SplCnt, SplCntInd);
                                    
                                    Part->cells += SplCntInd-1;
                                    
                                    /* Split the cell C and update the information for sizes of new cells */
                                    /* Put the new cells into the stack */
                                    i = ind4;
                                    if (StackMarkers[i] != tv->stackmark) {
                                        BigCellSize = 0;
                                    }
                                    for (k = 0; k < SplCntInd; k++) {
                                        value = SplPos[SplCnt[k]];
                                        cls[i] = value;
                                        
                                        if ((StackMarkers[ind4] != tv->stackmark) && (value > BigCellSize)) {
                                            BigCell = i;
                                            BigCellPos = CStackInd;
                                            BigCellSize = cls[i];
                                        }
                                        SplPos[SplCnt[k]] = i;
                                        i += value;
                                        if (i < ind1) {
                                            CStack[++CStackInd] = i;
                                            StackMarkers[i] = tv->stackmark;
                                            TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
                                        }
                                    }
                                    if ((StackMarkers[ind4] != tv->stackmark) && (BigCell != ind4)) {
                                        CStack[BigCellPos] = ind4;
                                        StackMarkers[BigCell] = 0;
                                        StackMarkers[ind4] = tv->stackmark;
                                    }
                                    
                                    /* Permute elements of the cell C */
                                    i = ind4;
                                    do {
                                        SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                                    }
                                    while(++i < ind1);
                                    
                                    /* Reconstruct the cell C and update the inverse partition */
                                    newcell = ind4;
                                    i = ind4;
                                    ind2 = newcell+cls[newcell]-1;
                                    do {
                                        lab[i] = SplCnt[i];
                                        InvLab[lab[i]] = i;
                                        Part->inv[i] = newcell;
                                        if (i == ind2) {
                                            newcell = i+1;
                                            if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                        }
                                    }
                                    while (++i < ind1);
                                    
                                    for (i = ind4; i < ind1; i+=cls[i]) {
                                        if (cls[i] == 1) {
                                            Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                            if (newtrace) Singletons[SingInd] = i;
                                            SingInd++;
                                        }
                                    }
                                }
                            }
                            ind4 = ind1;
                        }
                        
                        /* DENSE-DENSE CASE */
                        if (SplInd) {
                            if (!TraceCell) {
                                newtrace = TRUE;
                            }
                            TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+3*n, &Traceccend)
                            if (newtrace) {
                                TheTraceCC[TraceCCInd-1] = ind0;
                            }
                        }
                        else {
                            if ((!newtrace) && TraceCell) {
                                return 0;
                            }
                        }
                        
                    }
                }
            }
        }
        while (weightend < iend1int);
    }  /* end while (CStackInd > 0) */
    
    tv->augmented_cells = Part->cells - tv->augmented_cells;
    
    if (make_code) {
        for (i=SpineTL->trcstart; i < TraceInd; i++) {
            ind0 = TheTrace[i];
            longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
            labi = lab[ind0];
            iend1int = TheGraph[labi].d;
            nghb = TheGraph[labi].e;
            for (j1int = 0; j1int < iend1int; ++j1int) {
                k = nghb[j1int];
                value = Part->inv[InvLab[k]];
                longcode = MASHCOMM(longcode, value);
            }
        }
    }
    Part->code = Cand->code = CLEANUP(longcode);
    tlp1 = tv->tolevel+1;
    if (newtrace) {
        if ((tlp1 < n) && (tlp1 > tv->treedepth)) {
            tv->treedepth = tlp1;
            if (tv->strategy) {
                Spine[tlp1].part = NewPartition(n);
            }
            else {
                NewPartSpine(tlp1,n);
            }
            Spine[tlp1].liststart = Spine[tlp1].listend = NULL;
            Spine[tlp1].listcounter = 0;
        }
        *TraceEnd = TraceInd;
        SpineTL->ccend = TraceCCInd;
        SpineTL->singend = SingInd;
        *Tracestpend = TraceStepsInd;
        
        SpineTL->thetracexists = TRUE;
        if (tlp1 < n) {
            Spine[tlp1].ccstart = TraceCCInd;
            Spine[tlp1].singstart = Spine[tlp1].singend = SingInd;
            Spine[tlp1].stpstart = TraceStepsInd;
            Spine[tlp1].thetracexists = FALSE;
            Spine[tlp1].part->code = -1;
        }
        return 2;
    }
    else {
        if (TraceInd < *TraceEnd) {
            return 0;
        }
        if (Cand->code > SpineTL->part->code) {
            return 2;
        }
        else {
            if (Cand->code < SpineTL->part->code) {
                return 0;
            }
        }
        return 1;
    }
}

void traces_refine_notrace(Candidate *Cand,
                           int n,
                           Partition *Part,
                           struct TracesVars* tv,
                           struct TracesInfo *ti) {
    int i, j, k, jk, sc, ind0, ind1, ind2, ind3, labi, auxcells;
    int value, iend, newcell;
    int HitClsInd, SplInd, SplCntInd, CStackInd;
    int iend1int, j1int;
    unsigned int longcode;
    int Split = 0;
    int Sparse = TRUE;
    int *lab, *cls, *InvLab, *SplitCell, *LabCell;
    int BigCell, BigCellPos, BigCellSize;
    int *nghb;
    const int variation = 1;
    int currentweight, weightstart, weightend, currentcell, currentsize;
    
    if (tv->stackmark > (NAUTY_INFINITY-2)) {
        memset(StackMarkers, 0, n*sizeof(int));
        tv->stackmark = 0;
    }
    tv->stackmark++;
    
    tv->augmented_cells = Part->cells;
    
    lab = Cand->lab;
    InvLab = Cand->invlab;
    cls = Part->cls;
    
    CStackInd = 1;
    CStack[1] = tv->tcellexpath+cls[tv->tcellexpath];
    
    for (i = 1; i <= CStackInd; i++) {
        StackMarkers[CStack[i]] = tv->stackmark;
    }
    
    longcode = Part->cells;
    
    while (CStackInd > 0) {
        
        weightend = 0;
        
        if (tv->mark > (NAUTY_INFINITY-2)) {
            memset(Markers, 0, n*sizeof(int));
            memset(MarkHitVtx, 0, n*sizeof(int));
            tv->mark = 0;
        }
        tv->mark++;
        
        k = Select_from_CStack(cls, CStackInd);
        
        currentcell = CStack[k];
        currentsize = currentcell+cls[currentcell];
        CStack[k] = CStack[CStackInd--];		/* Current Cell */
        longcode = MASHNONCOMM(longcode, currentcell);
        StackMarkers[currentcell] = 0;
        
        labi = lab[currentcell];
        iend1int = TheGraph[labi].d;
        nghb =  TheGraph[labi].e;
        
        do {
            ind0 = currentcell;
            ind2 = currentsize;
            weightstart = weightend;
            
            if (tv->options->weighted) {
                currentweight = (TheGraph[labi].w)[weightstart];
                while ((iend1int > weightend) && ((TheGraph[labi].w)[weightend] == currentweight)) {
                    weightend++;
                }
            } else {
                weightend = TheGraph[labi].d;
            }
            
            /* Analysis of occurrences of neighbors of the current cell */
            /* The list of cells with neighbors in the current cell is  built */
            
            if (cls[ind0] == 1) {   /* SINGLETON CURRENT CELL CASE */
                
                /* NEIGHCOUNT_SING_MULT */
                HitClsInd = 0;
                for (j1int = weightstart; j1int < weightend; ++j1int) {
                    k = nghb[j1int];
                    value = Part->inv[InvLab[k]];
                    if (cls[value] > 1) {
                        if (Markers[value] != tv->mark) {
                            HitCls[HitClsInd++] = value;
                            Markers[value] = tv->mark;
                            ElmHitCll[value] = value;
                        }
                        HitVtx[ElmHitCll[value]++] = k;
                    } else {
                        switch (variation) {
                            case 1:
                                longcode = MASHCOMM(longcode, value);
                                break;
                            default:
                                break;
                        }
                    }
                }
                /* end NEIGHCOUNT_SING_MULT */
                
                tv->mark++;
                FIND_SPLIT_CELLS;
                
                /* Sorting the cells to be split */
                sort_Split_Array(SplCls, SplInd);
                
                /* REARRANGE THE CELLS */
                for (j = 0; j < SplInd; j++) {
                    /* REARRANGE_CELLS */
                    ind1 = SplCls[j];
                    cls[ind1] = cls[ind1]-ElmHitCll[ind1];
                    newcell = ind1+cls[ind1];
                    cls[newcell] = ElmHitCll[ind1];
                    Part->cells++;
                    if (StackMarkers[ind1] != tv->stackmark) {
                        if (cls[newcell] < cls[ind1]) {
                            CStack[++CStackInd] = newcell;
                            StackMarkers[newcell] = tv->stackmark;
                        }
                        else {
                            CStack[++CStackInd] = ind1;
                            StackMarkers[ind1] = tv->stackmark;
                        }
                    }
                    else {
                        CStack[++CStackInd] = newcell;
                        StackMarkers[newcell] = tv->stackmark;
                    }
                    SplitCell = HitVtx+ind1;
                    ind3 = cls[newcell];
                    LabCell = lab+newcell;
                    for (jk = 0; jk < ind3; jk++) {
                        k = SplitCell[jk];
                        i = LabCell[jk];
                        Part->inv[newcell+jk] = newcell;
                        lab[InvLab[k]] = i;
                        InvLab[i] = InvLab[k];
                        LabCell[jk] = k;
                        InvLab[k] = newcell+jk;
                    }
                    /* END REARRANGE_CELLS */
                    
                    if (cls[ind1] == 1) {
                        Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[ind1]);
                    }
                    if (cls[newcell] == 1) {
                        Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[newcell]);
                    }
                }
            }
            else {
                if (ti->thegraphisparse) {
                    
                    /* NEIGHCOUNT_SPARSE_MULT */
                    if (cls[ind0] != n) {
                        HitClsInd = 0;
                        for (i = ind0; i < ind2; i++) {
                            labi = lab[i];
                            nghb = TheGraph[labi].e;
                            for (j = weightstart; j < weightend; j++) {
                                k = nghb[j];
                                if (MarkHitVtx[k] == tv->mark) {
                                    NghCounts[k]++;
                                }
                                else {
                                    value = Part->inv[InvLab[k]];
                                    if (cls[value] > 1) {
                                        MarkHitVtx[k] = tv->mark;
                                        NghCounts[k] = 1;
                                        if (Markers[value] != tv->mark) {
                                            HitCls[HitClsInd++] = value;
                                            Markers[value] = tv->mark;
                                            HitVtx[value] = k;
                                            ElmHitCll[value] = 1;
                                        }
                                        else {
                                            HitVtx[value+ElmHitCll[value]++] = k;
                                        }
                                    }
                                    else {
                                        switch (variation) {
                                            case 1:
                                                longcode = MASHCOMM(longcode, value);
                                                break;
                                            default:
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /* End NEIGHCOUNT_SPARSE_MULT */
                    
                    tv->mark++;
                    
                    SplInd = 0;
                    SplCls[0] = n;
                    for (j = 0; j < HitClsInd; j++) {
                        ind1 = HitCls[j];
                        if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
                            SplCls[SplInd++] = ind1;
                        }
                        else {
                            ind2 = ind1+cls[ind1];
                            Split = FALSE;
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    Split = TRUE;
                                    break;
                                }
                            }
                            if (!Split) {
                                longcode = MASHCOMM(longcode, ind1);
                            }
                        }
                    }
                    
                    /* Sorting the cells to be split */
                    sort_Split_Array(SplCls, SplInd);
                    
                    for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
                        ind0 = SplCls[sc];
                        ind1 = ind0 + cls[ind0];
                        SplCntInd = 0;
                        if (ElmHitCll[ind0] < cls[ind0]) {
                            SplCnt[SplCntInd++] = 0;
                            SplPos[0] = cls[ind0] - ElmHitCll[ind0];
                        }
                        
                        /* According to the numbers of neighbors of C into the current cell */
                        /* compute how many vertices in C will be placed into the same new cell */
                        iend = ind0 + ElmHitCll[ind0];
                        for (i = ind0; i < iend; i++) {
                            value = NghCounts[HitVtx[i]];
                            if (Markers[value] != tv->mark) {
                                Markers[value] = tv->mark;
                                SplCnt[SplCntInd++] = value;
                                SplPos[value] = 1;
                            }
                            else {
                                SplPos[value]++;
                            }
                        }
                        tv->mark++;
                        
                        /* Sort the values deriving from the previous step */
                        sort_Split_Array(SplCnt, SplCntInd);
                        Part->cells += SplCntInd-1;
                        
                        /* Split the cell C and update the information for sizes of new cells */
                        /* Put the new cells into the stack */
                        i = ind0;
                        if (StackMarkers[i] != tv->stackmark) {
                            BigCellSize = 0;
                        }
                        for (k = 0; k < SplCntInd; k++) {
                            value = SplPos[SplCnt[k]];
                            cls[i] = value;
                            if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                BigCell = i;
                                BigCellPos = CStackInd;
                                BigCellSize = cls[i];
                            }
                            SplPos[SplCnt[k]] = i;
                            i += value;
                            if (i < ind1) {
                                CStack[++CStackInd] = i;
                                StackMarkers[i] = tv->stackmark;
                            }
                        }
                        
                        if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                            CStack[BigCellPos] = ind0;
                            StackMarkers[BigCell] = 0;
                            StackMarkers[ind0] = tv->stackmark;
                        }
                        
                        /* Permute elements of the cell C */
                        iend = ind0 + ElmHitCll[ind0];
                        for (i = ind0; i < iend; i++) {
                            value = HitVtx[i];
                            j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
                            k = InvLab[value];				/* where HitVtx[i] is in lab */
                            lab[k] = lab[j];
                            lab[j] = value;
                            InvLab[value] = j;
                            InvLab[lab[k]] = k;
                            NghCounts[value] = 0;
                        }
                        
                        /* Reconstruct the cell C and update the inverse partition */
                        newcell = ind1 - ElmHitCll[ind0];
                        i = newcell;
                        ind2 = newcell+cls[newcell]-1;
                        do {
                            Part->inv[i] = newcell;
                            if (i == ind2) {
                                newcell = i+1;
                                if (newcell < n) ind2 = newcell+cls[newcell]-1;
                            }
                        }
                        while (++i < ind1);
                        
                        for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
                            if (cls[i] == 1) {
                                Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[i]);
                            }
                        }
                        
                    }
                }
                else {
                    if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
                        Sparse = FALSE;
                    }
                    else {
                        Sparse = TRUE;
                    }
                    if (Sparse) {
                        /* Counting occurrences of neighbors of the current cell */
                        /* The list of cells with neighbors in the current cell is also built */
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                            HitCls[0] = 0;
                            HitClsInd = 1;
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_SPARSE_MULT */
                            HitClsInd = 0;
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb =  TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                    value = Part->inv[InvLab[k]];
                                    if (Markers[value] != tv->mark) {
                                        if (cls[value] > 1) HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                    }
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_SPARSE_MULT */
                            
                        }
                        
                        tv->mark++;
                        
                        SplInd = 0;
                        for (j = 0; j < HitClsInd; j++) {
                            ind1 = HitCls[j];
                            ind2 = ind1+cls[ind1];
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    break;
                                }
                            }
                        }
                        
                        /* Sorting the cells to be split */
                        sort_Split_Array(SplCls, SplInd);
                        
                        for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
                            ind0 = SplCls[j];
                            ind1 = ind0+cls[ind0];
                            SplCntInd = 0;
                            
                            /* According to the numbers of neighbors of C into the current cell */
                            /* compute how many vertices in C will be placed into the same new cell */
                            for (i = ind0; i < ind1; i++) {
                                value = NghCounts[lab[i]];
                                if (Markers[value] != tv->mark) {
                                    Markers[value] = tv->mark;
                                    SplCnt[SplCntInd++] = value;
                                    SplPos[value] = 1;
                                }
                                else {
                                    SplPos[value]++;
                                }
                            }
                            
                            tv->mark++;
                            
                            /* Sort the values deriving from the previous step */
                            sort_Split_Array(SplCnt, SplCntInd);
                            
                            Part->cells += SplCntInd-1;
                            
                            /* Split the cell C and update the information for sizes of new cells */
                            /* Put the new cells into the stack */
                            i = ind0;
                            if (StackMarkers[i] != tv->stackmark) {
                                BigCellSize = 0;
                            }
                            for (k = 0; k < SplCntInd; k++) {
                                value = SplPos[SplCnt[k]];
                                cls[i] = value;
                                if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                    BigCell = i;
                                    BigCellPos = CStackInd;
                                    BigCellSize = cls[i];
                                }
                                SplPos[SplCnt[k]] = i;
                                i += value;
                                if (i < ind1) {
                                    CStack[++CStackInd] = i;
                                    StackMarkers[i] = tv->stackmark;
                                }
                            }
                            if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                CStack[BigCellPos] = ind0;
                                StackMarkers[BigCell] = 0;
                                StackMarkers[ind0] = tv->stackmark;
                            }
                            
                            /* Permute elements of the cell C */
                            i = ind0;
                            do {
                                SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                            }
                            while(++i < ind1);
                            
                            /* Reconstruct the cell C and update the inverse partition */
                            newcell = ind0;
                            i = ind0;
                            ind2 = newcell+cls[newcell]-1;
                            do {
                                lab[i] = SplCnt[i];
                                InvLab[lab[i]] = i;
                                Part->inv[i] = newcell;
                                
                                if (i == ind2) {
                                    newcell = i+1;
                                    if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                }
                            }
                            while (++i < ind1);
                            
                            for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
                                if (cls[i] == 1) {
                                    Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[i]);
                                }
                            }
                            
                        }
                    }
                    else {
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_DENSE_MULT */
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb = TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_DENSE_MULT */
                            
                        }
                        
                        ind0 = 0;
                        while (ind0 < n) {	/* For each cell C with size(C) > 1 */
                            ind1 = ind0+cls[ind0];
                            if (cls[ind0] > 1) {
                                
                                /* Determine whether C must be split */
                                SplCntInd = 0;
                                value = NghCounts[lab[ind0]];
                                for (i = ind0+1; i < ind1; i++)
                                {
                                    if (NghCounts[lab[i]] != value)
                                    {
                                        Markers[value] = tv->mark;
                                        SplCnt[SplCntInd++] = value;
                                        SplPos[value] = i-ind0;
                                        do {
                                            value = NghCounts[lab[i++]];
                                            if (Markers[value] != tv->mark) {
                                                Markers[value] = tv->mark;
                                                SplCnt[SplCntInd++] = value;
                                                SplPos[value] = 1;
                                            }
                                            else {
                                                SplPos[value]++;
                                            }
                                        }
                                        while(i != ind1);
                                        break;
                                    }
                                }
                                
                                if (SplCntInd) {
                                    tv->mark++;
                                    
                                    /* Sort the values deriving from the previous step */
                                    sort_Split_Array(SplCnt, SplCntInd);
                                    Part->cells += SplCntInd-1;
                                    
                                    /* Split the cell C and update the information for sizes of new cells */
                                    /* Put the new cells into the stack */
                                    i = ind0;
                                    if (StackMarkers[i] != tv->stackmark) {
                                        BigCellSize = 0;
                                    }
                                    for (k = 0; k < SplCntInd; k++) {
                                        value = SplPos[SplCnt[k]];
                                        cls[i] = value;
                                        if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                            BigCell = i;
                                            BigCellPos = CStackInd;
                                            BigCellSize = cls[i];
                                        }
                                        SplPos[SplCnt[k]] = i;
                                        i += value;
                                        if (i < ind1) {
                                            CStack[++CStackInd] = i;
                                            StackMarkers[i] = tv->stackmark;
                                        }
                                    }
                                    if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                        CStack[BigCellPos] = ind0;
                                        StackMarkers[BigCell] = 0;
                                        StackMarkers[ind0] = tv->stackmark;
                                    }
                                    
                                    /* Permute elements of the cell C */
                                    i = ind0;
                                    do {
                                        SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                                    }
                                    while(++i < ind1);
                                    
                                    /* Reconstruct the cell C and update the inverse partition */
                                    newcell = ind0;
                                    i = ind0;
                                    ind2 = newcell+cls[newcell]-1;
                                    do {
                                        lab[i] = SplCnt[i];
                                        InvLab[lab[i]] = i;
                                        Part->inv[i] = newcell;
                                        if (i == ind2) {
                                            newcell = i+1;
                                            if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                        }
                                    }
                                    while (++i < ind1);
                                    
                                    for (i = ind0; i < ind1; i+=cls[i]) {
                                        if (cls[i] == 1) {
                                            Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[i]);
                                        }
                                    }
                                    
                                }
                            }
                            ind0 = ind1;
                        }
                    }
                }
            }
        }
        while (weightend < iend1int);
    }  /* end while (CStackInd > 0) */
    
    tv->augmented_cells = Part->cells - tv->augmented_cells;
    Cand->code = CLEANUP(longcode);
    return;
}

void traces_refine_maketrie(Candidate *Cand,
                            int n,
                            Partition *Part,
                            struct TracesVars* tv,
                            struct TracesInfo *ti) {
    int i, j, k, jk, sc, ind0, ind1, ind2, ind3, labi;
    int value, iend, newcell;
    int HitClsInd, SplInd, SplCntInd, CStackInd;
    int j1int, iend1int;
    unsigned int longcode;
    int Split = 0;
    int Sparse = TRUE;
    int *lab, *cls, *InvLab, *SplitCell, *LabCell;
    int BigCell, BigCellPos, BigCellSize;
    int *nghb;
    const int variation = 1;
    int currentweight, weightstart, weightend, currentcell, currentsize;
    
    if (tv->stackmark > (NAUTY_INFINITY-2)) {
        memset(StackMarkers, 0, n*sizeof(int));
        tv->stackmark = 0;
    }
    tv->stackmark++;
    
    tv->augmented_cells = Part->cells;
    
    lab = Cand->lab;
    InvLab = Cand->invlab;
    cls = Part->cls;
    
    CStack[1] = Spine[tv->tolevel_tl].tgtpos;
    CStackInd = 1;
    for (i = 1; i <= CStackInd; i++) {
        StackMarkers[CStack[i]] = tv->stackmark;
    }
    
    longcode = Part->cells;
    
    while (CStackInd > 0)
    {
        weightend = 0;
        
        if (tv->mark > (NAUTY_INFINITY-2)) {
            memset(Markers, 0, n*sizeof(int));
            memset(MarkHitVtx, 0, n*sizeof(int));
            tv->mark = 0;
        }
        tv->mark++;
        
        if (Part->cells == n) break;
        
        k = Select_from_CStack(cls, CStackInd);
        
        currentcell = CStack[k];
        currentsize = currentcell+cls[currentcell];
        CStack[k] = CStack[CStackInd--];
        longcode = MASHNONCOMM(longcode, currentcell);
        StackMarkers[currentcell] = 0;
        
        labi = lab[currentcell];
        iend1int = TheGraph[labi].d;
        nghb =  TheGraph[labi].e;
        
        do {
            ind0 = currentcell;
            ind2 = currentsize;
            weightstart = weightend;
            
            if (tv->options->weighted) {
                currentweight = (TheGraph[labi].w)[weightstart];
                while ((iend1int > weightend) && ((TheGraph[labi].w)[weightend] == currentweight)) {
                    weightend++;
                }
            } else {
                weightend = TheGraph[labi].d;
            }
            
            /* Analysis of occurrences of neighbors of the current cell */
            /* The list of cells with neighbors in the current cell is  built */
            if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
                
                /* NEIGHCOUNT_SING_MULT */
                HitClsInd = 0;
                for (j1int = weightstart; j1int < weightend; ++j1int) {
                    k = nghb[j1int];
                    value = Part->inv[InvLab[k]];
                    if (cls[value] > 1) {
                        if (Markers[value] != tv->mark) {
                            HitCls[HitClsInd++] = value;
                            Markers[value] = tv->mark;
                            ElmHitCll[value] = value;
                        }
                        HitVtx[ElmHitCll[value]++] = k;
                    } else {
                        switch (variation) {
                            case 1:
                                longcode = MASHCOMM(longcode, value);
                                break;
                            default:
                                break;
                        }
                    }
                }
                /* end NEIGHCOUNT_SING_MULT */
                
                tv->mark++;
                FIND_SPLIT_CELLS;
                
                /* Sorting the cells to be split */
                sort_Split_Array(SplCls, SplInd);
                
                for (j = 0; j < SplInd; j++) {
                    ind1 = SplCls[j];
                    i = ind1+cls[ind1]-ElmHitCll[ind1];
                    trieref = trie_make(trieref, i, n, tv);
                }
                
                /* REARRANGE THE CELLS */
                for (j = 0; j < SplInd; j++) {
                    /* REARRANGE_CELLS */
                    ind1 = SplCls[j];
                    cls[ind1] = cls[ind1]-ElmHitCll[ind1];
                    newcell = ind1+cls[ind1];
                    cls[newcell] = ElmHitCll[ind1];
                    Part->cells++;
                    if (StackMarkers[ind1] != tv->stackmark) {
                        if (cls[newcell] < cls[ind1]) {
                            CStack[++CStackInd] = newcell;
                            StackMarkers[newcell] = tv->stackmark;
                        }
                        else {
                            CStack[++CStackInd] = ind1;
                            StackMarkers[ind1] = tv->stackmark;
                        }
                    }
                    else {
                        CStack[++CStackInd] = newcell;
                        StackMarkers[newcell] = tv->stackmark;
                    }
                    SplitCell = HitVtx+ind1;
                    ind3 = cls[newcell];
                    LabCell = lab+newcell;
                    for (jk = 0; jk < ind3; jk++) {
                        k = SplitCell[jk];
                        i = LabCell[jk];
                        Part->inv[newcell+jk] = newcell;
                        lab[InvLab[k]] = i;
                        InvLab[i] = InvLab[k];
                        LabCell[jk] = k;
                        InvLab[k] = newcell+jk;
                    }
                    /* END REARRANGE_CELLS */
                }
            }
            else {
                if (ti->thegraphisparse) {
                    
                    /* NEIGHCOUNT_SPARSE_MULT */
                    if (cls[ind0] != n) {
                        HitClsInd = 0;
                        for (i = ind0; i < ind2; i++) {
                            labi = lab[i];
                            nghb = TheGraph[labi].e;
                            for (j = weightstart; j < weightend; j++) {
                                k = nghb[j];
                                if (MarkHitVtx[k] == tv->mark) {
                                    NghCounts[k]++;
                                }
                                else {
                                    value = Part->inv[InvLab[k]];
                                    if (cls[value] > 1) {
                                        MarkHitVtx[k] = tv->mark;
                                        NghCounts[k] = 1;
                                        if (Markers[value] != tv->mark) {
                                            HitCls[HitClsInd++] = value;
                                            Markers[value] = tv->mark;
                                            HitVtx[value] = k;
                                            ElmHitCll[value] = 1;
                                        }
                                        else {
                                            HitVtx[value+ElmHitCll[value]++] = k;
                                        }
                                    }
                                    else {
                                        switch (variation) {
                                            case 1:
                                                longcode = MASHCOMM(longcode, value);
                                                break;
                                            default:
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /* End NEIGHCOUNT_SPARSE_MULT */
                    
                    tv->mark++;
                    
                    SplInd = 0;
                    SplCls[0] = n;
                    for (j = 0; j < HitClsInd; j++) {
                        ind1 = HitCls[j];
                        if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
                            SplCls[SplInd++] = ind1;
                        }
                        else {
                            ind2 = ind1+cls[ind1];
                            Split = FALSE;
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    Split = TRUE;
                                    break;
                                }
                            }
                            if (!Split) {
                                longcode = MASHCOMM(longcode, ind1);
                            }
                        }
                    }
                    
                    /* Sorting the cells to be split */
                    sort_Split_Array(SplCls, SplInd);
                    
                    for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
                        ind0 = SplCls[sc];
                        ind1 = ind0 + cls[ind0];
                        SplCntInd = 0;
                        if (ElmHitCll[ind0] < cls[ind0]) {
                            SplCnt[SplCntInd++] = 0;
                            SplPos[0] = cls[ind0] - ElmHitCll[ind0];
                        }
                        
                        /* According to the numbers of neighbors of C into the current cell */
                        /* compute how many vertices in C will be placed into the same new cell */
                        iend = ind0 + ElmHitCll[ind0];
                        for (i = ind0; i < iend; i++) {
                            value = NghCounts[HitVtx[i]];
                            if (Markers[value] != tv->mark) {
                                Markers[value] = tv->mark;
                                SplCnt[SplCntInd++] = value;
                                SplPos[value] = 1;
                            }
                            else {
                                SplPos[value]++;
                            }
                        }
                        
                        tv->mark++;
                        
                        /* Sort the values deriving from the previous step */
                        sort_Split_Array(SplCnt, SplCntInd);
                        Part->cells += SplCntInd-1;
                        
                        /* Split the cell C and update the information for sizes of new cells */
                        /* Put the new cells into the stack */
                        i = ind0;
                        if (StackMarkers[i] != tv->stackmark) {
                            BigCellSize = 0;
                        }
                        for (k = 0; k < SplCntInd; k++) {
                            value = SplPos[SplCnt[k]];
                            cls[i] = value;
                            if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                BigCell = i;
                                BigCellPos = CStackInd;
                                BigCellSize = cls[i];
                            }
                            SplPos[SplCnt[k]] = i;
                            i += value;
                            if (i < ind1) {
                                CStack[++CStackInd] = i;
                                StackMarkers[i] = tv->stackmark;
                                trieref = trie_make(trieref, i, n, tv);
                            }
                        }
                        
                        if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                            CStack[BigCellPos] = ind0;
                            StackMarkers[BigCell] = 0;
                            StackMarkers[ind0] = tv->stackmark;
                        }
                        /* Permute elements of the cell C */
                        iend = ind0 + ElmHitCll[ind0];
                        for (i = ind0; i < iend; i++) {
                            value = HitVtx[i];
                            j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
                            k = InvLab[value];				/* where HitVtx[i] is in lab */
                            lab[k] = lab[j];
                            lab[j] = value;
                            InvLab[value] = j;
                            InvLab[lab[k]] = k;
                            NghCounts[value] = 0;
                        }
                        
                        /* Reconstruct the cell C and update the inverse partition */
                        newcell = ind1 - ElmHitCll[ind0];
                        i = newcell;
                        ind2 = newcell+cls[newcell]-1;
                        do {
                            Part->inv[i] = newcell;
                            if (i == ind2) {
                                newcell = i+1;
                                if (newcell < n) ind2 = newcell+cls[newcell]-1;
                            }
                        }
                        while (++i < ind1);
                    }
                }
                else {
                    if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
                        Sparse = FALSE;
                    }
                    else {
                        Sparse = TRUE;
                    }
                    if (Sparse) {
                        /* Counting occurrences of neighbors of the current cell */
                        /* The list of cells with neighbors in the current cell is also built */
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                            HitCls[0] = 0;
                            HitClsInd = 1;
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_SPARSE_MULT */
                            HitClsInd = 0;
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb =  TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                    value = Part->inv[InvLab[k]];
                                    if (Markers[value] != tv->mark) {
                                        if (cls[value] > 1) HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                    }
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_SPARSE_MULT */
                            
                        }
                        
                        tv->mark++;
                        
                        SplInd = 0;
                        for (j = 0; j < HitClsInd; j++) {
                            ind1 = HitCls[j];
                            ind2 = ind1+cls[ind1];
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    break;
                                }
                            }
                        }
                        
                        /* Sorting the cells to be split */
                        sort_Split_Array(SplCls, SplInd);
                        
                        for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
                            ind0 = SplCls[j];
                            ind1 = ind0+cls[ind0];
                            SplCntInd = 0;
                            
                            /* According to the numbers of neighbors of C into the current cell */
                            /* compute how many vertices in C will be placed into the same new cell */
                            for (i = ind0; i < ind1; i++) {
                                value = NghCounts[lab[i]];
                                if (Markers[value] != tv->mark) {
                                    Markers[value] = tv->mark;
                                    SplCnt[SplCntInd++] = value;
                                    SplPos[value] = 1;
                                }
                                else {
                                    SplPos[value]++;
                                }
                            }
                            
                            tv->mark++;
                            
                            /* Sort the values deriving from the previous step */
                            sort_Split_Array(SplCnt, SplCntInd);
                            Part->cells += SplCntInd-1;
                            
                            /* Split the cell C and update the information for sizes of new cells */
                            /* Put the new cells into the stack */
                            i = ind0;
                            if (StackMarkers[i] != tv->stackmark) {
                                BigCellSize = 0;
                            }
                            for (k = 0; k < SplCntInd; k++) {
                                value = SplPos[SplCnt[k]];
                                cls[i] = value;
                                if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                    BigCell = i;
                                    BigCellPos = CStackInd;
                                    BigCellSize = cls[i];
                                }
                                SplPos[SplCnt[k]] = i;
                                i += value;
                                if (i < ind1) {
                                    CStack[++CStackInd] = i;
                                    StackMarkers[i] = tv->stackmark;
                                    trieref = trie_make(trieref, i, n, tv);
                                }
                            }
                            if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                CStack[BigCellPos] = ind0;
                                StackMarkers[BigCell] = 0;
                                StackMarkers[ind0] = tv->stackmark;
                            }
                            
                            /* Permute elements of the cell C */
                            i = ind0;
                            do {
                                SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                            }
                            while(++i < ind1);
                            
                            /* Reconstruct the cell C and update the inverse partition */
                            newcell = ind0;
                            i = ind0;
                            ind2 = newcell+cls[newcell]-1;
                            do {
                                lab[i] = SplCnt[i];
                                InvLab[lab[i]] = i;
                                Part->inv[i] = newcell;
                                if (i == ind2) {
                                    newcell = i+1;
                                    if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                }
                            }
                            while (++i < ind1);
                        }
                    }
                    else {
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_DENSE_MULT */
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb = TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_DENSE_MULT */
                            
                        }
                        
                        ind0 = 0;
                        while (ind0 < n) {	/* For each cell C with size(C) > 1 */
                            ind1 = ind0+cls[ind0];
                            if (cls[ind0] > 1) {
                                
                                /* Determine whether C must be split */
                                SplCntInd = 0;
                                value = NghCounts[lab[ind0]];
                                for (i = ind0+1; i < ind1; i++)
                                {
                                    if (NghCounts[lab[i]] != value)
                                    {
                                        Markers[value] = tv->mark;
                                        SplCnt[SplCntInd++] = value;
                                        SplPos[value] = i-ind0;
                                        do {
                                            value = NghCounts[lab[i++]];
                                            if (Markers[value] != tv->mark) {
                                                Markers[value] = tv->mark;
                                                SplCnt[SplCntInd++] = value;
                                                SplPos[value] = 1;
                                            }
                                            else {
                                                SplPos[value]++;
                                            }
                                        }
                                        while(i != ind1);
                                        break;
                                    }
                                }
                                
                                if (SplCntInd) {
                                    tv->mark++;
                                    
                                    /* Sort the values deriving from the previous step */
                                    sort_Split_Array(SplCnt, SplCntInd);
                                    
                                    Part->cells += SplCntInd-1;
                                    
                                    /* Split the cell C and update the information for sizes of new cells */
                                    /* Put the new cells into the stack */
                                    i = ind0;
                                    if (StackMarkers[i] != tv->stackmark) {
                                        BigCellSize = 0;
                                    }
                                    for (k = 0; k < SplCntInd; k++) {
                                        value = SplPos[SplCnt[k]];
                                        cls[i] = value;
                                        if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                            BigCell = i;
                                            BigCellPos = CStackInd;
                                            BigCellSize = cls[i];
                                        }
                                        SplPos[SplCnt[k]] = i;
                                        i += value;
                                        if (i < ind1) {
                                            CStack[++CStackInd] = i;
                                            StackMarkers[i] = tv->stackmark;
                                        }
                                    }
                                    if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                        CStack[BigCellPos] = ind0;
                                        StackMarkers[BigCell] = 0;
                                        StackMarkers[ind0] = tv->stackmark;
                                        trieref = trie_make(trieref, i, n, tv);
                                    }
                                    
                                    /* Permute elements of the cell C */
                                    i = ind0;
                                    do {
                                        SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                                    }
                                    while(++i < ind1);
                                    
                                    /* Reconstruct the cell C and update the inverse partition */
                                    newcell = ind0;
                                    i = ind0;
                                    ind2 = newcell+cls[newcell]-1;
                                    do {
                                        lab[i] = SplCnt[i];
                                        InvLab[lab[i]] = i;
                                        Part->inv[i] = newcell;
                                        if (i == ind2) {
                                            newcell = i+1;
                                            if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                        }
                                    }
                                    while (++i < ind1);
                                }
                            }
                            ind0 = ind1;
                        }
                    }
                }
            }
        }
        while (weightend < iend1int);
    } /* end while (CStackInd > 0) */
    
    tv->augmented_cells = Part->cells - tv->augmented_cells;
    
    Cand->code = CLEANUP(longcode);
    return;
}

int traces_refine_comptrie(Candidate *Cand,
                           int n,
                           Partition *Part,
                           struct TracesVars* tv,
                           struct TracesInfo *ti) {
    int i, j, k, jk, sc, ind0, ind1, ind2, ind3, labi;
    int value, iend, newcell;
    int HitClsInd, SplInd, SplCntInd, CStackInd;
    int j1int, iend1int;
    unsigned int longcode;
    int Split = 0;
    int Sparse = TRUE;
    int *lab, *cls, *InvLab, *SplitCell, *LabCell;
    int BigCell, BigCellPos, BigCellSize;
    int *nghb;
    const int variation = 1;
    int currentweight, weightstart, weightend, currentcell, currentsize;
    
    if (tv->stackmark > (NAUTY_INFINITY-2)) {
        memset(StackMarkers, 0, n*sizeof(int));
        tv->stackmark = 0;
    }
    tv->stackmark++;
    
    tv->augmented_cells = Part->cells;
    
    lab = Cand->lab;
    InvLab = Cand->invlab;
    cls = Part->cls;
    
    CStack[1] = Spine[tv->tolevel_tl].tgtpos;
    CStackInd = 1;
    for (i = 1; i <= CStackInd; i++) {
        StackMarkers[CStack[i]] = tv->stackmark;
    }
    
    longcode = Part->cells;
    while (CStackInd > 0)
    {
        
        weightend = 0;
        
        if (tv->mark > (NAUTY_INFINITY-2)) {
            memset(Markers, 0, n*sizeof(int));
            memset(MarkHitVtx, 0, n*sizeof(int));
            tv->mark = 0;
        }
        tv->mark++;
        
        if (Part->cells == n) break;
        
        k = Select_from_CStack(cls, CStackInd);
        
        currentcell = CStack[k];
        currentsize = currentcell+cls[currentcell];
        CStack[k] = CStack[CStackInd--];		/* Current Cell */
        longcode = MASHNONCOMM(longcode, currentcell);
        StackMarkers[currentcell] = 0;
        
        labi = lab[currentcell];
        iend1int = TheGraph[labi].d;
        nghb =  TheGraph[labi].e;
        
        do {
            ind0 = currentcell;
            ind2 = currentsize;
            weightstart = weightend;
            
            if (tv->options->weighted) {
                currentweight = (TheGraph[labi].w)[weightstart];
                while ((iend1int > weightend) && ((TheGraph[labi].w)[weightend] == currentweight)) {
                    weightend++;
                }
            } else {
                weightend = TheGraph[labi].d;
            }
            
            /* Analysis of occurrences of neighbors of the current cell */
            /* The list of cells with neighbors in the current cell is  built */
            if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
                
                /* NEIGHCOUNT_SING_MULT */
                HitClsInd = 0;
                for (j1int = weightstart; j1int < weightend; ++j1int) {
                    k = nghb[j1int];
                    value = Part->inv[InvLab[k]];
                    if (cls[value] > 1) {
                        if (Markers[value] != tv->mark) {
                            HitCls[HitClsInd++] = value;
                            Markers[value] = tv->mark;
                            ElmHitCll[value] = value;
                        }
                        HitVtx[ElmHitCll[value]++] = k;
                    } else {
                        switch (variation) {
                            case 1:
                                longcode = MASHCOMM(longcode, value);
                                break;
                            default:
                                break;
                        }
                    }
                }
                /* end NEIGHCOUNT_SING_MULT */
                
                tv->mark++;
                FIND_SPLIT_CELLS;
                
                /* Sorting the cells to be split */
                sort_Split_Array(SplCls, SplInd);
                
                for (j = 0; j < SplInd; j++) {
                    ind1 = SplCls[j];
                    i = ind1+cls[ind1]-ElmHitCll[ind1];
                    trieref = trie_comp(trieref, i);
                    if (trieref == NULL) return 0;
                }
                
                /* REARRANGE THE CELLS */
                for (j = 0; j < SplInd; j++) {
                    /* REARRANGE_CELLS */
                    ind1 = SplCls[j];
                    cls[ind1] = cls[ind1]-ElmHitCll[ind1];
                    newcell = ind1+cls[ind1];
                    cls[newcell] = ElmHitCll[ind1];
                    Part->cells++;
                    if (StackMarkers[ind1] != tv->stackmark) {
                        if (cls[newcell] < cls[ind1]) {
                            CStack[++CStackInd] = newcell;
                            StackMarkers[newcell] = tv->stackmark;
                        }
                        else {
                            CStack[++CStackInd] = ind1;
                            StackMarkers[ind1] = tv->stackmark;
                        }
                    }
                    else {
                        CStack[++CStackInd] = newcell;
                        StackMarkers[newcell] = tv->stackmark;
                    }
                    SplitCell = HitVtx+ind1;
                    ind3 = cls[newcell];
                    LabCell = lab+newcell;
                    for (jk = 0; jk < ind3; jk++) {
                        k = SplitCell[jk];
                        i = LabCell[jk];
                        Part->inv[newcell+jk] = newcell;
                        lab[InvLab[k]] = i;
                        InvLab[i] = InvLab[k];
                        LabCell[jk] = k;
                        InvLab[k] = newcell+jk;
                    }
                    /* END REARRANGE_CELLS */
                }
            }
            else {
                if (ti->thegraphisparse) {
                    
                    /* NEIGHCOUNT_SPARSE_MULT */
                    if (cls[ind0] != n) {
                        HitClsInd = 0;
                        for (i = ind0; i < ind2; i++) {
                            labi = lab[i];
                            nghb = TheGraph[labi].e;
                            for (j = weightstart; j < weightend; j++) {
                                k = nghb[j];
                                if (MarkHitVtx[k] == tv->mark) {
                                    NghCounts[k]++;
                                }
                                else {
                                    value = Part->inv[InvLab[k]];
                                    if (cls[value] > 1) {
                                        MarkHitVtx[k] = tv->mark;
                                        NghCounts[k] = 1;
                                        if (Markers[value] != tv->mark) {
                                            HitCls[HitClsInd++] = value;
                                            Markers[value] = tv->mark;
                                            HitVtx[value] = k;
                                            ElmHitCll[value] = 1;
                                        }
                                        else {
                                            HitVtx[value+ElmHitCll[value]++] = k;
                                        }
                                    }
                                    else {
                                        switch (variation) {
                                            case 1:
                                                longcode = MASHCOMM(longcode, value);
                                                break;
                                            default:
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /* End NEIGHCOUNT_SPARSE_MULT */
                    
                    tv->mark++;
                    
                    SplInd = 0;
                    SplCls[0] = n;
                    for (j = 0; j < HitClsInd; j++) {
                        ind1 = HitCls[j];
                        if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
                            SplCls[SplInd++] = ind1;
                        }
                        else {
                            ind2 = ind1+cls[ind1];
                            Split = FALSE;
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    Split = TRUE;
                                    break;
                                }
                            }
                            if (!Split) {
                                longcode = MASHCOMM(longcode, ind1);
                            }
                        }
                    }
                    
                    /* Sorting the cells to be split */
                    sort_Split_Array(SplCls, SplInd);
                    
                    for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
                        ind0 = SplCls[sc];
                        ind1 = ind0 + cls[ind0];
                        SplCntInd = 0;
                        if (ElmHitCll[ind0] < cls[ind0]) {
                            SplCnt[SplCntInd++] = 0;
                            SplPos[0] = cls[ind0] - ElmHitCll[ind0];
                        }
                        
                        /* According to the numbers of neighbors of C into the current cell */
                        /* compute how many vertices in C will be placed into the same new cell */
                        iend = ind0 + ElmHitCll[ind0];
                        for (i = ind0; i < iend; i++) {
                            value = NghCounts[HitVtx[i]];
                            if (Markers[value] != tv->mark) {
                                Markers[value] = tv->mark;
                                SplCnt[SplCntInd++] = value;
                                SplPos[value] = 1;
                            }
                            else {
                                SplPos[value]++;
                            }
                        }
                        
                        tv->mark++;
                        
                        /* Sort the values deriving from the previous step */
                        sort_Split_Array(SplCnt, SplCntInd);
                        Part->cells += SplCntInd-1;
                        
                        /* Split the cell C and update the information for sizes of new cells */
                        /* Put the new cells into the stack */
                        i = ind0;
                        if (StackMarkers[i] != tv->stackmark) {
                            BigCellSize = 0;
                        }
                        for (k = 0; k < SplCntInd; k++) {
                            value = SplPos[SplCnt[k]];
                            cls[i] = value;
                            if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                BigCell = i;
                                BigCellPos = CStackInd;
                                BigCellSize = cls[i];
                            }
                            SplPos[SplCnt[k]] = i;
                            i += value;
                            if (i < ind1) {
                                CStack[++CStackInd] = i;
                                StackMarkers[i] = tv->stackmark;
                                trieref = trie_comp(trieref, i);
                                if (trieref == NULL) return 0;
                            }
                        }
                        
                        if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                            CStack[BigCellPos] = ind0;
                            StackMarkers[BigCell] = 0;
                            StackMarkers[ind0] = tv->stackmark;
                        }
                        /* Permute elements of the cell C */
                        iend = ind0 + ElmHitCll[ind0];
                        for (i = ind0; i < iend; i++) {
                            value = HitVtx[i];
                            j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
                            k = InvLab[value];				/* where HitVtx[i] is in lab */
                            lab[k] = lab[j];
                            lab[j] = value;
                            InvLab[value] = j;
                            InvLab[lab[k]] = k;
                            NghCounts[value] = 0;
                        }
                        
                        /* Reconstruct the cell C and update the inverse partition */
                        newcell = ind1 - ElmHitCll[ind0];
                        i = newcell;
                        ind2 = newcell+cls[newcell]-1;
                        do {
                            Part->inv[i] = newcell;
                            if (i == ind2) {
                                newcell = i+1;
                                if (newcell < n) ind2 = newcell+cls[newcell]-1;
                            }
                        }
                        while (++i < ind1);
                    }
                }
                else {
                    if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
                        Sparse = FALSE;
                    }
                    else {
                        Sparse = TRUE;
                    }
                    if (Sparse) {
                        /* Counting occurrences of neighbors of the current cell */
                        /* The list of cells with neighbors in the current cell is also built */
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                            HitCls[0] = 0;
                            HitClsInd = 1;
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_SPARSE_MULT */
                            HitClsInd = 0;
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb =  TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                    value = Part->inv[InvLab[k]];
                                    if (Markers[value] != tv->mark) {
                                        if (cls[value] > 1) HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                    }
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_SPARSE_MULT */
                            ;
                        }
                        
                        tv->mark++;
                        
                        SplInd = 0;
                        for (j = 0; j < HitClsInd; j++) {
                            ind1 = HitCls[j];
                            ind2 = ind1+cls[ind1];
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    break;
                                }
                            }
                        }
                        
                        /* Sorting the cells to be split */
                        sort_Split_Array(SplCls, SplInd);
                        
                        for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
                            ind0 = SplCls[j];
                            ind1 = ind0+cls[ind0];
                            SplCntInd = 0;
                            
                            /* According to the numbers of neighbors of C into the current cell */
                            /* compute how many vertices in C will be placed into the same new cell */
                            for (i = ind0; i < ind1; i++) {
                                value = NghCounts[lab[i]];
                                if (Markers[value] != tv->mark) {
                                    Markers[value] = tv->mark;
                                    SplCnt[SplCntInd++] = value;
                                    SplPos[value] = 1;
                                }
                                else {
                                    SplPos[value]++;
                                }
                            }
                            
                            tv->mark++;
                            
                            /* Sort the values deriving from the previous step */
                            sort_Split_Array(SplCnt, SplCntInd);
                            
                            Part->cells += SplCntInd-1;
                            
                            /* Split the cell C and update the information for sizes of new cells */
                            /* Put the new cells into the stack */
                            i = ind0;
                            if (StackMarkers[i] != tv->stackmark) {
                                BigCellSize = 0;
                            }
                            for (k = 0; k < SplCntInd; k++) {
                                value = SplPos[SplCnt[k]];
                                cls[i] = value;
                                if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                    BigCell = i;
                                    BigCellPos = CStackInd;
                                    BigCellSize = cls[i];
                                }
                                SplPos[SplCnt[k]] = i;
                                i += value;
                                if (i < ind1) {
                                    CStack[++CStackInd] = i;
                                    StackMarkers[i] = tv->stackmark;
                                    trieref = trie_comp(trieref, i);
                                    if (trieref == NULL) return 0;
                                }
                            }
                            if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                CStack[BigCellPos] = ind0;
                                StackMarkers[BigCell] = 0;
                                StackMarkers[ind0] = tv->stackmark;
                            }
                            
                            /* Permute elements of the cell C */
                            i = ind0;
                            do {
                                SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                            }
                            while(++i < ind1);
                            
                            /* Reconstruct the cell C and update the inverse partition */
                            newcell = ind0;
                            i = ind0;
                            ind2 = newcell+cls[newcell]-1;
                            do {
                                lab[i] = SplCnt[i];
                                InvLab[lab[i]] = i;
                                Part->inv[i] = newcell;
                                
                                if (i == ind2) {
                                    newcell = i+1;
                                    if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                }
                            }
                            while (++i < ind1);
                        }
                    }
                    else {
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_DENSE_MULT */
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb = TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_DENSE_MULT */
                            
                        }
                        
                        ind0 = 0;
                        while (ind0 < n) {	/* For each cell C with size(C) > 1 */
                            ind1 = ind0+cls[ind0];
                            if (cls[ind0] > 1) {
                                
                                /* Determine whether C must be split */
                                SplCntInd = 0;
                                value = NghCounts[lab[ind0]];
                                for (i = ind0+1; i < ind1; i++)
                                {
                                    if (NghCounts[lab[i]] != value)
                                    {
                                        Markers[value] = tv->mark;
                                        SplCnt[SplCntInd++] = value;
                                        SplPos[value] = i-ind0;
                                        do {
                                            value = NghCounts[lab[i++]];
                                            if (Markers[value] != tv->mark) {
                                                Markers[value] = tv->mark;
                                                SplCnt[SplCntInd++] = value;
                                                SplPos[value] = 1;
                                            }
                                            else {
                                                SplPos[value]++;
                                            }
                                        }
                                        while(i != ind1);
                                        break;
                                    }
                                }
                                
                                if (SplCntInd) {
                                    tv->mark++;
                                    
                                    /* Sort the values deriving from the previous step */
                                    sort_Split_Array(SplCnt, SplCntInd);
                                    
                                    Part->cells += SplCntInd-1;
                                    
                                    /* Split the cell C and update the information for sizes of new cells */
                                    /* Put the new cells into the stack */
                                    i = ind0;
                                    if (StackMarkers[i] != tv->stackmark) {
                                        BigCellSize = 0;
                                    }
                                    for (k = 0; k < SplCntInd; k++) {
                                        value = SplPos[SplCnt[k]];
                                        cls[i] = value;
                                        if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                            BigCell = i;
                                            BigCellPos = CStackInd;
                                            BigCellSize = cls[i];
                                        }
                                        SplPos[SplCnt[k]] = i;
                                        i += value;
                                        if (i < ind1) {
                                            CStack[++CStackInd] = i;
                                            StackMarkers[i] = tv->stackmark;
                                        }
                                    }
                                    if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                        CStack[BigCellPos] = ind0;
                                        StackMarkers[BigCell] = 0;
                                        StackMarkers[ind0] = tv->stackmark;
                                        trieref = trie_comp(trieref, i);
                                        if (trieref == NULL) return 0;
                                    }
                                    
                                    /* Permute elements of the cell C */
                                    i = ind0;
                                    do {
                                        SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                                    }
                                    while(++i < ind1);
                                    
                                    /* Reconstruct the cell C and update the inverse partition */
                                    newcell = ind0;
                                    i = ind0;
                                    ind2 = newcell+cls[newcell]-1;
                                    do {
                                        lab[i] = SplCnt[i];
                                        InvLab[lab[i]] = i;
                                        Part->inv[i] = newcell;
                                        if (i == ind2) {
                                            newcell = i+1;
                                            if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                        }
                                    }
                                    while (++i < ind1);
                                }
                            }
                            ind0 = ind1;
                        }
                    }
                }
            }
        }
        while (weightend < iend1int);
    }  /* end while (CStackInd > 0) */
    
    tv->augmented_cells = Part->cells - tv->augmented_cells;
    
    Cand->code = CLEANUP(longcode);
    return 1;
}

int traces_refine_sametrace(Candidate *Cand,
                            int n,
                            Partition *Part,
                            struct TracesVars* tv,
                            struct TracesInfo *ti) {
    int i, j, k, jk, sc, ind0, ind1, ind2, ind3, ind4, labi;
    int value, iend, newcell;
    int HitClsInd, SplInd, SplCntInd, CStackInd, TraceInd, TraceCCInd, TraceStepsInd;
    int j1int, iend1int;
    unsigned int longcode;
    int Sparse = TRUE;
    int *lab, *cls, *InvLab, *TracePos, *SplitCell, *LabCell, *TraceEnd, Traceccend, *Tracestpend;
    int BigCell, BigCellPos, BigCellSize;
    boolean TraceCell = FALSE;
    int *nghb;
    const int variation = 0;
    int currentweight, weightstart, weightend, currentcell, currentsize;
    
    if (tv->stackmark > (NAUTY_INFINITY-2)) {
        memset(StackMarkers, 0, n*sizeof(int));
        tv->stackmark = 0;
    }
    tv->stackmark++;
    
    tv->augmented_cells = Part->cells;
    
    SpineTL = Spine+tv->tolevel;
    TraceEnd = &(SpineTL->trcend);
    Traceccend = SpineTL->ccend;
    Tracestpend = &(SpineTL->stpend);
    TraceCCInd = SpineTL->ccstart;
    TraceStepsInd = SpineTL->stpstart;
    
    lab = Cand->lab;
    InvLab = Cand->invlab;
    cls = Part->cls;
    
    UPDATEMIN(Part->active, n-1);
    memcpy(CStack+1, TheTrace+SpineTL->trcstart, (Part->active)*sizeof(int));
    CStackInd = Part->active;
    for (i = 1; i <= CStackInd; i++) {
        StackMarkers[CStack[i]] = tv->stackmark;
    }
    
    longcode = Part->cells;
    TraceInd = SpineTL->trcstart+Part->active;
    
    while (CStackInd > 0)
    {
        
        weightend = 0;
        
        if (tv->mark > (NAUTY_INFINITY-2)) {
            memset(Markers, 0, n*sizeof(int));
            memset(MarkHitVtx, 0, n*sizeof(int));
            tv->mark = 0;
        }
        tv->mark++;
        
        if (Part->cells == n) break;
        
        k = Select_from_CStack(cls, CStackInd);
        
        currentcell = CStack[k];
        currentsize = currentcell+cls[currentcell];
        CStack[k] = CStack[CStackInd--];		/* Current Cell */
        StackMarkers[currentcell] = 0;
        
        labi = lab[currentcell];
        iend1int = TheGraph[labi].d;
        nghb =  TheGraph[labi].e;
        
        do {
            ind0 = currentcell;
            ind2 = currentsize;
            weightstart = weightend;
            
            if (tv->options->weighted) {
                currentweight = (TheGraph[labi].w)[weightstart];
                while ((iend1int > weightend) && ((TheGraph[labi].w)[weightend] == currentweight)) {
                    weightend++;
                }
            } else {
                weightend = TheGraph[labi].d;
            }
            
            TraceCell = ((ind0 == TheTraceCC[TraceCCInd]) && (TraceCCInd < Traceccend));
            
            /* Analysis of occurrences of neighbors of the current cell */
            /* The list of cells with neighbors in the current cell is  built */
            if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
                
                /* NEIGHCOUNT_SING_MULT */
                HitClsInd = 0;
                for (j1int = weightstart; j1int < weightend; ++j1int) {
                    k = nghb[j1int];
                    value = Part->inv[InvLab[k]];
                    if (cls[value] > 1) {
                        if (Markers[value] != tv->mark) {
                            HitCls[HitClsInd++] = value;
                            Markers[value] = tv->mark;
                            ElmHitCll[value] = value;
                        }
                        HitVtx[ElmHitCll[value]++] = k;
                    } else {
                        switch (variation) {
                            case 1:
                                longcode = MASHCOMM(longcode, value);
                                break;
                            default:
                                break;
                        }
                    }
                }
                /* end NEIGHCOUNT_SING_MULT */
                
                tv->mark++;
                FIND_SPLIT_CELLS;
                
                /* SINGLETON CC CASE */
                if (SplInd) {
                    SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
                    
                    /* Sorting the cells to be split */
                    sort_Split_Array(SplCls, SplInd);
                    
                    for (j = 0; j < SplInd; j++) {
                        ind1 = SplCls[j];
                        i = ind1+cls[ind1]-ElmHitCll[ind1];
                        SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
                    }
                    
                    /* REARRANGE THE CELLS */
                    for (j = 0; j < SplInd; j++) {
                        /* REARRANGE_CELLS */
                        ind1 = SplCls[j];
                        cls[ind1] = cls[ind1]-ElmHitCll[ind1];
                        newcell = ind1+cls[ind1];
                        cls[newcell] = ElmHitCll[ind1];
                        Part->cells++;
                        if (StackMarkers[ind1] != tv->stackmark) {
                            if (cls[newcell] < cls[ind1]) {
                                CStack[++CStackInd] = newcell;
                                StackMarkers[newcell] = tv->stackmark;
                            }
                            else {
                                CStack[++CStackInd] = ind1;
                                StackMarkers[ind1] = tv->stackmark;
                            }
                        }
                        else {
                            CStack[++CStackInd] = newcell;
                            StackMarkers[newcell] = tv->stackmark;
                        }
                        SplitCell = HitVtx+ind1;
                        ind3 = cls[newcell];
                        LabCell = lab+newcell;
                        for (jk = 0; jk < ind3; jk++) {
                            k = SplitCell[jk];
                            i = LabCell[jk];
                            Part->inv[newcell+jk] = newcell;
                            lab[InvLab[k]] = i;
                            InvLab[i] = InvLab[k];
                            LabCell[jk] = k;
                            InvLab[k] = newcell+jk;
                        }
                        /* END REARRANGE_CELLS */
                        
                        if (cls[ind1] == 1) {
                            Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[ind1]);
                        }
                        if (cls[newcell] == 1) {
                            Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[newcell]);
                        }
                    }
                }
                else {
                    if (TraceCell) {
                        return 0;
                    }
                }
                
            }
            else {
                if (ti->thegraphisparse) {
                    
                    /* NEIGHCOUNT_SPARSE_MULT */
                    if (cls[ind0] != n) {
                        HitClsInd = 0;
                        for (i = ind0; i < ind2; i++) {
                            labi = lab[i];
                            nghb = TheGraph[labi].e;
                            for (j = weightstart; j < weightend; j++) {
                                k = nghb[j];
                                if (MarkHitVtx[k] == tv->mark) {
                                    NghCounts[k]++;
                                }
                                else {
                                    value = Part->inv[InvLab[k]];
                                    if (cls[value] > 1) {
                                        MarkHitVtx[k] = tv->mark;
                                        NghCounts[k] = 1;
                                        if (Markers[value] != tv->mark) {
                                            HitCls[HitClsInd++] = value;
                                            Markers[value] = tv->mark;
                                            HitVtx[value] = k;
                                            ElmHitCll[value] = 1;
                                        }
                                        else {
                                            HitVtx[value+ElmHitCll[value]++] = k;
                                        }
                                    }
                                    else {
                                        switch (variation) {
                                            case 1:
                                                longcode = MASHCOMM(longcode, value);
                                                break;
                                            default:
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /* End NEIGHCOUNT_SPARSE_MULT */
                    
                    tv->mark++;
                    
                    SplInd = 0;
                    SplCls[0] = n;
                    for (j = 0; j < HitClsInd; j++) {
                        ind1 = HitCls[j];
                        if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
                            SplCls[SplInd++] = ind1;
                        }
                        else {
                            ind2 = ind1+cls[ind1];
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    break;
                                }
                            }
                        }
                    }
                    
                    /* SPARSE CASE */
                    if (SplInd) {
                        SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+n, &Traceccend)
                        
                        /* Sorting the cells to be split */
                        sort_Split_Array(SplCls, SplInd);
                        
                        for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
                            ind0 = SplCls[sc];
                            ind1 = ind0 + cls[ind0];
                            SplCntInd = 0;
                            if (ElmHitCll[ind0] < cls[ind0]) {
                                SplCnt[SplCntInd++] = 0;
                                SplPos[0] = cls[ind0] - ElmHitCll[ind0];
                            }
                            
                            /* According to the numbers of neighbors of C into the current cell */
                            /* compute how many vertices in C will be placed into the same new cell */
                            iend = ind0 + ElmHitCll[ind0];
                            for (i = ind0; i < iend; i++) {
                                value = NghCounts[HitVtx[i]];
                                if (Markers[value] != tv->mark) {
                                    Markers[value] = tv->mark;
                                    SplCnt[SplCntInd++] = value;
                                    SplPos[value] = 1;
                                }
                                else {
                                    SplPos[value]++;
                                }
                            }
                            tv->mark++;
                            
                            if (SplCntInd) {
                                SAMETRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+n, Tracestpend)
                            }
                            
                            /* Sort the values deriving from the previous step */
                            sort_Split_Array(SplCnt, SplCntInd);
                            
                            Part->cells += SplCntInd-1;
                            
                            /* Split the cell C and update the information for sizes of new cells */
                            /* Put the new cells into the stack */
                            i = ind0;
                            if (StackMarkers[i] != tv->stackmark) {
                                BigCellSize = 0;
                            }
                            for (k = 0; k < SplCntInd; k++) {
                                value = SplPos[SplCnt[k]];
                                cls[i] = value;
                                if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                    BigCell = i;
                                    BigCellPos = CStackInd;
                                    BigCellSize = cls[i];
                                }
                                SplPos[SplCnt[k]] = i;
                                i += value;
                                if (i < ind1) {
                                    CStack[++CStackInd] = i;
                                    StackMarkers[i] = tv->stackmark;
                                    SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
                                }
                            }
                            
                            if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                CStack[BigCellPos] = ind0;
                                StackMarkers[BigCell] = 0;
                                StackMarkers[ind0] = tv->stackmark;
                            }
                            /* Permute elements of the cell C */
                            iend = ind0 + ElmHitCll[ind0];
                            for (i = ind0; i < iend; i++) {
                                value = HitVtx[i];
                                j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
                                k = InvLab[value];				/* where HitVtx[i] is in lab */
                                lab[k] = lab[j];
                                lab[j] = value;
                                InvLab[value] = j;
                                InvLab[lab[k]] = k;
                                NghCounts[value] = 0;
                            }
                            
                            /* Reconstruct the cell C and update the inverse partition */
                            newcell = ind1 - ElmHitCll[ind0];
                            i = newcell;
                            ind2 = newcell+cls[newcell]-1;
                            do {
                                Part->inv[i] = newcell;
                                if (i == ind2) {
                                    newcell = i+1;
                                    if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                }
                            }
                            while (++i < ind1);
                            
                            for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
                                if (cls[i] == 1) {
                                    Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                }
                            }
                            
                        }
                    }
                    else {
                        if (TraceCell) {
                            return 0;
                        }
                    }
                    
                }
                else {
                    if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
                        Sparse = FALSE;
                    }
                    else {
                        Sparse = TRUE;
                    }
                    if (Sparse) {
                        /* Counting occurrences of neighbors of the current cell */
                        /* The list of cells with neighbors in the current cell is also built */
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                            HitCls[0] = 0;
                            HitClsInd = 1;
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_SPARSE_MULT */
                            HitClsInd = 0;
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb =  TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                    value = Part->inv[InvLab[k]];
                                    if (Markers[value] != tv->mark) {
                                        if (cls[value] > 1) HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                    }
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_SPARSE_MULT */
                            
                        }
                        
                        tv->mark++;
                        
                        SplInd = 0;
                        for (j = 0; j < HitClsInd; j++) {
                            ind1 = HitCls[j];
                            ind2 = ind1+cls[ind1];
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    break;
                                }
                            }
                        }
                        
                        /* DENSE-SPARSE CASE */
                        if (SplInd) {
                            SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+2*n, &Traceccend)
                            
                            /* Sorting the cells to be split */
                            sort_Split_Array(SplCls, SplInd);
                            
                            for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
                                ind0 = SplCls[j];
                                ind1 = ind0+cls[ind0];
                                SplCntInd = 0;
                                
                                /* According to the numbers of neighbors of C into the current cell */
                                /* compute how many vertices in C will be placed into the same new cell */
                                for (i = ind0; i < ind1; i++) {
                                    value = NghCounts[lab[i]];
                                    if (Markers[value] != tv->mark) {
                                        Markers[value] = tv->mark;
                                        SplCnt[SplCntInd++] = value;
                                        SplPos[value] = 1;
                                    }
                                    else {
                                        SplPos[value]++;
                                    }
                                }
                                tv->mark++;
                                
                                if (SplCntInd) {
                                    SAMETRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+2*n, Tracestpend)
                                }
                                
                                /* Sort the values deriving from the previous step */
                                sort_Split_Array(SplCnt, SplCntInd);
                                Part->cells += SplCntInd-1;
                                
                                /* Split the cell C and update the information for sizes of new cells */
                                /* Put the new cells into the stack */
                                i = ind0;
                                if (StackMarkers[i] != tv->stackmark) {
                                    BigCellSize = 0;
                                }
                                for (k = 0; k < SplCntInd; k++) {
                                    value = SplPos[SplCnt[k]];
                                    cls[i] = value;
                                    if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                        BigCell = i;
                                        BigCellPos = CStackInd;
                                        BigCellSize = cls[i];
                                    }
                                    SplPos[SplCnt[k]] = i;
                                    i += value;
                                    if (i < ind1) {
                                        CStack[++CStackInd] = i;
                                        StackMarkers[i] = tv->stackmark;
                                        SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
                                    }
                                }
                                
                                if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                    CStack[BigCellPos] = ind0;
                                    StackMarkers[BigCell] = 0;
                                    StackMarkers[ind0] = tv->stackmark;
                                }
                                
                                /* Permute elements of the cell C */
                                i = ind0;
                                do {
                                    SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                                }
                                while(++i < ind1);
                                
                                /* Reconstruct the cell C and update the inverse partition */
                                newcell = ind0;
                                i = ind0;
                                ind2 = newcell+cls[newcell]-1;
                                do {
                                    lab[i] = SplCnt[i];
                                    InvLab[lab[i]] = i;
                                    Part->inv[i] = newcell;
                                    if (i == ind2) {
                                        newcell = i+1;
                                        if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                    }
                                }
                                while (++i < ind1);
                                
                                for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
                                    if (cls[i] == 1) {
                                        Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                    }
                                }
                                
                            }
                        }
                        else {
                            if (TraceCell) {
                                return 0;
                            }
                        }
                        
                    }
                    else {
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_DENSE_MULT */
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb = TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_DENSE_MULT */
                            
                        }
                        SplInd = 0;
                        ind4 = 0;
                        while (ind4 < n) {	/* For each cell C with size(C) > 1 */
                            ind1 = ind4+cls[ind4];
                            if (cls[ind4] > 1) {
                                
                                /* Determine whether C must be split */
                                SplCntInd = 0;
                                value = NghCounts[lab[ind4]];
                                for (i = ind4+1; i < ind1; i++) {
                                    if (NghCounts[lab[i]] != value)
                                    {
                                        SplInd++;
                                        Markers[value] = tv->mark;
                                        SplCnt[SplCntInd++] = value;
                                        SplPos[value] = i-ind4;
                                        do {
                                            value = NghCounts[lab[i++]];
                                            if (Markers[value] != tv->mark) {
                                                Markers[value] = tv->mark;
                                                SplCnt[SplCntInd++] = value;
                                                SplPos[value] = 1;
                                            }
                                            else {
                                                SplPos[value]++;
                                            }
                                        }
                                        while(i != ind1);
                                        break;
                                    }
                                }
                                tv->mark++;
                                
                                if (SplCntInd) {
                                    SAMETRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+3*n, Tracestpend)
                                    
                                    /* Sort the values deriving from the previous step */
                                    sort_Split_Array(SplCnt, SplCntInd);
                                    
                                    Part->cells += SplCntInd-1;
                                    
                                    /* Split the cell C and update the information for sizes of new cells */
                                    /* Put the new cells into the stack */
                                    i = ind4;
                                    if (StackMarkers[i] != tv->stackmark) {
                                        BigCellSize = 0;
                                    }
                                    for (k = 0; k < SplCntInd; k++) {
                                        value = SplPos[SplCnt[k]];
                                        cls[i] = value;
                                        if ((StackMarkers[ind4] != tv->stackmark) && (value > BigCellSize)) {
                                            BigCell = i;
                                            BigCellPos = CStackInd;
                                            BigCellSize = cls[i];
                                        }
                                        SplPos[SplCnt[k]] = i;
                                        i += value;
                                        if (i < ind1) {
                                            CStack[++CStackInd] = i;
                                            StackMarkers[i] = tv->stackmark;
                                            SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
                                        }
                                    }
                                    if ((StackMarkers[ind4] != tv->stackmark) && (BigCell != ind4)) {
                                        CStack[BigCellPos] = ind4;
                                        StackMarkers[BigCell] = 0;
                                        StackMarkers[ind4] = tv->stackmark;
                                    }
                                    
                                    /* Permute elements of the cell C */
                                    i = ind4;
                                    do {
                                        SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                                    }
                                    while(++i < ind1);
                                    
                                    /* Reconstruct the cell C and update the inverse partition */
                                    newcell = ind4;
                                    i = ind4;
                                    ind2 = newcell+cls[newcell]-1;
                                    do {
                                        lab[i] = SplCnt[i];
                                        InvLab[lab[i]] = i;
                                        Part->inv[i] = newcell;
                                        if (i == ind2) {
                                            newcell = i+1;
                                            if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                        }
                                    }
                                    while (++i < ind1);
                                    
                                    for (i = ind4, k = 0; k < SplCntInd; i+=cls[i], k++) {
                                        if (cls[i] == 1) {
                                            Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                        }
                                    }
                                    
                                }
                            }
                            ind4 = ind1;
                        }
                        
                        /* DENSE-DENSE CASE */
                        if (SplInd) {
                            SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+3*n, &Traceccend)
                        }
                        else {
                            if (TraceCell) {
                                return 0;
                            }
                        }
                        
                    }
                }
            }
        }
        while (weightend < iend1int);
    }  /* end while (CStackInd > 0) */
    
    for (i=SpineTL->trcstart; i < TraceInd; i++) {
        ind0 = TheTrace[i];
        longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
        labi = lab[ind0];
        iend1int = TheGraph[labi].d;
        nghb = TheGraph[labi].e;
        for (j1int = 0; j1int < iend1int; ++j1int) {
            k = nghb[j1int];
            value = Part->inv[InvLab[k]];
            longcode = MASHCOMM(longcode, value);
        }
    }
    
    tv->augmented_cells = Part->cells - tv->augmented_cells;
    
    Part->code = Cand->code = CLEANUP(longcode);
    if ((Cand->code != SpineTL->part->code) || (TraceInd != *TraceEnd)) return FALSE;
    return TRUE;
}

void refine_tr(sparsegraph *sg, int *lab, int *ptn, int *numcells, int *code, TracesOptions *options_arg) {
    
    const int n = sg->nv;
    const int m = SETWORDSNEEDED(n);
    
    int i, j, ord;
    Partition *CurrPart;
    Candidate *CurrCand;
    struct TracesVars tvar, *tv;
    struct TracesInfo tinf, *ti;
    
    if (n > (NAUTY_INFINITY-2))
    {
        fprintf(ERRFILE, "Traces: need n <= %d, but n=%d\n\n",
                NAUTY_INFINITY-2, n);
        return;
    }
    
    Allocate_refine_Structures(n);
    
    tv = &tvar;
    ti = &tinf;
    
    tv->options = options_arg;
    tv->mark = tv->stackmark = NAUTY_INFINITY-1;
    tv->maxdeg = 0;
    tv->mindeg = n;
    
    outfile = (tv->options->outfile == NULL ? stdout : tv->options->outfile);
    
    for (i = 0; i < n; i++) {
        IDENTITY_PERM[i] = i;
    }
    
    copy_sg_structure(&redgraph, sg);
    
    tv->graph = &redgraph;
    if (tv->options->weighted) {
        tv->graph->w = malloc(tv->graph->wlen*sizeof(int));
        if (tv->graph->w == NULL) {
            fprintf(ERRFILE, "\nError, memory not allocated.\n");
            exit(1);
        }
        memcpy(tv->graph->w, sg->w, tv->graph->wlen*sizeof(int));
    }
    memcpy(tv->graph->e, sg->e, tv->graph->elen*sizeof(int));
    
    for (i=0; i<n; i++) {
        TheGraph[i].d = sg->d[i];
        if (TheGraph[i].d > tv->maxdeg) {
            tv->maxdeg = TheGraph[i].d;
        }
        if (TheGraph[i].d < tv->mindeg) {
            tv->mindeg = TheGraph[i].d;
        }
        TheGraph[i].e = tv->graph->e + sg->v[i];
        if (sg->w)
            TheGraph[i].w = tv->graph->w + sg->v[i];
        else
            TheGraph[i].w = NULL;
        TheGraph[i].one = FALSE;
    }
    
    ord = 0;
    
    /*----------- WEIGHTS --------------*/
    if (tv->options->weighted) {
        WeightCodes(n);
        ord = trie_classify(n,tv);
    }
    /*----------------------------------*/
    
    if ((tv->maxdeg == tv->mindeg) && (ord == 0)) ti->regular = TRUE; else ti->regular = FALSE;
    
    /* The graph is sparse? */
    if (sg->nde < n || sg->nde / n < n / (sg->nde / n)) {
        ti->thegraphisparse = TRUE;
    }
    else {
        ti->thegraphisparse = FALSE;
    }
    
    /* Initialize candidate, partition, cells, orbits */
    CurrCand = NewCandidate(n, &GarbList, TRUE);
    CurrPart = NewPartition(n);
    memset(CurrPart->inv, 0, n*sizeof(int));
    
    CurrCand->singcode = 0;
    
    if (ti->regular) {
        memcpy(CurrCand->lab, lab, n*sizeof(int));
        CurrPart->cells = 0;
        j = 0;
        for (i = 0; i < n; i++) {
            if (j) CurrPart->inv[i] = j;
            CurrCand->invlab[CurrCand->lab[i]] = i;
            if (!ptn[i]) {
                CurrPart->cls[j] = i-j+1;
                if (CurrPart->cls[j] == 1) {
                    CurrCand->singcode = MASHCOMM(CurrCand->singcode, CurrCand->lab[j]);
                }
                TheTrace[CurrPart->cells++] = j;
                j = i+1;
            }
        }
    } else {
        if (tv->options->weighted) {
            CurrPart->cells = traces_vertexclass_refine (n, lab, ptn, CurrCand, CurrPart, WeightsSeq);
        }
        else {
            CurrPart->cells = traces_vertexclass_refine (n, lab, ptn, CurrCand, CurrPart, sg->d);
        }
    }
    
    /* First refinement */
    refine_tr_refine(CurrCand, n, CurrPart, tv, ti);
    
    for (i = CurrPart->cls[0]; i < n; i+=CurrPart->cls[i]) {
        ptn[i-1] = 0;
    }
    ptn[n-1] = 0;
    
    memcpy(lab, CurrCand->lab, n*sizeof(int));
    *code = CurrCand->code;
    *numcells = CurrPart->cells;
    
    FREECAND(CurrCand)
    FREEPART(CurrPart)
    
    if (tv->graph != sg) {
        SG_FREE(redgraph);
    }
    
#if !MAXN
    DYNFREE(CStack, CStack_sz);
    DYNFREE(IDENTITY_PERM, IDENTITY_PERM_sz);
    DYNFREE(Markers, Markers_sz);
    DYNFREE(MarkHitVtx, MarkHitVtx_sz);
    DYNFREE(TreeMarkers, TreeMarkers_sz);
    DYNFREE(NghCounts, NghCounts_sz);
    DYNFREE(Singletons, Singletons_sz);
    DYNFREE(SplPos, SplPos_sz);
    DYNFREE(SplCls, SplCls_sz);
    DYNFREE(SplCnt, SplCnt_sz);
    DYNFREE(StackMarkers, StackMarkers_sz);
    DYNFREE(TheTrace, TheTrace_sz);
    DYNFREE(TheTraceSteps, TheTraceSteps_sz);
    DYNFREE(TheTraceCC, TheTraceCC_sz);
    DYNFREE(TheTraceSplNum, TheTraceSplNum_sz);
    DYNFREE(WeightsSeq, WeightsSeq_sz);
    DYNFREE(WorkArray1, WorkArray1_sz);
    DYNFREE(WorkArray2, WorkArray2_sz);
    DYNFREE(WorkArray3, WorkArray3_sz);
    DYNFREE(WorkArray4, WorkArray4_sz);
    DYNFREE(WorkArray5, WorkArray5_sz);
    DYNFREE(TreeStack, TreeStack_sz);
    DYNFREE(Spine, Spine_sz);
    DYNFREE(TrieArray, TrieArray_sz);
    DYNFREE(TheGraph, TheGraph_sz);
#endif
}

int refine_tr_refine(Candidate *Cand,
                     int n,
                     Partition *Part,
                     struct TracesVars* tv,
                     struct TracesInfo *ti) {
    
    int i, j, k, jk, sc, ind0, ind1, ind2, ind3, ind4, labi;
    int value, iend, newcell;
    int HitClsInd, SplInd, SplCntInd, CStackInd, TraceInd, TraceCCInd, TraceStepsInd, SingInd;
    int j1int, iend1int;
    unsigned int longcode;
    int newtrace = FALSE;
    int Sparse = TRUE;
    int *lab, *cls, *InvLab, *TracePos, *SplitCell, *LabCell, TraceEnd, Traceccend, Tracestpend;
    int BigCell, BigCellPos, BigCellSize;
    boolean TraceCell = FALSE;
    int *nghb;
    const int variation = 0;
    int currentweight, weightstart, weightend, currentcell, currentsize;
    
    HitClsInd = 0;
    if (tv->stackmark > (NAUTY_INFINITY-2)) {
        memset(StackMarkers, 0, n*sizeof(int));
        tv->stackmark = 0;
    }
    tv->stackmark++;
    
    lab = Cand->lab;
    InvLab = Cand->invlab;
    cls = Part->cls;
    
    UPDATEMIN(Part->active, n-1);
    memcpy(CStack+1, TheTrace, (Part->cells)*sizeof(int));
    CStackInd = Part->cells;
    
    for (i = 1; i <= CStackInd; i++) {
        StackMarkers[CStack[i]] = tv->stackmark;
    }
    
    longcode = Part->cells;
    newtrace = TRUE;
    while (CStackInd > 0) {
        
        weightend = 0;
        
        if (tv->mark > (NAUTY_INFINITY-2)) {
            memset(Markers, 0, n*sizeof(int));
            memset(MarkHitVtx, 0, n*sizeof(int));
            tv->mark = 0;
        }
        tv->mark++;
        
        if (Part->cells == n) break;
        
        k = Select_from_CStack(cls, CStackInd);
        
        currentcell = CStack[k];
        currentsize = currentcell+cls[currentcell];
        CStack[k] = CStack[CStackInd--];
        StackMarkers[currentcell] = 0;
        
        labi = lab[currentcell];
        iend1int = TheGraph[labi].d;
        nghb =  TheGraph[labi].e;
        
        do {
            
            ind0 = currentcell;
            ind2 = currentsize;
            weightstart = weightend;
            
            if (tv->options->weighted) {
                currentweight = (TheGraph[labi].w)[weightstart];
                while ((iend1int > weightend) && ((TheGraph[labi].w)[weightend] == currentweight)) {
                    weightend++;
                }
            } else {
                weightend = TheGraph[labi].d;
            }
            
            /* Analysis of occurrences of neighbors of the current cell */
            /* The list of cells with neighbors in the current cell is  built */
            if (cls[ind0] == 1) {			/* SINGLETON CURRENT CELL CASE */
                
                /* NEIGHCOUNT_SING_MULT */
                HitClsInd = 0;
                for (j1int = weightstart; j1int < weightend; ++j1int) {
                    k = nghb[j1int];
                    value = Part->inv[InvLab[k]];
                    if (cls[value] > 1) {
                        if (Markers[value] != tv->mark) {
                            HitCls[HitClsInd++] = value;
                            Markers[value] = tv->mark;
                            ElmHitCll[value] = value;
                        }
                        HitVtx[ElmHitCll[value]++] = k;
                    } else {
                        switch (variation) {
                            case 1:
                                longcode = MASHCOMM(longcode, value);
                                break;
                            default:
                                break;
                        }
                    }
                }
                /* end NEIGHCOUNT_SING_MULT */
                
                tv->mark++;
                FIND_SPLIT_CELLS;
                
                /* SINGLETON CC CASE */
                if (SplInd) {
                    
                    /* Sorting the cells to be split */
                    sort_Split_Array(SplCls, SplInd);
                    
                    for (j = 0; j < SplInd; j++) {
                        ind1 = SplCls[j];
                        i = ind1+cls[ind1]-ElmHitCll[ind1];
                    }
                    
                    /* REARRANGE THE CELLS */
                    for (j = 0; j < SplInd; j++) {
                        /* REARRANGE_CELLS */
                        ind1 = SplCls[j];
                        cls[ind1] = cls[ind1]-ElmHitCll[ind1];
                        newcell = ind1+cls[ind1];
                        cls[newcell] = ElmHitCll[ind1];
                        Part->cells++;
                        if (StackMarkers[ind1] != tv->stackmark) {
                            if (cls[newcell] < cls[ind1]) {
                                CStack[++CStackInd] = newcell;
                                StackMarkers[newcell] = tv->stackmark;
                            }
                            else {
                                CStack[++CStackInd] = ind1;
                                StackMarkers[ind1] = tv->stackmark;
                            }
                        }
                        else {
                            CStack[++CStackInd] = newcell;
                            StackMarkers[newcell] = tv->stackmark;
                        }
                        SplitCell = HitVtx+ind1;
                        ind3 = cls[newcell];
                        LabCell = lab+newcell;
                        for (jk = 0; jk < ind3; jk++) {
                            k = SplitCell[jk];
                            i = LabCell[jk];
                            Part->inv[newcell+jk] = newcell;
                            lab[InvLab[k]] = i;
                            InvLab[i] = InvLab[k];
                            LabCell[jk] = k;
                            InvLab[k] = newcell+jk;
                        }
                        /* END REARRANGE_CELLS */
                        
                    }
                }
                else {
                }
            }
            else {
                if (ti->thegraphisparse) {
                    
                    /* NEIGHCOUNT_SPARSE_MULT */
                    if (cls[ind0] != n) {
                        HitClsInd = 0;
                        for (i = ind0; i < ind2; i++) {
                            labi = lab[i];
                            nghb = TheGraph[labi].e;
                            for (j = weightstart; j < weightend; j++) {
                                k = nghb[j];
                                if (MarkHitVtx[k] == tv->mark) {
                                    NghCounts[k]++;
                                }
                                else {
                                    value = Part->inv[InvLab[k]];
                                    if (cls[value] > 1) {
                                        MarkHitVtx[k] = tv->mark;
                                        NghCounts[k] = 1;
                                        if (Markers[value] != tv->mark) {
                                            HitCls[HitClsInd++] = value;
                                            Markers[value] = tv->mark;
                                            HitVtx[value] = k;
                                            ElmHitCll[value] = 1;
                                        }
                                        else {
                                            HitVtx[value+ElmHitCll[value]++] = k;
                                        }
                                    }
                                    else {
                                        switch (variation) {
                                            case 1:
                                                longcode = MASHCOMM(longcode, value);
                                                break;
                                            default:
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /* End NEIGHCOUNT_SPARSE_MULT */
                    
                    tv->mark++;
                    
                    SplInd = 0;
                    SplCls[0] = n;
                    for (j = 0; j < HitClsInd; j++) {
                        ind1 = HitCls[j];
                        if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
                            SplCls[SplInd++] = ind1;
                        }
                        else {
                            ind2 = ind1+cls[ind1];
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    break;
                                }
                            }
                        }
                    }
                    
                    /* SPARSE CASE */
                    if (SplInd) {
                        
                        /* Sorting the cells to be split */
                        sort_Split_Array(SplCls, SplInd);
                        
                        for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
                            ind0 = SplCls[sc];
                            ind1 = ind0 + cls[ind0];
                            SplCntInd = 0;
                            if (ElmHitCll[ind0] < cls[ind0]) {
                                SplCnt[SplCntInd++] = 0;
                                SplPos[0] = cls[ind0] - ElmHitCll[ind0];
                            }
                            
                            /* According to the numbers of neighbors of C into the current cell */
                            /* compute how many vertices in C will be placed into the same new cell */
                            iend = ind0 + ElmHitCll[ind0];
                            for (i = ind0; i < iend; i++) {
                                value = NghCounts[HitVtx[i]];
                                if (Markers[value] != tv->mark) {
                                    Markers[value] = tv->mark;
                                    SplCnt[SplCntInd++] = value;
                                    SplPos[value] = 1;
                                }
                                else {
                                    SplPos[value]++;
                                }
                            }
                            tv->mark++;
                            
                            sort_Split_Array(SplCnt,SplCntInd);
                            Part->cells += SplCntInd-1;
                            
                            /* Split the cell C and update the information for sizes of new cells */
                            /* Put the new cells into the stack */
                            i = ind0;
                            if (StackMarkers[i] != tv->stackmark) {
                                BigCellSize = 0;
                            }
                            for (k = 0; k < SplCntInd; k++) {
                                value = SplPos[SplCnt[k]];
                                cls[i] = value;
                                if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                    BigCell = i;
                                    BigCellPos = CStackInd;
                                    BigCellSize = cls[i];
                                }
                                SplPos[SplCnt[k]] = i;
                                i += value;
                                if (i < ind1) {
                                    CStack[++CStackInd] = i;
                                    StackMarkers[i] = tv->stackmark;
                                }
                            }
                            
                            if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                CStack[BigCellPos] = ind0;
                                StackMarkers[BigCell] = 0;
                                StackMarkers[ind0] = tv->stackmark;
                            }
                            /* Permute elements of the cell C */
                            iend = ind0 + ElmHitCll[ind0];
                            for (i = ind0; i < iend; i++) {
                                value = HitVtx[i];
                                j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
                                k = InvLab[value];				/* where HitVtx[i] is in lab */
                                lab[k] = lab[j];
                                lab[j] = value;
                                InvLab[value] = j;
                                InvLab[lab[k]] = k;
                                NghCounts[value] = 0;
                            }
                            
                            /* Reconstruct the cell C and update the inverse partition */
                            newcell = ind1 - ElmHitCll[ind0];
                            i = newcell;
                            ind2 = newcell+cls[newcell]-1;
                            do {
                                Part->inv[i] = newcell;
                                if (i == ind2) {
                                    newcell = i+1;
                                    if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                }
                            }
                            while (++i < ind1);
                        }
                    }
                }
                else {
                    if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
                        Sparse = FALSE;
                    }
                    else {
                        Sparse = TRUE;
                    }
                    if (Sparse) {
                        /* Counting occurrences of neighbors of the current cell */
                        /* The list of cells with neighbors in the current cell is also built */
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                            HitCls[0] = 0;
                            HitClsInd = 1;
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_SPARSE_MULT */
                            HitClsInd = 0;
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb =  TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                    value = Part->inv[InvLab[k]];
                                    if (Markers[value] != tv->mark) {
                                        if (cls[value] > 1) HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                    }
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_SPARSE_MULT */
                            
                        }
                        
                        tv->mark++;
                        
                        
                        SplInd = 0;
                        for (j = 0; j < HitClsInd; j++) {
                            ind1 = HitCls[j];
                            ind2 = ind1+cls[ind1];
                            value = NghCounts[lab[ind1++]];
                            for (i = ind1; i < ind2; i++)
                            {
                                if (NghCounts[lab[i]] != value)
                                {
                                    SplCls[SplInd++] = HitCls[j];
                                    break;
                                }
                            }
                        }
                        
                        /* DENSE-SPARSE CASE */
                        if (SplInd) {
                            
                            /* Sorting the cells to be split */
                            sort_Split_Array(SplCls, SplInd);
                            
                            for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
                                ind0 = SplCls[j];
                                ind1 = ind0+cls[ind0];
                                SplCntInd = 0;
                                
                                /* According to the numbers of neighbors of C into the current cell */
                                /* compute how many vertices in C will be placed into the same new cell */
                                for (i = ind0; i < ind1; i++) {
                                    value = NghCounts[lab[i]];
                                    if (Markers[value] != tv->mark) {
                                        Markers[value] = tv->mark;
                                        SplCnt[SplCntInd++] = value;
                                        SplPos[value] = 1;
                                    }
                                    else {
                                        SplPos[value]++;
                                    }
                                }
                                tv->mark++;
                                
                                /* Sort the values deriving from the previous step */
                                sort_Split_Array(SplCnt, SplCntInd);
                                Part->cells += SplCntInd-1;
                                
                                /* Split the cell C and update the information for sizes of new cells */
                                /* Put the new cells into the stack */
                                i = ind0;
                                if (StackMarkers[i] != tv->stackmark) {
                                    BigCellSize = 0;
                                }
                                for (k = 0; k < SplCntInd; k++) {
                                    value = SplPos[SplCnt[k]];
                                    cls[i] = value;
                                    
                                    if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
                                        BigCell = i;
                                        BigCellPos = CStackInd;
                                        BigCellSize = cls[i];
                                    }
                                    SplPos[SplCnt[k]] = i;
                                    i += value;
                                    if (i < ind1) {
                                        CStack[++CStackInd] = i;
                                        StackMarkers[i] = tv->stackmark;
                                    }
                                }
                                
                                if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
                                    CStack[BigCellPos] = ind0;
                                    StackMarkers[BigCell] = 0;
                                    StackMarkers[ind0] = tv->stackmark;
                                }
                                
                                /* Permute elements of the cell C */
                                i = ind0;
                                do {
                                    SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                                }
                                while(++i < ind1);
                                
                                /* Reconstruct the cell C and update the inverse partition */
                                newcell = ind0;
                                i = ind0;
                                ind2 = newcell+cls[newcell]-1;
                                do {
                                    lab[i] = SplCnt[i];
                                    InvLab[lab[i]] = i;
                                    
                                    Part->inv[i] = newcell;
                                    if (i == ind2) {
                                        newcell = i+1;
                                        if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                    }
                                }
                                while (++i < ind1);
                            }
                        }
                    }
                    else {
                        if (cls[ind0] == n) {
                            for (i = 0; i < n; i++) {
                                NghCounts[i] = TheGraph[i].d;
                            }
                        }
                        else {
                            memset(NghCounts, 0, n*sizeof(int));
                            
                            /* NEIGHCOUNT_DENSE_DENSE_MULT */
                            for (i = ind0; i < ind2; i++) {
                                labi = lab[i];
                                nghb = TheGraph[labi].e;
                                for (j1int = weightstart; j1int < weightend; ++j1int) {
                                    k = nghb[j1int];
                                    (NghCounts[k])++;
                                }
                            }
                            /* End NEIGHCOUNT_DENSE_DENSE_MULT */
                            
                        }
                        SplInd = 0;
                        ind4 = 0;
                        while (ind4 < n) {	/* For each cell C with size(C) > 1 */
                            ind1 = ind4+cls[ind4];
                            if (cls[ind4] > 1) {
                                
                                
                                /* Determine whether C must be split */
                                SplCntInd = 0;
                                value = NghCounts[lab[ind4]];
                                for (i = ind4+1; i < ind1; i++) {
                                    if (NghCounts[lab[i]] != value)
                                    {
                                        SplInd++;
                                        Markers[value] = tv->mark;
                                        SplCnt[SplCntInd++] = value;
                                        SplPos[value] = i-ind4;
                                        do {
                                            value = NghCounts[lab[i++]];
                                            if (Markers[value] != tv->mark) {
                                                Markers[value] = tv->mark;
                                                SplCnt[SplCntInd++] = value;
                                                SplPos[value] = 1;
                                            }
                                            else {
                                                SplPos[value]++;
                                            }
                                        }
                                        while(i != ind1);
                                        break;
                                    }
                                }
                                tv->mark++;
                                
                                if (SplCntInd) {
                                    
                                    /* Sort the values deriving from the previous step */
                                    sort_Split_Array(SplCnt, SplCntInd);
                                    Part->cells += SplCntInd-1;
                                    
                                    /* Split the cell C and update the information for sizes of new cells */
                                    /* Put the new cells into the stack */
                                    i = ind4;
                                    if (StackMarkers[i] != tv->stackmark) {
                                        BigCellSize = 0;
                                    }
                                    for (k = 0; k < SplCntInd; k++) {
                                        value = SplPos[SplCnt[k]];
                                        cls[i] = value;
                                        
                                        if ((StackMarkers[ind4] != tv->stackmark) && (value > BigCellSize)) {
                                            BigCell = i;
                                            BigCellPos = CStackInd;
                                            BigCellSize = cls[i];
                                        }
                                        SplPos[SplCnt[k]] = i;
                                        i += value;
                                        if (i < ind1) {
                                            CStack[++CStackInd] = i;
                                            StackMarkers[i] = tv->stackmark;
                                        }
                                    }
                                    if ((StackMarkers[ind4] != tv->stackmark) && (BigCell != ind4)) {
                                        CStack[BigCellPos] = ind4;
                                        StackMarkers[BigCell] = 0;
                                        StackMarkers[ind4] = tv->stackmark;
                                    }
                                    
                                    /* Permute elements of the cell C */
                                    i = ind4;
                                    do {
                                        SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
                                    }
                                    while(++i < ind1);
                                    
                                    /* Reconstruct the cell C and update the inverse partition */
                                    newcell = ind4;
                                    i = ind4;
                                    ind2 = newcell+cls[newcell]-1;
                                    do {
                                        lab[i] = SplCnt[i];
                                        InvLab[lab[i]] = i;
                                        Part->inv[i] = newcell;
                                        if (i == ind2) {
                                            newcell = i+1;
                                            if (newcell < n) ind2 = newcell+cls[newcell]-1;
                                        }
                                    }
                                    while (++i < ind1);
                                }
                            }
                            ind4 = ind1;
                        }
                        
                        /* DENSE-DENSE CASE */
                    }
                }
            }
        }
        while (weightend < iend1int);
    }  /* end while (CStackInd > 0) */
    
    Cand->code = CLEANUP(longcode);
    for (ind0=Part->cls[0]; ind0 < n; ind0+=Part->cls[ind0]) {
        longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
        labi = lab[ind0];
        iend1int = TheGraph[labi].d;
        nghb = TheGraph[labi].e;
        for (j1int = 0; j1int < iend1int; ++j1int) {
            k = nghb[j1int];
            value = Part->inv[InvLab[k]];
            longcode = MASHCOMM(longcode, value);
        }
    }
    Part->code = Cand->code = CLEANUP(longcode);
    return TraceInd;
}

void Allocate_Traces_Structures(int n) {
#if !MAXN
    DYNALLOC1(int, AUTPERM, AUTPERM_sz, n, "Traces");
    DYNALLOC1(int, BreakSteps, BreakSteps_sz, n, "Traces");
    DYNALLOC1(int, CurrOrbSize, CurrOrbSize_sz, n, "Traces");
    DYNALLOC1(int, CurrRefCells, CurrRefCells_sz, n, "Traces");
    DYNALLOC1(boolean, Diff, Diff_sz, n, "Traces");
    DYNALLOC1(int, CStack, CStack_sz, n, "Traces");
    DYNALLOC1(boolean, Factorials, Factorials_sz, n, "Traces");
    DYNALLOC1(int, fix, fix_sz, n, "Traces");
    DYNALLOC1(int, IDENTITY_PERM, IDENTITY_PERM_sz, n, "Traces");
    DYNALLOC1(int, Markers, Markers_sz, n, "Traces");
    DYNALLOC1(int, TreeMarkers, TreeMarkers_sz, n, "Traces");
    DYNALLOC1(int, AutMarkers, AutMarkers_sz, n, "Traces");
    DYNALLOC1(int, MarkHitVtx, MarkHitVtx_sz, n, "Traces");
    DYNALLOC1(int, MultRefCells, MultRefCells_sz, n, "Traces");
    DYNALLOC1(int, NghCounts, NghCounts_sz, n, "Traces");
    DYNALLOC1(int, OrbSize, OrbSize_sz, n, "Traces");
    DYNALLOC1(int, OrbList, OrbList_sz, n, "Traces");
    DYNALLOC1(pair, PrmPairs, PrmPairs_sz, n, "Traces");
    DYNALLOC1(int, TempOrbList, TempOrbList_sz, n, "Traces");
    DYNALLOC1(int, RefCells, RefCells_sz, n, "Traces");
    DYNALLOC1(int, Singletons, Singletons_sz, n, "Traces");
    DYNALLOC1(int, SplCls, SplCls_sz, n, "Traces");
    DYNALLOC1(int, SplCnt, SplCnt_sz, n, "Traces");
    DYNALLOC1(int, SplPos, SplPos_sz, n, "Traces");
    DYNALLOC1(int, StackMarkers, StackMarkers_sz, n, "Traces");
    DYNALLOC1(int, TheTrace, TheTrace_sz, n+10, "Traces");
    DYNALLOC1(int, TheTraceCC, TheTraceCC_sz, n, "Traces");
    DYNALLOC1(int, TheTraceSplNum, TheTraceSplNum_sz, n, "Traces");
    DYNALLOC1(int, TheTraceSteps, TheTraceSteps_sz, n+10, "Traces");
    DYNALLOC1(int, TEMPLAB, TEMPLAB_sz, n, "Traces");
    DYNALLOC1(int, TEMPINVLAB, TEMPINVLAB_sz, n, "Traces");
    DYNALLOC1(int, WeightsSeq, WeightsSeq_sz, n, "Traces");
    DYNALLOC1(int, WorkArray, WorkArray_sz, n, "Traces");
    DYNALLOC1(int, WorkArray0, WorkArray0_sz, n, "Traces");
    DYNALLOC1(int, WorkArray1, WorkArray1_sz, n, "Traces");
    DYNALLOC1(int, WorkArray2, WorkArray2_sz, n, "Traces");
    DYNALLOC1(int, WorkArray3, WorkArray3_sz, n, "Traces");
    DYNALLOC1(int, WorkArray4, WorkArray4_sz, n, "Traces");
    DYNALLOC1(int, WorkArray5, WorkArray5_sz, n, "Traces");
    DYNALLOC1(int, WorkArray6, WorkArray6_sz, n, "Traces");
    DYNALLOC1(int, WorkArray7, WorkArray7_sz, n, "Traces");
    DYNALLOC1(int, TreeStack, TreeStack_sz, n, "Traces");
    DYNALLOC1(TracesSpine, Spine, Spine_sz, n, "Traces");
    DYNALLOC1(trie*, TrieArray, TrieArray_sz, n, "Traces");
    DYNALLOC1(grph_strct, TheGraph, TheGraph_sz, n, "Traces");
    DYNALLOC1(ExpPathInfo, EPCodes, EPCodes_sz, n, "Traces");
#endif
    return;
}

void Allocate_refine_Structures(int n) {
#if !MAXN
    DYNALLOC1(int, CStack, CStack_sz, n, "refine_tr");
    DYNALLOC1(int, IDENTITY_PERM, IDENTITY_PERM_sz, n, "refine_tr");
    DYNALLOC1(int, Markers, Markers_sz, n, "refine_tr");
    DYNALLOC1(int, MarkHitVtx, MarkHitVtx_sz, n, "refine_tr");
    DYNALLOC1(int, TreeMarkers, TreeMarkers_sz, n, "refine_tr");
    DYNALLOC1(int, NghCounts, NghCounts_sz, n, "refine_tr");
    DYNALLOC1(int, Singletons, Singletons_sz, n, "refine_tr");
    DYNALLOC1(int, SplPos, SplPos_sz, n, "refine_tr");
    DYNALLOC1(int, SplCls, SplCls_sz, n, "refine_tr");
    DYNALLOC1(int, SplCnt, SplCnt_sz, n, "refine_tr");
    DYNALLOC1(int, StackMarkers, StackMarkers_sz, n, "refine_tr");
    DYNALLOC1(int, TheTrace, TheTrace_sz, n+10, "refine_tr");
    DYNALLOC1(int, TheTraceSteps, TheTraceSteps_sz, n+10, "refine_tr");
    DYNALLOC1(int, TheTraceCC, TheTraceCC_sz, n, "refine_tr");
    DYNALLOC1(int, TheTraceSplNum, TheTraceSplNum_sz, n, "refine_tr");
    DYNALLOC1(int, WeightsSeq, WeightsSeq_sz, n, "refine_tr");
    DYNALLOC1(int, WorkArray1, WorkArray1_sz, n, "refine_tr");
    DYNALLOC1(int, WorkArray2, WorkArray2_sz, n, "refine_tr");
    DYNALLOC1(int, WorkArray3, WorkArray3_sz, n, "refine_tr");
    DYNALLOC1(int, WorkArray4, WorkArray4_sz, n, "refine_tr");
    DYNALLOC1(int, WorkArray5, WorkArray5_sz, n, "refine_tr");
    DYNALLOC1(int, TreeStack, TreeStack_sz, n, "refine_tr");
    DYNALLOC1(TracesSpine, Spine, Spine_sz, n, "refine_tr");
    DYNALLOC1(trie*, TrieArray, TrieArray_sz, n, "refine_tr");
    DYNALLOC1(grph_strct, TheGraph, TheGraph_sz, n, "refine_tr");
#endif
    
#define HitCls WorkArray2
#define HitVtx WorkArray3
#define ElmHitCll WorkArray5
    return;
}

struct Candidate *NewCandidate(int n, Candidate **GarbList, int Mrk) {
    struct Candidate *Cand;
    
    if (*GarbList) {
        Cand = *GarbList;
        *GarbList = (*GarbList)->next;
    }
    else {
        Cand = malloc(sizeof(*Cand));
        if (Cand == NULL) {
            fprintf(ERRFILE, "\nError, memory not allocated.\n");
            exit(1);
        }
        Cand->lab = malloc(n*sizeof(*Cand->lab));
        if (Cand->lab == NULL) {
            fprintf(ERRFILE, "\nError, memory not allocated.\n");
            exit(1);
        }
        Cand->invlab = malloc(n*sizeof(*Cand->invlab));
        if (Cand->invlab == NULL) {
            fprintf(ERRFILE, "\nError, memory not allocated.\n");
            exit(1);
        }
    }
    Cand->do_it = Mrk;
    Cand->indnum = 0;
    Cand->code = 0;
    Cand->next = NULL;
    Cand->stnode = NULL;
    Cand->sortedlab = FALSE;
    return Cand;
}

int Check_degree_one(sparsegraph *sg, Candidate *Cand, Partition *Part, int n) {
    
    int i;
    
    for (i=0; i<n; i += Part->cls[i]) {
        if (sg->d[Cand->lab[i]] == 1) {
            return TRUE;
        }
    }
    return FALSE;
}

int CheckForAutomorphisms(Candidate *CurrCand, Candidate *NextCand,
                          struct TracesVars* tv, struct TracesInfo* ti,
                          int m, int n, Partition* Part) {
    Candidate *CheckAutList;
    int i, j, k, tgt_level, numtemporbits;
    int CheckLevel, CheckLevelEnd;
    int temp, tmp, tmp1, arg, arg1, val, val1;
    searchtrie *TrieCandFrom, *TrieCheckFrom;
    
    CheckLevel = 0;
    CheckLevelEnd = 0;
    temp = 0;
    tv->gotonode = NULL;
    tv->conta6++;
    
    switch (tv->compstage) {
        case 0:
            if (tv->strategy) {
                CheckLevel = CheckLevelEnd = tv->maxtreelevel;
            }
            else {
                if ((Spine[tv->tolevel].part)->cells == tv->finalnumcells) {
                    CheckLevel = CheckLevelEnd = tv->maxtreelevel;
                }
                else {
                    CheckLevel = 1;
                    if ((Spine[tv->maxtreelevel].part)->cells == tv->finalnumcells) {
                        CheckLevelEnd = tv->maxtreelevel - 1;
                    }
                    else {
                        CheckLevelEnd = tv->maxtreelevel;
                    }
                }
            }
            break;
        case 1:
            CheckLevel = CheckLevelEnd = tv->tolevel;
            break;
        case 2:
            if (m || (tv->tolevel == tv->maxtreelevel+1)) {
                CheckLevel = CheckLevelEnd = tv->maxtreelevel+1;
            }
            else {
                CheckLevel = 1;
                if ((Spine[tv->maxtreelevel].part)->cells == tv->finalnumcells) {
                    CheckLevelEnd = tv->maxtreelevel - 1;
                }
                else {
                    CheckLevelEnd = tv->maxtreelevel;
                }
            }
            break;
        default:
            break;
    }
    
    while (CheckLevel <= CheckLevelEnd) {
        CheckAutList = Spine[CheckLevel].liststart;
        while (CheckAutList) {
            if (CheckAutList->do_it && lookup(CheckAutList->stnode) && (CheckAutList != NextCand) && (CheckAutList != CurrCand)) {
                if (CheckAutList->code == NextCand->code) {
                    SETMARK(Markers, tv->mark)
                    if (Part->cells == n) {
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        for (i = 0; i < n; i++) {
                            arg = NextCand->lab[i];
                            val = CheckAutList->lab[i];
                            SETPAIRSAUT(arg, val)
                        }
                    }
                    else {
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        SETMARK(CellMarkers1, tv->markcell1)
                        SETMARK(CellMarkers2, tv->markcell2)
                        for (i=0; i<n; i+=Part->cls[i]) {
                            if (Part->cls[i] == 1) {
                                arg = NextCand->lab[i];
                                val = CheckAutList->lab[i];
                                SETPAIRSAUT(arg, val)
                                if (tv->preprocessed && Diff[arg])
                                    MakeTree(arg, val, tv->input_graph, n, tv, TRUE);
                            }
                            else {
                                k = i;
                                for (j=i; j<i+Part->cls[i]; j++) {
                                    arg = arg1 = NextCand->lab[j];
                                    if (CellMarkers1[arg] != tv->markcell1) {
                                        CellMarkers1[arg] = tv->markcell1;
                                        while ((CellMarkers2[CheckAutList->lab[k]] == tv->markcell2) && (k < i+Part->cls[i])) {
                                            k++;
                                        }
                                        if (k < i+Part->cls[i]) {
                                            val = val1 = CheckAutList->lab[k];
                                            CellMarkers2[val] = tv->markcell2;
                                            SETPAIRSAUT(arg, val)
                                            if (tv->preprocessed && Diff[arg])
                                                MakeTree(arg, val, tv->input_graph, n, tv, TRUE);
                                            tmp = FirstNeighbour(arg, NextCand, Part, CellMarkers1, tv->markcell1, &arg, n);
                                            if (tmp) {
                                                CellMarkers1[arg] = tv->markcell1;
                                                tmp = FirstNeighbour(val, CheckAutList, Part, CellMarkers2, tv->markcell2, &val, n);
                                                CellMarkers2[val] = tv->markcell2;
                                                SETPAIRSAUT(arg, val)
                                                if (tv->preprocessed && Diff[arg])
                                                    MakeTree(arg, val, tv->input_graph, n, tv, TRUE);
                                                while (tmp) {
                                                    tmp = NextNeighbour(arg, NextCand, Part, CellMarkers1, tv->markcell1, &arg, n);
                                                    if (tmp) {
                                                        CellMarkers1[arg] = tv->markcell1;
                                                        tmp = NextNeighbour(val, CheckAutList, Part, CellMarkers2, tv->markcell2, &val, n);
                                                        if (tmp) {
                                                            CellMarkers2[val] = tv->markcell2;
                                                            SETPAIRSAUT(arg, val)
                                                            if (tv->preprocessed && Diff[arg])
                                                                MakeTree(arg, val, tv->input_graph, n, tv, TRUE);
                                                        }
                                                    }
                                                }
                                                arg = arg1;
                                                val = val1;
                                                do {
                                                    tmp = NextNeighbour(arg, NextCand, Part, CellMarkers1, tv->markcell1, &arg, n);
                                                    if (tmp) {
                                                        CellMarkers1[arg] = tv->markcell1;
                                                        tmp = NextNeighbour(val, CheckAutList, Part, CellMarkers2, tv->markcell2, &val, n);
                                                        if (tmp) {
                                                            CellMarkers2[val] = tv->markcell2;
                                                            SETPAIRSAUT(arg, val)
                                                            if (tv->preprocessed && Diff[arg])
                                                                MakeTree(arg, val, tv->input_graph, n, tv, TRUE);
                                                        }
                                                    }
                                                } while (tmp);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    if (isautom_sg_pair((graph*)tv->input_graph, AUTPERM, tv->options->digraph, m, n, tv)) {
                        if (!findperm(gensB, AUTPERM, n)) {
                            if (tv->options->verbosity >= 2) tv->schreier3 -= CPUTIME;
                            addgenerator(&gpB, &gensB, AUTPERM, n);
                            if (tv->options->verbosity >= 2) tv->schreier3 += CPUTIME;
                            if (tv->options->verbosity >= 2) {
                                fprintf(outfile, "[A(%d,%d)] ", CheckLevel, CheckAutList->name);
                            }
                            tv->lev_of_lastauto = tv->tolevel;
                            tv->stats->numgenerators++;
                            orbjoin_sp_perm(tv->orbits, AUTPERM, OrbList, n, &tv->stats->numorbits);
                            ti->thegrouphaschanged = TRUE;
                            ti->identitygroup = FALSE;
                            if (tv->options->verbosity >= 2 && tv->options->writeautoms) {
                                PRINT_RETURN
                            }
                            if (tv->options->writeautoms) {
                                fprintf(outfile, "Gen(A) #%d: ", tv->stats->numgenerators);
                                writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n);
                            }
                            if (tv->options->userautomproc) {
                                (*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n);
                            }
                        }
                        else {
                            if (tv->options->verbosity >= 2) {
                                fprintf(outfile, "[A*(%d,%d)] ", CheckLevel, CheckAutList->name);
                            }
                        }
                        TrieCandFrom = NULL;
                        TrieCheckFrom = CheckAutList->stnode;
                        if (CurrCand->stnode->level <= 1) {
                            tgt_level = CurrCand->stnode->level + 1;
                            while (TrieCheckFrom->level > tgt_level) {
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if (tv->tolevel <= TrieCheckFrom->level) {
                                tgt_level = tv->tolevel;
                                while (TrieCheckFrom->level != tgt_level) {
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            }
                            else {
                                TrieCandFrom = CurrCand->stnode;
                                tgt_level = TrieCheckFrom->level;
                                while (TrieCandFrom->level != tgt_level) {
                                    TrieCandFrom = TrieCandFrom->father;
                                }
                            }
                        }
                        if (TrieCandFrom) {
                            while (TrieCandFrom->father != TrieCheckFrom->father) {
                                TrieCandFrom = TrieCandFrom->father;
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if ((TrieCheckFrom->level > 1) && (TrieCheckFrom->father != CurrCand->stnode)) {
                                TrieCandFrom = CurrCand->stnode;
                                TrieCheckFrom = TrieCheckFrom->father;
                                while (TrieCandFrom->father != TrieCheckFrom->father) {
                                    TrieCandFrom = TrieCandFrom->father;
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            }
                        }
                        
                        while (TrieCheckFrom->goes_to) {
                            TrieCheckFrom = TrieCheckFrom->goes_to;
                        }
                        
                        for (temp=1; temp<=tv->tolevel; temp++) {
                            if (CheckAutList->lab[Spine[temp].tgtpos] != NextCand->lab[Spine[temp].tgtpos]) {
                                break;
                            }
                        }
                        
                        if (temp == tv->tolevel) {
                            if (TempOrbits) {
                                if (tv->compstage == 0) {
                                    for (j=0; j<tv->permInd; j++) {
                                        orbjoin_sp_pair(TempOrbits, TempOrbList, n,
                                                        PrmPairs[j].arg, PrmPairs[j].val, &numtemporbits);
                                    }
                                }
                                else {
                                    orbjoin(TempOrbits, AUTPERM, n);
                                }
                            } else {
                                orbjoin(tv->currorbit, AUTPERM, n);
                            }
                        }
                        
                        switch (tv->compstage) {
                            case 0:
                                if (tv->strategy && (tv->steps == 1)) {
                                    RemoveFromLevel(temp, tv->maxtreelevel-1, tv->strategy, FALSE);
                                    if (TrieCandFrom) {
                                        TrieCheckFrom->index += TrieCandFrom->index;
                                        TrieCandFrom->goes_to = TrieCheckFrom;
                                        PRINT_INDEX(TrieCheckFrom,4,1)
                                    }
                                    else {
                                        TrieCheckFrom->index++;
                                        PRINT_INDEX(TrieCheckFrom,4,2)
                                    }
                                    NextCand->do_it = FALSE;
                                }
                                else {
                                    if (CheckAutList->lab[Spine[temp].tgtpos] >= NextCand->lab[Spine[temp].tgtpos]) {
                                        CheckAutList->do_it = FALSE;
                                        if (TrieCandFrom) {
                                            TrieCandFrom->index += TrieCheckFrom->index;
                                            tv->newindex = 0;
                                            TrieCheckFrom->goes_to = TrieCandFrom;
                                            PRINT_INDEX(TrieCandFrom,4,3)
                                        }
                                        else {
                                            if (CurrCand->stnode->level > 1) {
                                                tv->newgotonode = TrieCheckFrom;
                                                tv->newindex = TrieCheckFrom->index;
                                            }
                                            else {
                                                tv->newgotonode = NULL;
                                                tv->newindex = 0;
                                            }
                                        }
                                    }
                                    else {
                                        if (TrieCandFrom) {
                                            TrieCheckFrom->index += TrieCandFrom->index;
                                            TrieCandFrom->goes_to = TrieCheckFrom;
                                            PRINT_INDEX(TrieCheckFrom,4,4)
                                        }
                                        else {
                                            TrieCheckFrom->index++;
                                            PRINT_INDEX(TrieCheckFrom,4,5)
                                        }
                                        NextCand->do_it = FALSE;
                                    }
                                }
                                break;
                            case 1:
                                TrieCheckFrom->index ++;
                                PRINT_INDEX(TrieCheckFrom,4,6)
                                tv->gotonode = TrieCheckFrom;
                                break;
                            case 2:
                                if (TrieCandFrom) {
                                    TrieCheckFrom->index += TrieCandFrom->index;
                                    TrieCandFrom->goes_to = TrieCheckFrom;
                                    PRINT_INDEX(TrieCheckFrom,4,7)
                                }
                                else {
                                    TrieCheckFrom->index++;
                                    PRINT_INDEX(TrieCheckFrom,4,8)
                                }
                                if (temp == tv->maxtreelevel) {
                                    tmp1 = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                    for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == tmp1) break;
                                    if (i == AutomCount[0]) AutomCount[AutomCount[0]++] = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                }
                                break;
                            default:
                                break;
                        }
                        return temp;
                    }
                }
            }
            CheckAutList = CheckAutList->next;
        }
        CheckLevel++;
    }
    return FALSE;
}


int CheckForSingAutomorphisms(Candidate *CurrCand, Partition *NextPart, Candidate *NextCand,
                              struct TracesVars* tv, struct TracesInfo* ti,
                              int m, int n) {
    int i, j, temp, tmp, tmp1, result, tgt_level, numtemporbits;
    TracesSpine *SpineTL;
    Candidate *CheckAutList;
    searchtrie *TrieCandFrom, *TrieCheckFrom;
    
    SpineTL = Spine+tv->tolevel;
    CheckAutList = SpineTL->liststart;
    tv->gotonode = NULL;
    
    result = 0;
    while (CheckAutList != NULL) {
        if (CheckAutList->do_it && (CheckAutList->stnode->father == CurrCand->stnode)) {
            if (CheckAutList->firstsingcode == NextCand->firstsingcode) {
                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                if ((tv->tolevel == 1) && (Spine[0].part->cells == 1)) tmp = 2; else tmp = tv->tolevel;
                if (TreeFyTwo(tmp, CheckAutList, NextCand, NextPart, n, tv, ti)) {
                    if (isautom_sg((graph*)tv->input_graph, AUTPERM, tv->options->digraph, m, n)) {
                        if (!findperm(gensB, AUTPERM, n)) {
                            if (tv->options->verbosity >= 2) tv->schreier3 -= CPUTIME;
                            addgenerator(&gpB, &gensB, AUTPERM, n);
                            if (tv->options->verbosity >= 2) tv->schreier3 += CPUTIME;
                            result = CheckAutList->name;
                            if (TempOrbits) {
                                if (tv->compstage == 0) {
                                    for (j=0; j<tv->permInd; j++) {
                                        orbjoin_sp_pair(TempOrbits, TempOrbList, n,
                                                        PrmPairs[j].arg, PrmPairs[j].val, &numtemporbits);
                                    }
                                }
                                else {
                                    orbjoin(TempOrbits, AUTPERM, n);
                                }
                            }
                            tv->stats->numgenerators++;
                            orbjoin_sp_perm(tv->orbits, AUTPERM, OrbList, n, &tv->stats->numorbits);
                            
                            ti->thegrouphaschanged = TRUE;
                            ti->identitygroup = FALSE;
                            tv->lev_of_lastauto = tv->tolevel;
                            if (tv->options->verbosity >= 2) fprintf(outfile, "[a(%d)] ", CheckAutList->name);
                            if (tv->options->verbosity >= 2 && tv->options->writeautoms) {
                                PRINT_RETURN
                            }
                            if (tv->options->writeautoms) {
                                fprintf(outfile, "Gen(a) #%d: ", tv->stats->numgenerators);
                                writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n);
                            }
                            if (tv->options->userautomproc) {
                                (*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n);
                            }
                        }
                        else {
                            if (tv->options->verbosity >= 2) {
                                fprintf(outfile, "[a*]");
                            }
                            if (TempOrbits) {
                                if (tv->compstage == 0) {
                                    for (j=0; j<tv->permInd; j++) {
                                        orbjoin_sp_pair(TempOrbits, TempOrbList, n,
                                                        PrmPairs[j].arg, PrmPairs[j].val, &numtemporbits);
                                    }
                                }
                                else {
                                    orbjoin(TempOrbits, AUTPERM, n);
                                }
                            }
                            else {
                                orbjoin(tv->currorbit, AUTPERM, n);
                            }
                            result = -CheckAutList->name;
                        }
                        
                        TrieCandFrom = NULL;
                        TrieCheckFrom = CheckAutList->stnode;
                        if (CurrCand->stnode->level <= 1) {
                            tgt_level = CurrCand->stnode->level + 1;
                            while (TrieCheckFrom->level > tgt_level) {
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if (tv->tolevel <= TrieCheckFrom->level) {
                                tgt_level = tv->tolevel;
                                while (TrieCheckFrom->level != tgt_level) {
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            }
                            else {
                                TrieCandFrom = CurrCand->stnode;
                                tgt_level = TrieCheckFrom->level;
                                while (TrieCandFrom->level != tgt_level) {
                                    TrieCandFrom = TrieCandFrom->father;
                                }
                            }
                        }
                        if (TrieCandFrom) {
                            while (TrieCandFrom->father != TrieCheckFrom->father) {
                                TrieCandFrom = TrieCandFrom->father;
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if ((TrieCheckFrom->level > 1) && (TrieCheckFrom->father != CurrCand->stnode)) {
                                TrieCandFrom = CurrCand->stnode;
                                TrieCheckFrom = TrieCheckFrom->father;
                                while (TrieCandFrom->father != TrieCheckFrom->father) {
                                    TrieCandFrom = TrieCandFrom->father;
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            }
                        }
                        
                        while (TrieCheckFrom->goes_to) {
                            TrieCheckFrom = TrieCheckFrom->goes_to;
                        }
                        
                        for (temp=1; temp<=tv->tolevel; temp++) {
                            if (CheckAutList->lab[Spine[temp].tgtpos] != NextCand->lab[Spine[temp].tgtpos]) {
                                break;
                            }
                        }
                        
                        switch (tv->compstage) {
                            case 0:
                                if (tv->strategy && (tv->steps == 1)) {
                                    RemoveFromLevel(temp, tv->maxtreelevel-1, tv->strategy, FALSE);
                                    if (TrieCandFrom) {
                                        TrieCheckFrom->index += TrieCandFrom->index;
                                        TrieCandFrom->goes_to = TrieCheckFrom;
                                        PRINT_INDEX(TrieCheckFrom,4,9)
                                    }
                                    else {
                                        TrieCheckFrom->index++;
                                        PRINT_INDEX(TrieCheckFrom,4,10)
                                    }
                                    NextCand->do_it = FALSE;
                                }
                                else {
                                    if (CheckAutList->lab[Spine[temp].tgtpos] >= NextCand->lab[Spine[temp].tgtpos]) {
                                        CheckAutList->do_it = FALSE;
                                        if (TrieCandFrom) {
                                            TrieCandFrom->index += TrieCheckFrom->index;
                                            tv->newindex = 0;
                                            TrieCheckFrom->goes_to = TrieCandFrom;
                                            PRINT_INDEX(TrieCandFrom,4,11)
                                        }
                                        else {
                                            if (CurrCand->stnode->level > 1) {
                                                tv->newgotonode = TrieCheckFrom;
                                                tv->newindex = TrieCheckFrom->index;
                                            }
                                            else {
                                                tv->newgotonode = NULL;
                                                tv->newindex = 0;
                                            }
                                        }
                                    }
                                    else {
                                        if (TrieCandFrom) {
                                            TrieCheckFrom->index += TrieCandFrom->index;
                                            TrieCandFrom->goes_to = TrieCheckFrom;
                                            PRINT_INDEX(TrieCheckFrom,4,12)
                                        }
                                        else {
                                            TrieCheckFrom->index++;
                                            PRINT_INDEX(TrieCheckFrom,4,13)
                                        }
                                        NextCand->do_it = FALSE;
                                    }
                                }
                                break;
                            case 1:
                                TrieCheckFrom->index++;
                                PRINT_INDEX(TrieCheckFrom,4,14)
                                tv->gotonode = TrieCheckFrom;
                                break;
                            case 2:
                                if (TrieCandFrom) {
                                    TrieCheckFrom->index += TrieCandFrom->index;
                                    TrieCandFrom->goes_to = TrieCheckFrom;
                                    PRINT_INDEX(TrieCheckFrom,4,15)
                                }
                                else {
                                    TrieCheckFrom->index++;
                                    PRINT_INDEX(TrieCheckFrom,4,16)
                                }
                                if (temp == tv->maxtreelevel) {
                                    tmp1 = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                    for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == tmp1) break;
                                    if (i == AutomCount[0]) AutomCount[AutomCount[0]++] = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                }
                                break;
                            default:
                                break;
                        }
                        return result;
                    }
                }
            }
        }
        CheckAutList = CheckAutList->next;
    }
    return result;
}

int CheckForMatching(Candidate *CurrCand, Candidate *NextCand, Partition *Part, struct TracesVars* tv, struct TracesInfo* ti, int m, int n) {
    int i, j, vtx, vtx1, temp, tmp1, tgt_level, numtemporbits, pos;
    TracesSpine *SpineTL;
    Candidate *CheckAutList;
    int *cls;
    searchtrie *TrieCandFrom, *TrieCheckFrom;
    boolean CodeVerify;
    
    SpineTL = Spine+tv->tolevel;
    CheckAutList = SpineTL->liststart;
    cls = Part->cls;
    numtemporbits = 0;
    tv->gotonode = NULL;
    
    while (CheckAutList != NULL) {
        if (CheckAutList->do_it && (CheckAutList->singcode == NextCand->singcode)) {
            TrieCheckFrom = CheckAutList->stnode->father;
            TrieCandFrom = CurrCand->stnode;
            while (TrieCandFrom != TrieCheckFrom) {
                TrieCandFrom = TrieCandFrom->father;
                TrieCheckFrom = TrieCheckFrom->father;
            }
            
            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
            
            CodeVerify = TRUE;
            for (i=Spine[TrieCheckFrom->level+1].singstart; i<SpineTL->singend; i++) {
                Markers[NextCand->lab[Singletons[i]]] = tv->mark;
            }
            for (i=Spine[TrieCheckFrom->level+1].singstart; i<SpineTL->singend; i++) {
                pos = Singletons[i];
                vtx1 = CheckAutList->lab[pos];
                if (Markers[vtx1] != tv->mark) {
                    CodeVerify = FALSE;
                    break;
                }
                vtx = NextCand->lab[pos];
                SETPAIRSAUT(vtx, vtx1)
                
                if (tv->preprocessed && Diff[vtx1]) {
                    MakeTree(vtx, vtx1, tv->input_graph, n, tv, TRUE);
                }
            }
            tv->conta7++;
            
            if (CodeVerify) {
                if (isautom_sg_pair((graph*)tv->input_graph, AUTPERM, tv->options->digraph, m, n, tv)) {
                    
                    if (!findperm(gensB, AUTPERM, n)) {
                        if (tv->options->verbosity >= 2) tv->schreier3 -= CPUTIME;
                        if (tv->options->generators) addpermutation(&gensB, AUTPERM, n);
                        
                        if (tv->options->verbosity >= 2) tv->schreier3 += CPUTIME;
                        
                        if (CheckAutList->stnode->father == CurrCand->stnode) {
                            if (TempOrbits) {
                                if (tv->compstage == 0) {
                                    for (j=0; j<tv->permInd; j++) {
                                        orbjoin_sp_pair(TempOrbits, TempOrbList, n, PrmPairs[j].arg, PrmPairs[j].val, &numtemporbits);
                                    }
                                }
                                else {
                                    orbjoin(TempOrbits, AUTPERM, n);
                                }
                            } else {
                                orbjoin(tv->currorbit, AUTPERM, n);
                            }
                        }
                        
                        tv->stats->numgenerators++;
                        
                        for (j=0; j<tv->permInd; j++) {
                            orbjoin_sp_pair(tv->orbits, OrbList, n, PrmPairs[j].arg, PrmPairs[j].val, &tv->stats->numorbits);
                        }
                        
                        ti->thegrouphaschanged = TRUE;
                        ti->first_matching = TRUE;
                        tv->lev_of_lastauto = tv->tolevel;
                        if (tv->options->verbosity >= 2) fprintf(outfile, "[M(%d)] ", CheckAutList->name);
                        if (tv->options->verbosity >= 2 && tv->options->writeautoms) {
                            PRINT_RETURN
                        }
                        if (tv->options->writeautoms) {
                            fprintf(outfile, "Gen(M) #%d: ", tv->stats->numgenerators);
                            writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n);
                        }
                        if (tv->options->userautomproc) {
                            (*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n);
                        }
                        
                    }
                    else {
                        if (tv->options->verbosity >= 2) {
                            fprintf(outfile, "[M*]");
                        }
                        if (TempOrbits) {
                            if (tv->compstage == 0) {
                                for (j=0; j<tv->permInd; j++) {
                                    orbjoin_sp_pair(TempOrbits, TempOrbList, n,
                                                    PrmPairs[j].arg, PrmPairs[j].val, &numtemporbits);
                                }
                            }
                            else {
                                orbjoin(TempOrbits, AUTPERM, n);
                            }
                        }
                    }
                    
                    TrieCandFrom = NULL;
                    TrieCheckFrom = CheckAutList->stnode;
                    if (CurrCand->stnode->level <= 1) {
                        tgt_level = CurrCand->stnode->level + 1;
                        while (TrieCheckFrom->level > tgt_level) {
                            TrieCheckFrom = TrieCheckFrom->father;
                        }
                    }
                    else {
                        if (tv->tolevel <= TrieCheckFrom->level) {
                            tgt_level = tv->tolevel;
                            while (TrieCheckFrom->level != tgt_level) {
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            TrieCandFrom = CurrCand->stnode;
                            tgt_level = TrieCheckFrom->level;
                            while (TrieCandFrom->level != tgt_level) {
                                TrieCandFrom = TrieCandFrom->father;
                            }
                        }
                    }
                    if (TrieCandFrom) {
                        while (TrieCandFrom->father != TrieCheckFrom->father) {
                            TrieCandFrom = TrieCandFrom->father;
                            TrieCheckFrom = TrieCheckFrom->father;
                        }
                    }
                    else {
                        if ((TrieCheckFrom->level > 1) && (TrieCheckFrom->father != CurrCand->stnode)) {
                            TrieCandFrom = CurrCand->stnode;
                            TrieCheckFrom = TrieCheckFrom->father;
                            while (TrieCandFrom->father != TrieCheckFrom->father) {
                                TrieCandFrom = TrieCandFrom->father;
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                    }
                    
                    while (TrieCheckFrom->goes_to) {
                        TrieCheckFrom = TrieCheckFrom->goes_to;
                    }
                    
                    for (temp=1; temp<=tv->tolevel; temp++) {
                        if (CheckAutList->lab[Spine[temp].tgtpos] != NextCand->lab[Spine[temp].tgtpos]) {
                            break;
                        }
                    }
                    switch (tv->compstage) {
                        case 0:
                            if (tv->strategy && (tv->steps == 1)) {
                                RemoveFromLevel(temp, tv->maxtreelevel-1, tv->strategy, FALSE);
                                if (TrieCandFrom) {
                                    TrieCheckFrom->index += TrieCandFrom->index;
                                    TrieCandFrom->goes_to = TrieCheckFrom;
                                    PRINT_INDEX(TrieCheckFrom,4,17)
                                }
                                else {
                                    TrieCheckFrom->index++;
                                    PRINT_INDEX(TrieCheckFrom,4,18)
                                }
                                NextCand->do_it = FALSE;
                            }
                            else {
                                if (CheckAutList->lab[Spine[temp].tgtpos] >= NextCand->lab[Spine[temp].tgtpos]) {
                                    CheckAutList->do_it = FALSE;
                                    if (TrieCandFrom) {
                                        TrieCandFrom->index += TrieCheckFrom->index;
                                        tv->newindex = 0;
                                        TrieCheckFrom->goes_to = TrieCandFrom;
                                        PRINT_INDEX(TrieCandFrom,4,19)
                                    }
                                    else {
                                        if (CurrCand->stnode->level > 1) {
                                            tv->newgotonode = TrieCheckFrom;
                                            tv->newindex = TrieCheckFrom->index;
                                        }
                                        else {
                                            tv->newgotonode = NULL;
                                            tv->newindex = 0;
                                        }
                                    }
                                }
                                else {
                                    if (TrieCandFrom) {
                                        TrieCheckFrom->index += TrieCandFrom->index;
                                        TrieCandFrom->goes_to = TrieCheckFrom;
                                        PRINT_INDEX(TrieCheckFrom,4,20)
                                    }
                                    else {
                                        TrieCheckFrom->index++;
                                        PRINT_INDEX(TrieCheckFrom,4,21)
                                    }
                                    NextCand->do_it = FALSE;
                                }
                            }
                            break;
                        case 1:
                            TrieCheckFrom->index ++;
                            tv->gotonode = TrieCheckFrom;
                            PRINT_INDEX(TrieCheckFrom,4,22)
                            break;
                        case 2:
                            if (TrieCandFrom) {
                                TrieCheckFrom->index += TrieCandFrom->index;
                                TrieCandFrom->goes_to = TrieCheckFrom;
                                PRINT_INDEX(TrieCheckFrom,4,23)
                            }
                            else {
                                TrieCheckFrom->index++;
                                PRINT_INDEX(TrieCheckFrom,4,24)
                            }
                            if (temp == tv->maxtreelevel) {
                                tmp1 = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == tmp1) break;
                                if (i == AutomCount[0]) AutomCount[AutomCount[0]++] = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                            }
                            break;
                        default:
                            break;
                    }
                    return temp;
                }
            }
        }
        CheckAutList = CheckAutList->next;
    }
    return FALSE;
}

void CodeClassify(int Level, int code, int cell) {
    switch (EPCodes[Level].info) {
        case 0:
            EPCodes[Level].code = code;
            EPCodes[Level].cell = cell;
            EPCodes[Level].info = 1;
            break;
        case 1:
            if (EPCodes[Level].cell != cell) {
                EPCodes[Level].info = 3;
            } else {
                if (EPCodes[Level].code != code) {
                    EPCodes[Level].info = 2;
                }
            }
            break;
        case 2:
            if (EPCodes[Level].cell != cell) {
                EPCodes[Level].info = 3;
            }
            break;
        default:
            break;
    }
}

void Complete(sparsegraph *sg_orig, Candidate *Cand, Partition *Part, int cell, TracesVars *tv,
              double *grpsize1, int *grpsize2, permnode **ring, int n) {
    int i, j, k;
    int arg, val;
    int numtemporbits;
    k = cell + Part->cls[cell];
    
    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
    for (i = cell; i < k; i++) {
        tv->currorbit[Cand->lab[i]] = Cand->lab[k];
        arg = Cand->lab[i];
        val = Cand->lab[i+1];
        SETPAIRSAUTANDTREE(arg, val)
    }
    arg = Cand->lab[i];
    val = Cand->lab[cell];
    SETPAIRSAUTANDTREE(arg, val)
    SPECIALGENERATORS
    
    if (Part->cls[cell] > 1) {
        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
        arg = Cand->lab[cell];
        val = Cand->lab[cell+1];
        SETPAIRSAUTANDTREE(arg, val)
        arg = Cand->lab[cell+1];
        val = Cand->lab[cell];
        SETPAIRSAUTANDTREE(arg, val)
        SPECIALGENERATORS
    }
}

int CompStage0(Partition *CurrPart, Partition *NextPart, Candidate *CurrCand, Candidate *NextCand,
               int m, int n, struct TracesVars* tv, struct TracesInfo *ti) {
    int i, j, i1, j2, k, cu, cu1, num_indv;
    int temp, tmp, auxcode, search_vtx, gom_level;
    boolean closeloop, firstsing, has_nexttcell;
    Candidate *SpTLliststart, *AuxCand;
    searchtrie *TreeNode, *TreeNode1, *TreeNode2;
    
#ifdef NAUTY_IN_MAGMA
    if (main_seen_interrupt) return NAUTY_KILLED;
#else
    if (nauty_kill_request) return NAUTY_KILLED;
#endif
    
    PRINT_FROM_VERB(4,tv->tolevel)
    if (TargetCell(CurrCand, CurrPart, n, tv, tv->tolevel)) {
        ++tv->tolevel;
        SpineTL = Spine+tv->tolevel;
        SpineTL->tgtcell = tv->tcell;
        SpineTL->tgtsize = CurrPart->cls[tv->tcell];
        SpineTL->tgtend = tv->tcell+SpineTL->tgtsize;
        SpineTL->tgtpos = SpineTL->tgtend - 1;
    }
    else {
        tv->finalnumcells = min(CurrPart->cells,tv->finalnumcells);    /* 160712 */
        ti->thereisnextlevel = SelectNextLevel(n, tv, ti);
        return 0;
    }
    
    tv->newgotonode = NULL;
    
    /*  CANDIDATE */
    temp = CurrCand->lab[Spine[1].tgtpos];
    k = SpineTL->tgtend;
    
    TreeNode = CurrCand->stnode;
    while (TreeNode) {
        if (TreeNode->goes_to) {
            CurrCand->do_it = FALSE;
            break;
        }
        TreeNode = TreeNode->father;
    }
    
    if (CurrCand->do_it) {
        if ((tv->orbits[temp] == temp) || tv->tolevel == 1) {
            ti->minimalinorbits = TRUE;
            
            if ((tv->group_level >= tv->tolevel) && (FixedBase(fix, tv, CurrCand, 0, tv->fromlevel))) {
                tv->nfix = tv->fromlevel;
                tv->currorbit = findcurrorbits(gpB, tv->nfix);
            } else {
                if ((!ti->identitygroup) &&
                    (((Spine[tv->fromlevel].liststart != Spine[tv->fromlevel].listend)
                      && (CurrPart->cls[tv->tcell] > 10))
                     || tv->strategy
                     || (tv->expathlength <=10)
                     )) {
                        
                        TempOrbits = NULL;
                        tv->samepref = FixBase(fix, tv, CurrCand, 0, tv->fromlevel);
                        if ((tv->samepref != tv->nfix) || ti->thegrouphaschanged) {
                            if (tv->options->verbosity >= 2) tv->schreier1 -= CPUTIME;
                            gom_level = getorbitsmin(fix, tv->nfix, gpB, &gensB, &tv->currorbit, CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell], n, TRUE);
                            if (tv->options->verbosity >= 2) tv->schreier1 += CPUTIME;
                            ti->thegrouphaschanged = FALSE;
                            
                            if (gom_level < tv->nfix) {
                                PRINT_NOTMIN_VERB(4)
                                
                                TreeNode = CurrCand->stnode;
                                j2 = CurrCand->lab[Spine[gom_level+1].tgtpos];
                                i1 = tv->currorbit[j2];
                                for (j=0; j < tv->nfix - gom_level; j++) {
                                    TreeNode = TreeNode->father;
                                }
                                TreeNode1 = TreeNode->first_child;
                                while (TreeNode1) {
                                    if (TreeNode1->vtx == i1) {
                                        break;
                                    }
                                    TreeNode1 = TreeNode1->next_sibling;
                                }
                                if (TreeNode1) {
                                    while (TreeNode1->goes_to) {
                                        TreeNode1 = TreeNode1->goes_to;
                                    }
                                    TreeNode2 = TreeNode1->next_sibling;
                                    while (TreeNode2->vtx != j2) {
                                        TreeNode2 = TreeNode2->next_sibling;
                                    }
                                    TreeNode1->index += TreeNode2->index;
                                    TreeNode2->goes_to = TreeNode1;
                                    PRINT_INDEX(TreeNode1,4,25)
                                    
                                    ti->minimalinorbits = FALSE;
                                }
                                else {
                                    tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
                                }
                            }
                        }
                        else {
                            tv->currorbit = findcurrorbits(gpB, tv->nfix);
                        }
                    }
                else {
                    TempOrbits = WorkArray1;
                    memcpy(TempOrbits, IDENTITY_PERM, n*sizeof(int));
                    memcpy(TempOrbList, IDENTITY_PERM, n*sizeof(int));
                    
                    tv->conta1++;
                    tv->currorbit = TempOrbits;
                }
            }
            if (ti->minimalinorbits) {
                memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
                memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
                tv->conta2++;
                auxcode = CurrCand->code;
                SpineTL->trcstart = CurrPart->cells;
                TheTrace[SpineTL->trcstart] = SpineTL->tgtpos;
                if (!CurrCand->sortedlab) {
                    quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
                    for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
                        CurrCand->invlab[CurrCand->lab[i]] = i;
                    }
                    CurrCand->sortedlab = TRUE;
                }
                
                tv->indivstart = tv->tcell+CurrCand->indnum;
                tv->indivend = tv->indivstart+tv->steps;
                if (tv->indivend > SpineTL->tgtend) {
                    tv->indivend = SpineTL->tgtend;
                }
                
                temp = CurrCand->lab[tv->indivstart];
                for (k = tv->indivstart; k < tv->indivend; k++) {
                    CurrCand->indnum++;
                    NextCand->singcode = CurrCand->singcode;
                    NextCand->vertex = CurrCand->lab[k];
                    NextCand->name = ++tv->name;
                    if (NextCand->name == (NAUTY_INFINITY-2)) {
                        NextCand->name = tv->name = 1;
                    }
                    
                    PRINT_INDIV_VERB(4,tv->tolevel)
                    if (tv->currorbit[NextCand->vertex] != NextCand->vertex) {
                        PRINT_SKIPPED_VERB(4)
                        
                        search_vtx = tv->currorbit[NextCand->vertex];
                        TreeNode = CurrCand->stnode;
                        if (TreeNode->first_child) {
                            TreeNode = TreeNode->first_child;
                            while (TreeNode) {
                                if (TreeNode->vtx == search_vtx) {
                                    break;
                                }
                                TreeNode = TreeNode->next_sibling;
                            }
                            if (TreeNode) {
                                while (TreeNode->goes_to) {
                                    TreeNode = TreeNode->goes_to;
                                }
                                TreeNode->index++;
                                PRINT_INDEX(TreeNode,4,26)
                                continue;
                            }
                        }
                    }
                    else {
                        PRINT_REFINE_VERB(4,'a')
                        
                        memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
                        memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
                        
                        tv->conta3++;
                        
                        if (NextPart->cls[tv->tcell] == 2) {
                            num_indv = 2;
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell]);
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell+1]);
                            if (SpineTL->singstart == SpineTL->singend) {
                                Singletons[SpineTL->singend++] = tv->tcell;
                                Singletons[SpineTL->singend++] = tv->tcell+1;
                            }
                        }
                        else {
                            num_indv = 1;
                            NextCand->singcode = MASHCOMM(NextCand->singcode, NextCand->vertex);
                            if (SpineTL->singstart == SpineTL->singend) {
                                Singletons[SpineTL->singend++] = tv->tcell + NextPart->cls[tv->tcell] - 1;
                            }
                        }
                        
                        Individualize(NextPart, NextCand, NextCand->vertex, tv->tcell, CurrPart->cells, SpineTL->tgtpos);
                        tv->stats->numnodes++;
                        tv->answ = traces_refine(NextCand,
                                                 n,
                                                 NextPart, tv, ti, num_indv, TRUE);
                        
                        switch (tv->answ) {
                            case 0:				/* Interrupted refinement: do not add to the list */
                                tv->stats->interrupted++;
                                SpineTL->levelcounter++;
                                break;
                            case 1 :			/* The same trace has been found once more : add to the list */
                                SpineTL->levelcounter++;
                                
                                NextCand->do_it = TRUE;
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel);
                                
                                tv->tolevel_tl = tv->tolevel;
                                NextCand->pathsingcode = NextCand->singcode;
                                NextCand->firstsingcode = 0;
                                
                                if (tv->steps > 1) {
                                    if (tv->fromlevel <= tv->lev_of_lastauto) {
                                        
                                        closeloop = CheckForMatching(CurrCand, NextCand, NextPart, tv, ti, m, n);
                                    }
                                    if (NextCand->do_it) {
                                        firstsing = TRUE;
                                        
                                        /* EXPERIMENTAL PATH */
                                        if (NextPart->cells != tv->finalnumcells) {    /* 160712 */
                                            if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                            while (NextPart->cells < n) {
                                                if (firstsing && BreakSteps[tv->tolevel]) {
                                                    firstsing = FALSE;
                                                    NextCand->firstsingcode = NextCand->pathsingcode;
                                                    if (CheckForSingAutomorphisms(CurrCand, NextPart, NextCand, tv, ti, m, n))
                                                        if (!NextCand->do_it) {
                                                            break;
                                                        }
                                                }
                                                
                                                has_nexttcell = TargetCellExpPath(NextCand, NextPart, tv);
                                                
                                                if (!has_nexttcell) {
                                                    NextCand->firstsingcode = NextCand->pathsingcode;
                                                    if (tv->options->verbosity >= 2) {
                                                        if (tv->tolevel_tl-tv->tolevel >= 6) {
                                                            PRINT_EXPPATHSTEP(NextCand, TRUE)
                                                        }
                                                        else {
                                                            fprintf(outfile, "(%d) ", tv->tolevel_tl);
                                                        }
                                                    }
                                                    break;
                                                }
                                                ExperimentalStep(NextPart, NextCand, tv, ti, m, n);
                                                if (NextPart->cells == n) {
                                                    has_nexttcell = FALSE;
                                                }
                                                PRINT_EXPPATHSTEP(NextCand, TRUE)
                                            }
                                            if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                        }
                                    }
                                    else {
                                        if (closeloop < tv->tolevel) k = SpineTL->tgtend;
                                        PRINT_RETURN
                                        break;
                                    }
                                    
                                    if (!tv->strategy && !tv->options->getcanon && (tv->tolevel_tl == tv->tolevel + 1) && ((NextPart->cells != tv->finalnumcells) || (NextPart->cells == n))) {    /* 160717 */
                                        tv->levelfromCS0 = tv->tolevel;
                                        tv->maxtreelevel = tv->tolevel_tl;
                                        if (tv->tolevel == 1) {
                                            tv->newst_stage1 = searchtrie_make(CurrCand, NextCand, n, tv);
                                            EXITFROMSTAGE0EXPATH1
                                        }
                                        else {
                                            temp = 0;
                                            for (i=0; i<tv->tolevel; i++) {
                                                temp += Spine[i].listcounter;
                                            }
                                            if (temp > 5) {
                                                tv->newst_stage1 = searchtrie_make(CurrCand, NextCand, n, tv);
                                                EXITFROMSTAGE0EXPATH1
                                            }
                                        }
                                    }
                                    
                                    /* ANY AUTOMORPHISM? */
                                    if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                                    tv->newindex = 0;
                                    if (NextCand->do_it) {
                                        closeloop = CheckForAutomorphisms(CurrCand, NextCand, tv, ti, m, n, NextPart);
                                        if (!NextCand->do_it && closeloop < tv->tolevel) k = SpineTL->tgtend;
                                    }
                                    if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                                    
                                    if (NextCand->do_it) {
                                        ADDTONEXTLEVEL;
                                        SpineTL->keptcounter++;
                                        searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                    }
                                }
                                else {
                                    if (BreakSteps[tv->tolevel]) {
                                        NextCand->firstsingcode = NextCand->pathsingcode;
                                        if (CheckForSingAutomorphisms(CurrCand, NextPart, NextCand, tv, ti, m, n))
                                            if (!NextCand->do_it) {
                                                PRINT_RETURN
                                                break;
                                            }
                                    }
                                    
                                    /* ANY AUTOMORPHISM? */
                                    if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                                    tv->newindex = 0;
                                    if (NextCand->do_it) {
                                        closeloop = CheckForAutomorphisms(CurrCand, NextCand, tv, ti, m, n, NextPart);
                                        if (!NextCand->do_it && closeloop < tv->tolevel) k = SpineTL->tgtend;
                                    }
                                    if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                                    
                                    if (NextCand->do_it) {
                                        ADDTONEXTLEVEL;
                                        SpineTL->keptcounter++;
                                        searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                    }
                                }
                                PRINT_RETURN
                                break;
                            case 2 :	/* Delete the old list and start a new one: a better trace has been found */
                                
                                tv->tolevel_tl = tv->tolevel;
                                has_nexttcell = FALSE;
                                if (NextPart->cells == n) {
                                    tv->stats->canupdates++;
                                    if (tv->options->usercanonproc != NULL)
                                    {
                                        (*tv->options->usercanonproc)((graph*)tv->input_graph, NextCand->lab, (graph*)tv->cangraph, tv->stats->canupdates, NextCand->code, m, n);
                                    }
                                }
                                
                                if (tv->tolevel > tv->treedepth) {
                                    tv->treedepth = tv->tolevel;
                                    if (tv->strategy) {
                                        SpineTL->part = NewPartition(n);
                                    }
                                    else {
                                        NewPartSpine(tv->tolevel,n);
                                    }
                                }
                                
                                if (!tv->strategy && (tv->tolevel > 1) && !SpineTL->liststart) {
                                    /* First Candidate at current level */
                                    tv->maxtreelevel = tv->tolevel;
                                    
                                    SpineTL->liststart = NewCandidate(n, &GarbList, TRUE);
                                    SpineTL->listend = SpineTL->liststart;
                                    
                                    tv->conta0++;
                                    CopyCand(SpineTL->liststart, NextCand, n, TEMPLAB, TEMPINVLAB);
                                    if (NextPart->cells < tv->finalnumcells) SpineTL->liststart->code = auxcode;
                                    COPYPART(SpineTL->part, NextPart);
                                    tv->newindex = 0;
                                    tv->newst_stage1 = searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                    
                                    SpineTL->listcounter = 1;
                                    SpTLliststart = SpineTL->liststart;
                                    
                                    i = tv->tolevel;
                                    if (tv->brkstpcount) {
                                        while ((i<n) && !BreakSteps[i]) {
                                            i++;
                                        }
                                        if (i<n) SpineTL->liststart->firstsingcode = Spine[i].singcode;
                                    }
                                    
                                    SpineTL->updates = 1;
                                    SpineTL->levelcounter = 1;
                                    SpineTL->keptcounter = 1;
                                    
                                    PRINT_LINE_PLUS(tv->fromlevel)
                                    
                                    if (tv->options->verbosity >= 2) PRINT_CANDIDATE(SpineTL->liststart, tv->tolevel);
                                    PRINT_RETURN;
                                    
                                    if (!tv->strategy && !tv->options->getcanon && (tv->tolevel+1 == tv->firstpathlength) && ((NextPart->cells != tv->finalnumcells) || (NextPart->cells == n))) {
                                        if ((tv->tolevel == 1) && (CurrPart->cls[tv->tcell] > 5)) {
                                            EXITFROMSTAGE0EXPATH2;
                                        }
                                        else {
                                            temp = 0;
                                            for (i=0; i<tv->tolevel; i++) {
                                                temp += Spine[i].listcounter;
                                            }
                                            if (temp > 5) {
                                                EXITFROMSTAGE0EXPATH2;
                                            }
                                        }
                                    }
                                }
                                else {
                                    memset(WorkArray, 0, n*sizeof(int));
                                    
                                    tv->lastcell = tv->lastlev = -1;
                                    has_nexttcell = TargetCellFirstPath(NextCand, NextPart, tv);
                                    
                                    if (!has_nexttcell) {
                                        tv->stats->canupdates++;
                                        if (tv->options->usercanonproc != NULL) {
                                            (*tv->options->usercanonproc)((graph*)tv->input_graph, NextCand->lab, (graph*)tv->cangraph, tv->stats->canupdates, NextCand->code, m, n);
                                        }
                                    }
                                    
                                    tv->tcellevel = tv->maxtreelevel = tv->tolevel;
                                    SpineTL->levelcounter++;
                                    SpineTL->updates++;
                                    SpineTL->keptcounter = 1;
                                    
                                    RemoveFromLevel(tv->tolevel, tv->treedepth, tv->strategy, TRUE);
                                    SpineTL->liststart = NewCandidate(n, &GarbList, TRUE);
                                    SpineTL->listend = SpineTL->liststart;
                                    
                                    tv->conta0++;
                                    CopyCand(SpineTL->liststart, NextCand, n, NULL, NULL);
                                    COPYPART(SpineTL->part, NextPart);
                                    
                                    tv->newindex = 0;
                                    
                                    tv->newst_stage1 = searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                    
                                    SpineTL->listcounter = 1;
                                    SpTLliststart = SpineTL->liststart;
                                    
                                    SpTLliststart->pathsingcode = SpineTL->singcode = SpTLliststart->singcode;
                                    SpTLliststart->firstsingcode = 0;
                                    
                                    PRINT_LINE
                                    if (tv->options->verbosity >= 2) PRINT_CANDIDATE(SpTLliststart, tv->tolevel);
                                    
                                    memset(BreakSteps, 0, n*sizeof(int));
                                    tv->brkstpcount = 0;
                                    
                                    if (tv->steps > 1) {
                                        
                                        /* EXPERIMENTAL PATH */
                                        if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                        PRINTF2("CStage0 2: %d\n", tv->finalnumcells);
                                        tv->finalnumcells = n;
                                        
                                        while (has_nexttcell) {
                                            ExperimentalStep(NextPart, SpTLliststart, tv, ti, m, n);
                                            
                                            Spine[tv->tolevel_tl].singcode = SpTLliststart->pathsingcode;
                                            has_nexttcell = TargetCellFirstPath(SpTLliststart, NextPart, tv);
                                            PRINT_EXPPATHSTEP(SpTLliststart, TRUE)
                                        }
                                        if (NextPart->cells < n) {
                                            PRINTF2("CStage0 3: %d\n", tv->finalnumcells);
                                            tv->finalnumcells = min(NextPart->cells,tv->finalnumcells);    /* 160712 */
                                            PRINTF2("CStage0 3<: %d\n", tv->finalnumcells);
                                        }
                                        
                                        PRINTF2("CS0 2?: finalnumcells: %d\n", tv->finalnumcells);
                                        if (NextPart->cells == tv->finalnumcells) {
                                            UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                        }
                                        
                                        if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                        
                                        tv->firstpathlength = tv->tolevel_tl;
                                        PRINT_RETURN
                                        if (!tv->strategy && !tv->options->getcanon && (NextPart->cells == tv->finalnumcells) && (tv->tolevel_tl == tv->tolevel + 1) && ((NextPart->cells != tv->finalnumcells) || (NextPart->cells == n))) {
                                            tv->maxtreelevel = tv->tolevel_tl;
                                            if ((tv->tolevel == 1) && (CurrPart->cls[tv->tcell] > 5)) {
                                                EXITFROMSTAGE0EXPATH2
                                            }
                                            else {
                                                temp = 0;
                                                for (i=0; i<tv->tolevel; i++) {
                                                    temp += Spine[i].listcounter;
                                                }
                                                if (temp > 5) {
                                                    EXITFROMSTAGE0EXPATH2
                                                }
                                            }
                                        }
                                        memcpy(TEMPLAB, SpTLliststart->lab, n*sizeof(int));
                                        memcpy(TEMPINVLAB, SpTLliststart->invlab, n*sizeof(int));
                                        tv->conta5++;
                                    }
                                    else {
                                        PRINT_RETURN
                                    }
                                }
                                
                                break;
                            default:
                                break;
                        }
                    }
                } /* end for */
            }
        }
    }
    
    /* REMOVE CURRENT CANDIDATE */
    if (SpineFL->liststart && (k >= SpineTL->tgtend)) {
        SpineFL->liststart = CurrCand->next;
        if (CurrCand->next == NULL) {
            SpineFL->listend = NULL;
        }
        SpineFL->listcounter--;
        CurrCand->next = GarbList;
        GarbList = CurrCand;
    }
    ti->thereisnextlevel = SelectNextLevel(n, tv, ti);
    return 0;
}

int CompStage1(Partition *CurrPart, Partition *NextPart, Candidate *CurrCand, Candidate *NextCand,
               int m, int n,
               struct TracesVars* tv, struct TracesInfo *ti) {
    int i, k, cu, cu1, tmp, gom_level, search_vtx, temp;
    searchtrie *TreeNode, *TrieNode;
    
#ifdef NAUTY_IN_MAGMA
    if (main_seen_interrupt) return NAUTY_KILLED;
#else
    if (nauty_kill_request) return NAUTY_KILLED;
#endif
    
    CurrCand->stnode = tv->newst_stage1;
    
    tv->tolevel++;
    SpineTL = Spine+tv->tolevel;
    tv->tcell = SpineTL->tgtcell;
    SpineTL->levelcounter = 0;
    SpineTL->keptcounter = 0;
    SpineTL->updates = 1;
    
    
    if (tv->options->verbosity >= 2) {
        LINE(32, "=")
        NEXTLINE
    }
    
    memset(RefCells, 0, n*sizeof(int));
    memset(MultRefCells, 0, n*sizeof(int));
    ti->thegrouphaschanged = TRUE;
    
    /*  CANDIDATE */
    memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
    memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
    NextCand->do_it = TRUE;
    SpineTL->trcstart = CurrPart->cells;
    
    tv->indivstart = tv->tcell;
    tv->indivend = SpineTL->tgtend;
    if (TheGraph[CurrCand->lab[tv->indivstart]].d == 1) {
        tv->indivstart = SpineTL->tgtend-1;
    }
    
    FixBase(fix, tv, NextCand, 0, tv->fromlevel);
    
    if (!ti->identitygroup) {
        if (tv->options->verbosity >= 2) tv->schreier2 -= CPUTIME;
        tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
        if (tv->options->verbosity >= 2) tv->schreier2 += CPUTIME;
    }
    else {
        if (n / CurrPart->cls[tv->tcell] < 256) {
            memcpy(tv->currorbit, IDENTITY_PERM, n*sizeof(int));
        }
        else {
            for (k = tv->indivstart; k < tv->indivend; k++) {
                tv->currorbit[CurrCand->lab[k]] = CurrCand->lab[k];
            }
        }
    }
    
    if (!CurrCand->sortedlab) {
        quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
        for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
            CurrCand->invlab[CurrCand->lab[i]] = i;
        }
        CurrCand->sortedlab = TRUE;
    }
    for (k = tv->indivstart; k < tv->indivend; k++) {
        NextCand->vertex = CurrCand->lab[k];
        NextCand->name = ++tv->name;
        if (NextCand->name == (NAUTY_INFINITY-2)) {
            NextCand->name = tv->name = 1;
        }
        if (tv->currorbit[CurrCand->lab[k]] != CurrCand->lab[k]) {
            search_vtx = tv->currorbit[NextCand->vertex];
            TreeNode = CurrCand->stnode;
            if (TreeNode->first_child) {
                TreeNode = TreeNode->first_child;
                while (TreeNode) {
                    if (TreeNode->vtx == search_vtx) {
                        break;
                    }
                    TreeNode = TreeNode->next_sibling;
                }
                if (TreeNode) {
                    while (TreeNode->goes_to) {
                        TreeNode = TreeNode->goes_to;
                    }
                    TreeNode->index++;
                    PRINT_INDEX(TreeNode,4,27)
                }
            }
            continue;
        }
        PRINT_REFINE_VERB(4,'b')
        memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
        memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
        
        Individualize(NextPart, NextCand, CurrCand->lab[k], tv->tcell, CurrPart->cells, SpineTL->tgtpos);
        
        tv->stats->numnodes++;
        SpineTL->levelcounter++;
        tv->tolevel_tl = tv->tolevel;
        trieref = trieroot;
        SpineTL->levelcounter++;
        
        traces_refine_maketrie(NextCand,
                               n,
                               NextPart, tv, ti);
        
        RefCells[CurrCand->lab[k]] = NextPart->cells;
        PRINTF2("CS1 1?: finalnumcells: %d\n", tv->finalnumcells);
        if ((NextPart->cells == tv->finalnumcells) || (NextPart->cells == n)) {
            if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel);
            
            /* ANY AUTOMORPHISM? */
            if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
            
            PRINTF2("CS1 2?: finalnumcells: %d\n", tv->finalnumcells);
            CheckForAutomorphisms(CurrCand, NextCand, tv, ti, m, n, NextPart);
            if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
            
            PRINT_RETURN
            
            /* ADD TO NEXT LEVEL */
            SpineTL->keptcounter++;
            if (!Spine[tv->tolevel].listend) COPYPART(Spine[tv->tolevel].part, NextPart);
            ADDTONEXTLEVEL;
            searchtrie_make(CurrCand, SpineTL->listend, n, tv);
        }
    } /* end for */
    PRINTF2("CS1 3: finalnumcells: %d\n", tv->finalnumcells);
    for (k = tv->indivstart; k < tv->indivend; k++) {
        MultRefCells[RefCells[tv->currorbit[CurrCand->lab[k]]] % n]++;
    }
    
    if (tv->options->verbosity >= 2) {
        if (MultRefCells[0]) {
            fprintf(outfile, tv->digstring, n);
            fprintf(outfile, "cells: %d; ", MultRefCells[0]);
        }
        for (k=1; k<n; k++) {
            if (MultRefCells[k]) {
                fprintf(outfile, tv->digstring, k);
                fprintf(outfile, "cells: %d; ", MultRefCells[k]);
            }
        }
        NEXTLINE
        
    }
    
#if !MAXN
    DYNALLOC1(searchtrie*, RefPath, RefPath_sz, tv->tolevel, "Traces-CS1");
#endif
    
    TreeNode = CurrCand->stnode;
    while (TreeNode) {
        RefPath[TreeNode->level] = TreeNode;
        TreeNode = TreeNode->father;
    }
    
    /* REMOVE CURRENT CANDIDATE */
    SpineFL->liststart = CurrCand->next;
    if (CurrCand->next == NULL) {
        SpineFL->listend = NULL;
        SpineFL->listcounter = 1;
    }
    SpineFL->listcounter--;
    CurrCand->next = GarbList;
    GarbList = CurrCand;
    
    if (tv->options->verbosity >= 2) {
        LINE(32, "=")
        NEXTLINE
    }
    tv->compstage = 2;
    tv->steps = n;
    
    if (tv->options->verbosity >= 2) tv->schreier1 -= CPUTIME;
    
    gom_level = getorbitsmin(fix, tv->nfix, gpB, &gensB, &tv->currorbit,
                             CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell], n, TRUE);
    if (tv->options->verbosity >= 2) tv->schreier1 += CPUTIME;
    ORBITSIZES
    ti->thereisnextlevel = SelectNextLevel(n, tv, ti);
    PRINTF2("CS1 4: finalnumcells: %d\n", tv->finalnumcells);
    SpineTL->part->cells = tv->finalnumcells;
    
    AutomCount[0] = 2;
    AutomCount[1] = CurrCand->vertex;
    
    return 0;
}

int CompStage2(Partition *CurrPart, Partition *NextPart, Candidate *CurrCand, Candidate *NextCand,
               int m, int n,
               struct TracesVars* tv, struct TracesInfo *ti) {
    int i, j, i1, j2, k, cu, cu1, vertex, gom_level;
    int temp, tmp, autom;
    Candidate *AuxCand;
    searchtrie *TreeNode, *TreeNode1, *TreeNode2;
    int *CuOrb,*AuxOrb;
    boolean has_nexttcell = FALSE;
    searchtrie *TrieNode;
    boolean schreierwrong;
    
#ifdef NAUTY_IN_MAGMA
    if (main_seen_interrupt) return NAUTY_KILLED;
#else
    if (nauty_kill_request) return NAUTY_KILLED;
#endif
    
    autom = 0;
    schreierwrong = FALSE;
    
    TreeNode = CurrCand->stnode;
    tv->cand_level = 0;
    
    while (TreeNode) {
        if (TreeNode->goes_to) {
            CurrCand->do_it = FALSE;
        }
        if (!tv->cand_level && TreeNode == RefPath[TreeNode->level]) {
            tv->cand_level = TreeNode->level;
        }
        TreeNode = TreeNode->father;
    }
    if (tv->cand_level+1 == tv->maxtreelevel) {
        ti->useTempOrbits1 = TRUE;
    }
    else {
        ti->useTempOrbits1 = FALSE;
    }
    if (tv->cand_level == tv->fromlevel) {
        ti->useTempOrbits2 = TRUE;
    }
    else {
        ti->useTempOrbits2 = FALSE;
    }
    
    PRINT_FROM_VERB(4,tv->tolevel)
    
    if (CurrCand->do_it) {
        if (tv->tolevel == 0) {
            tv->fromlevel = tv->tolevel;
            SpineFL = Spine+tv->fromlevel;
            vertex = Spine[tv->maxtreelevel+1].liststart->lab[Spine[1].tgtpos];
            k = n;
            
            if (TargetCell(CurrCand, CurrPart, n, tv, tv->tolevel)) {
                ++tv->tolevel;
                SpineTL = Spine+tv->tolevel;
                SpineTL->tgtcell = tv->tcell;
                SpineTL->tgtsize = CurrPart->cls[tv->tcell];
                SpineTL->tgtend = tv->tcell+SpineTL->tgtsize;
                SpineTL->tgtpos = SpineTL->tgtend - 1;
            }
            else {
                PRINTF2("CStage2 1: %d\n", tv->finalnumcells);
                tv->finalnumcells = min(CurrPart->cells,tv->finalnumcells);    /* 160712 */
                return 0;
            }
            
            memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
            memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
            SpineTL->trcstart = CurrPart->cells;
            TheTrace[SpineTL->trcstart] = SpineTL->tgtpos;
            
            tv->indivstart = tv->tcell+CurrCand->indnum;
            tv->indivend = tv->indivstart+tv->steps;
            if (tv->indivend > SpineTL->tgtend) {
                tv->indivend = SpineTL->tgtend;
            }
            memset(CurrRefCells, 0, n*sizeof(int));
            ti->thegrouphaschanged = TRUE;
            
            if (!CurrCand->sortedlab) {
                quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
                for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
                    CurrCand->invlab[CurrCand->lab[i]] = i;
                }
                CurrCand->sortedlab = TRUE;
            }
            
            for (k = tv->indivstart; k < tv->indivend; k++) {
                if ((tv->orbits[CurrCand->lab[k]] == CurrCand->lab[k]) && ((tv->finalnumcells < n) || (OrbSize[tv->orbits[CurrCand->lab[k]]] >= OrbSize[tv->orbits[vertex]]))) {
                    
                    CurrCand->indnum++;
                    NextCand->singcode = CurrCand->singcode;
                    NextCand->vertex = CurrCand->lab[k];
                    NextCand->name = ++tv->name;
                    if (NextCand->name == (NAUTY_INFINITY-2)) {
                        NextCand->name = tv->name = 1;
                    }
                    
                    if (ti->thegrouphaschanged) {
                        if (tv->fromlevel == tv->maxtreelevel) {
                            CURRORBITSIZES
                        }
                        ti->thegrouphaschanged = FALSE;
                    }
                    
                    if (tv->currorbit[CurrCand->lab[k]] != CurrCand->lab[k]) {
                        continue;
                    }
                    
                    memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
                    memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
                    if (NextPart->cls[tv->tcell] == 2) {
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell]);
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell+1]);
                    }
                    else {
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[k]+labelorg);
                    }
                    
                    Individualize(NextPart, NextCand, CurrCand->lab[k], tv->tcell, CurrPart->cells, SpineTL->tgtpos);
                    
                    tv->stats->numnodes++;
                    Spine[tv->tolevel+1].levelcounter++;
                    if (tv->fromlevel == tv->maxtreelevel) {
                        tv->tolevel_tl = tv->tolevel;
                        trieref = trieroot;
                        
                        tv->answ = traces_refine_comptrie(NextCand,
                                                          n,
                                                          NextPart, tv, ti);
                        if (tv->answ) {
                            if (NextPart->cells != tv->finalnumcells) {
                                CurrRefCells[NextPart->cells % n] += CurrOrbSize[CurrCand->lab[k]];
                                if (CurrRefCells[NextPart->cells % n] > MultRefCells[NextPart->cells % n]) {
                                    k = n;
                                    break;
                                }
                                continue;
                            }
                            if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                                }
                    }
                    else {
                        tv->answ = traces_refine_sametrace(NextCand,
                                                           n,
                                                           NextPart, tv, ti);
                        
                        if (tv->answ) {
                            if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                                if (tv->tolevel == tv->maxtreelevel) {
                                    tv->tolevel_tl = tv->tolevel;
                                    if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                    TargetCellExpPath(NextCand, NextPart, tv);
                                    ExperimentalStep(NextPart, NextCand, tv, ti, m, n);
                                    PRINT_EXPPATHSTEP(NextCand, tv->answ)
                                    PRINTF2("CS2 1?: finalnumcells: %d\n", tv->finalnumcells);
                                    if ((NextPart->cells == tv->finalnumcells) || (NextPart->cells == n)) {
                                        UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                    }
                                    
                                    if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                    if (!tv->answ) {
                                        PRINT_RETURN
                                    }
                                }
                        }
                    }
                    if (tv->answ) {
                        PRINTF2("CS2 2?: finalnumcells: %d\n", tv->finalnumcells);
                        if ((NextPart->cells == tv->finalnumcells) || (NextPart->cells == n)) {
                            if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                            temp = (tv->tolevel_tl == tv->tolevel+1);
                            autom = CheckForAutomorphisms(CurrCand, NextCand,
                                                          tv, ti, temp, n, NextPart);
                            if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                            
                            if (ti->thegrouphaschanged) {
                                ORBITSIZES
                            }
                        }
                        PRINT_RETURN
                        
                        /* ADD TO NEXT LEVEL */
                        PRINTF2_2("CS2 3?: cells: %d, finalnumcells: %d\n", NextPart->cells, tv->finalnumcells);
                        if ((NextPart->cells != tv->finalnumcells) || (tv->tolevel != tv->maxtreelevel) || (tv->tolevel_tl != tv->tolevel+1)) {
                            ADDTONEXTLEVEL;
                            searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                        }
                    }
                    else {
                        tv->stats->interrupted++;
                    }
                    if (tv->fromlevel == tv->maxtreelevel) {
                        k = n;
                        break;
                    }
                }
            } /* end for */
        }
        else {
            
            temp = CurrCand->lab[Spine[1].tgtpos];
            vertex = Spine[tv->maxtreelevel+1].liststart->lab[Spine[1].tgtpos];
            k = n;
            
            if (tv->cand_level ||
                ((tv->orbits[temp] == temp) && ((tv->finalnumcells < n) || (OrbSize[tv->orbits[temp]] >= OrbSize[tv->orbits[vertex]])))) {
                tv->fromlevel = tv->tolevel;
                SpineFL = Spine+tv->fromlevel;
                
                if (TargetCell(CurrCand, CurrPart, n, tv, tv->tolevel)) {
                    tv->tcellevel = ++tv->tolevel;
                    SpineTL = Spine+tv->tolevel;
                    SpineTL->tgtcell = tv->tcell;
                    SpineTL->tgtsize = CurrPart->cls[tv->tcell];
                    SpineTL->tgtend = tv->tcell+SpineTL->tgtsize;
                    SpineTL->tgtpos = SpineTL->tgtend - 1;
                }
                else {
                    PRINTF2("CStage2 2: %d\n", tv->finalnumcells);
                    tv->finalnumcells = min(CurrPart->cells,tv->finalnumcells);    /* 160712 */
                    PRINTF2("CStage2 2<: %d\n", tv->finalnumcells);
                    return 0;
                }
                ti->minimalinorbits = TRUE;
                
                if (!ti->identitygroup) {
                    
                    if (ti->useTempOrbits1 && ti->useTempOrbits2) {
                        CuOrb = TempOrbits;
                    }
                    else {
                        FixBase(fix, tv, CurrCand, 0, tv->fromlevel);
                        if (ti->useTempOrbits1 && tv->fromlevel == tv->maxtreelevel) {
                            tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
                            CuOrb = tv->currorbit;
                        }
                        else {
                            if (tv->options->verbosity >= 2) tv->schreier1 -= CPUTIME;
                            
                            gom_level = getorbitsmin(fix, tv->nfix, gpB, &gensB, &tv->currorbit,
                                                     CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell], n, TRUE);
                            if (tv->options->verbosity >= 2) tv->schreier1 += CPUTIME;
                            
                            CuOrb = tv->currorbit;
                            if (gom_level < tv->nfix) {
                                PRINT_NOTMIN_VERB(4)
                                if (ti->useTempOrbits1) {
                                    for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == CuOrb[tv->currorbit[CurrCand->vertex]]) break;
                                    if (i < AutomCount[0]) {
                                        AutomCount[AutomCount[0]++] = CurrCand->vertex;
                                    }
                                    ti->minimalinorbits = FALSE;
                                }
                                else {
                                    TreeNode = CurrCand->stnode;
                                    j2 = CurrCand->lab[Spine[gom_level+1].tgtpos];
                                    i1 = tv->currorbit[j2];
                                    for (j=0; j < tv->nfix - gom_level; j++) {
                                        TreeNode = TreeNode->father;
                                    }
                                    TreeNode1 = TreeNode->first_child;
                                    while (TreeNode1) {
                                        if (TreeNode1->vtx == i1) {
                                            break;
                                        }
                                        TreeNode1 = TreeNode1->next_sibling;
                                    }
                                    schreierwrong = FALSE;
                                    if (TreeNode1) {
                                        while (TreeNode1->goes_to) {
                                            TreeNode1 = TreeNode1->goes_to;
                                        }
                                        TreeNode2 = TreeNode->first_child;
                                        while (TreeNode2->vtx != j2) {
                                            TreeNode2 = TreeNode2->next_sibling;
                                        }
                                        
                                        TreeNode1->index += TreeNode2->index;
                                        TreeNode2->goes_to = TreeNode1;
                                        PRINT_INDEX(TreeNode1,4,28)
                                        PRINT_INDEX(TreeNode2,4,29)
                                        ti->minimalinorbits = FALSE;
                                    }
                                    else {
                                        tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
                                        schreierwrong = TRUE;
                                    }
                                }
                            }
                        }
                    }
                    ti->thegrouphaschanged = FALSE;
                }
                else {
                    CuOrb = IDENTITY_PERM;
                }
                
                if (ti->minimalinorbits) {
                    memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
                    memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
                    SpineTL->trcstart = CurrPart->cells;
                    TheTrace[SpineTL->trcstart] = SpineTL->tgtpos;
                    
                    tv->indivstart = tv->tcell+CurrCand->indnum;
                    tv->indivend = tv->indivstart+tv->steps;
                    if (tv->indivend > SpineTL->tgtend) {
                        tv->indivend = SpineTL->tgtend;
                    }
                    memset(CurrRefCells, 0, n*sizeof(int));
                    ti->thegrouphaschanged = TRUE;
                    
                    if (!CurrCand->sortedlab) {
                        quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
                        for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
                            CurrCand->invlab[CurrCand->lab[i]] = i;
                        }
                        CurrCand->sortedlab = TRUE;
                    }
                    
                    for (k = tv->indivstart; k < tv->indivend; k++) {
                        CurrCand->indnum++;
                        NextCand->singcode = CurrCand->singcode;
                        NextCand->vertex = CurrCand->lab[k];
                        NextCand->name = ++tv->name;
                        if (NextCand->name == (NAUTY_INFINITY-2)) {
                            NextCand->name = tv->name = 1;
                        }
                        
                        if (ti->thegrouphaschanged) {
                            if (tv->fromlevel == tv->maxtreelevel) {
                                CURRORBITSIZES
                            }
                            ti->thegrouphaschanged = FALSE;
                        }
                        
                        if (!schreierwrong) {
                            if (CuOrb[CurrCand->lab[k]] != CurrCand->lab[k]) {
                                continue;
                            }
                        }
                        
                        memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
                        memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
                        if (NextPart->cls[tv->tcell] == 2) {
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell]);
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell+1]);
                        }
                        else {
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[k]);
                        }
                        
                        Individualize(NextPart, NextCand, CurrCand->lab[k], tv->tcell, CurrPart->cells, SpineTL->tgtpos);
                        
                        tv->stats->numnodes++;
                        Spine[tv->tolevel+1].levelcounter++;
                        if (tv->fromlevel == tv->maxtreelevel) {
                            tv->tolevel_tl = tv->tolevel;
                            trieref = trieroot;
                            
                            tv->answ = traces_refine_comptrie(NextCand,
                                                              n,
                                                              NextPart, tv, ti);
                            
                            if (tv->answ) {
                                PRINTF2("CS2 4?: finalnumcells: %d\n", tv->finalnumcells);
                                if (NextPart->cells != tv->finalnumcells) {
                                    CurrRefCells[NextPart->cells % n] += CurrOrbSize[CurrCand->lab[k]];
                                    if (CurrRefCells[NextPart->cells % n] > MultRefCells[NextPart->cells % n]) {
                                        k = n;
                                        break;
                                    }
                                    continue;
                                }
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel);
                            }
                        }
                        else
                        {
                            tv->answ = traces_refine_sametrace(NextCand,
                                                               n,
                                                               NextPart, tv, ti);
                            
                            if (tv->answ) {
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                                    if (tv->tolevel == tv->maxtreelevel) {
                                        tv->tolevel_tl = tv->tolevel;
                                        if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                        if (TargetCellExpPath(NextCand, NextPart, tv)) {
                                            ExperimentalStep(NextPart, NextCand, tv, ti, m, n);
                                            PRINT_EXPPATHSTEP(NextCand, tv->answ)
                                            PRINTF2("CS2 5?: finalnumcells: %d\n", tv->finalnumcells);
                                            if ((NextPart->cells == tv->finalnumcells) || (NextPart->cells == n)) {
                                                UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                            }
                                            
                                        }
                                        if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                        if (!tv->answ) {
                                            PRINT_RETURN
                                        }
                                    }
                            }
                        }
                        if (tv->answ) {
                            PRINTF2("CS2 6?: finalnumcells: %d\n", tv->finalnumcells);
                            if ((NextPart->cells == tv->finalnumcells) || (NextPart->cells == n)) {
                                if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                                temp = (tv->tolevel_tl == tv->tolevel+1);
                                autom = CheckForAutomorphisms(CurrCand, NextCand,
                                                              tv, ti, temp, n, NextPart);
                                if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                                if (autom) {
                                    for (i=autom; i<=tv->maxtreelevel; i++) {
                                        AuxCand = Spine[i].liststart;
                                        while (AuxCand && Prefix(AuxCand, NextCand, autom)) {
                                            AuxCand->do_it = FALSE;
                                            AuxCand = AuxCand->next;
                                        }
                                    }
                                    if (autom == tv->tolevel) {
                                        autom = 0;
                                    }
                                }
                                
                                if (ti->thegrouphaschanged) {
                                    ORBITSIZES
                                }
                            }
                            PRINT_RETURN
                            
                            /* ADD TO NEXT LEVEL */
                            PRINTF2("CS2 7?: finalnumcells: %d\n", tv->finalnumcells);
                            if ((NextPart->cells != tv->finalnumcells) || (tv->tolevel != tv->maxtreelevel) || (tv->tolevel_tl != tv->tolevel+1)) {
                                ADDTONEXTLEVEL;
                                searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                            }
                        }
                        else {
                            tv->stats->interrupted++;
                        }
                        if (autom) {
                            k = n;
                            autom = 0;
                            break;
                        }
                        if (tv->fromlevel == tv->maxtreelevel) {
                            k = n;
                            break;
                        }
                    } /* end for */
                    TreeNode = RefPath[tv->maxtreelevel];
                }
            }
            else SpineTL = &Spine[tv->tolevel+1];
        }
        
    }
    
    /* REMOVE CURRENT CANDIDATE */
    if (!CurrCand->do_it || k >= SpineTL->tgtend) {
        SpineFL->liststart = CurrCand->next;
        if (CurrCand->next == NULL) {
            SpineFL->listend = NULL;
        }
        CurrCand->next = GarbList;
        GarbList = CurrCand;
    }
    ti->thereisnextlevel = SelectNextLevel(n, tv, ti);
    return 0;
}

void CopyCand(Candidate *W, Candidate *V,int n, int *lab, int *invlab) {
    
    if (lab) {
        memcpy(W->lab, lab, n*sizeof(int));
        memcpy(W->invlab, invlab, n*sizeof(int));
    }
    else {
        memcpy(W->lab, V->lab, n*sizeof(int));
        memcpy(W->invlab, V->invlab, n*sizeof(int));
    }
    W->name = V->name;
    W->vertex = V->vertex;
    W->code = V->code;
    W->singcode = V->singcode;
    W->firstsingcode = V->firstsingcode;
    W->do_it = V->do_it;
    W->sortedlab = FALSE;
}

sparsegraph* copy_sg_structure(sparsegraph *sg2, sparsegraph *sg1) {
    int *d1, *e1, *d2, *e2;
    int i, n;
    size_t *v1, *v2, k;
    
    if (!sg2)
    {
        if ((sg2 = (sparsegraph*)ALLOCS(1, sizeof(sparsegraph))) == NULL)
        {
            fprintf(ERRFILE, "copy_sg: malloc failed\n");
            exit(1);
        }
        SG_INIT(*sg2);
    }
    
    SG_VDE(sg1, v1, d1, e1);
    
    n = sg1->nv;
    
    k = 0;
    for (i = 0; i < n; ++i)
        if (v1[i]+d1[i]>k) k = v1[i] + d1[i];
    SG_ALLOC(*sg2, n, k, "copy_sg malloc");
    
    sg2->nv = n;
    sg2->nde = sg1->nde;
    sg2->elen = k;
    /* sg2->wlen = k; */
    return sg2;
}

void Edge_Delete(int vertex, int sons, Candidate *Cand, TracesVars *tv) {
    int d_vtx, j1, temp;
    int *sge, *sgw;
    
    if (TheGraph[vertex].d <= 1) {
        return;
    }
    
    d_vtx = TheGraph[vertex].d = TheGraph[vertex].d - sons;
    sge = TheGraph[vertex].e;
    sgw = TheGraph[vertex].w;
    
    for (j1=0; j1<d_vtx; j1++) {
        if (TheGraph[sge[j1]].one) {
            while (TheGraph[sge[TheGraph[vertex].d]].d == -1) {
                (TheGraph[vertex].d)++;
            }
            temp = sge[j1];
            sge[j1] = sge[TheGraph[vertex].d];
            sge[TheGraph[vertex].d] = temp;
            if (sgw) {
                temp = sgw[j1];
                sgw[j1] = sgw[TheGraph[vertex].d];
                sgw[TheGraph[vertex].d] = temp;
            }
        }
    }
    TheGraph[vertex].d = d_vtx;
}

void ExperimentalStep(Partition *NextPart, Candidate *NextCand,
                      TracesVars *tv, TracesInfo *ti, int m, int n) {
    int i, iend, min, tmp;
    
    SpineTL_tl = Spine+tv->tolevel_tl;
    NextPart->active = 1;
    
    /* EXPERIMENTAL PATH INDIVIDUALIZATION AND REFINEMENT */
    if (tv->answ == 2) {
        min = NextCand->lab[tv->tcellexpath];
        tmp = tv->tcellexpath;
        iend = tv->tcellexpath + NextPart->cls[tv->tcellexpath];
        for (i=tv->tcellexpath + 1; i<iend ; i++) {
            if (NextCand->lab[i] < min) {
                min = NextCand->lab[i];
                tmp = i;
            }
        }
    }
    else {
        tmp = tv->tcellexpath+KRAN(NextPart->cls[tv->tcellexpath]);
    }
    if (NextPart->cls[tv->tcellexpath] == 2) {
        NextCand->pathsingcode = MASHCOMM(NextCand->pathsingcode, NextCand->lab[tv->tcellexpath]);
        NextCand->pathsingcode = MASHCOMM(NextCand->pathsingcode, NextCand->lab[tv->tcellexpath+1]);
    }
    else {
        NextCand->pathsingcode = MASHCOMM(NextCand->pathsingcode, NextCand->lab[tmp]);
    }
    
    tv->indiv_vtx = NextCand->lab[tmp];
    Individualize(NextPart, NextCand, NextCand->lab[tmp], tv->tcellexpath, NextPart->cells, tv->tcellexpath + NextPart->cls[tv->tcellexpath]-1);
    
    tv->stats->numnodes++;
    if (tv->compstage == 0) {
        traces_refine_notrace(NextCand,
                              n,
                              NextPart, tv, ti);
    }
    else {
        if (tv->tolevel_tl == tv->maxtreelevel+1) {
            trieref = trieroot;
            tv->answ = traces_refine_comptrie(NextCand,
                                              n,
                                              NextPart, tv, ti);
            if (tv->answ == 0 ) {
                tv->stats->interrupted++;
            }
        }
        else {
            traces_refine_notrace(NextCand,
                                  n,
                                  NextPart, tv, ti);
        }
    }
    
    CodeClassify(tv->tolevel_tl, NextCand->code, tv->tcellexpath);
    
}

void factorial(double *size1, int *size2, int k) {
    int i;
    
    for(i = k; i; i--) {
        MULTIPLY(*size1, *size2, i);
    }
}

void factorial2(double *size1, int *size2, int k) {
    int i;
    
    for(i = k; i > 0; i -= 2) {
        MULTIPLY(*size1, *size2, i);
    }
}

boolean findperm(permnode *pn, int *p, int n) {
    permnode *rn;
    
    if (!pn) {
        return FALSE;
    }
    rn = pn;
    do {
        if (!memcmp(rn->p, p, n*sizeof(int))) {
            return TRUE;
        }
        rn = rn->next;
    } while (rn != pn);
    return FALSE;
}

int *findcurrorbits(schreier *gp, int k) {
    int i;
    schreier *sh;
    
    sh = gp;
    for (i = 0; i < k; i++) {
        sh = sh->next;
    }
    return sh->orbits;
}

int FirstNeighbour(int vtx, Candidate *Cand, Partition *Part, int* Markers, int mark, int *ngh, int n) {
    int *e_vtx;
    int i, k, deg;
    int ngh1, ngh2, cell1, cell2;
    
    k = 0;
    
    deg = TheGraph[vtx].d;
    e_vtx = TheGraph[vtx].e;
    
    if (deg == n-1) {
        return 0;
    }
    
    for (i=0; i<deg; i++) {
        if (Markers[e_vtx[i]] != mark) {
            cell1 = Part->inv[Cand->invlab[e_vtx[i]]];
            if (Part->cls[cell1] > 1) {
                ngh1 = e_vtx[i++];
                k++;
                break;
            }
        }
    }
    for (; i<deg; i++) {
        if (Markers[e_vtx[i]] != mark) {
            cell2 = Part->inv[Cand->invlab[e_vtx[i]]];
            if (Part->cls[cell2] > 1) {
                ngh2 = e_vtx[i];
                k++;
                break;
            }
        }
    }
    switch (k) {
        case 0:
            break;
            
        case 1:
            *ngh = ngh1;
            break;
            
        case 2:
            if (cell1 < cell2) {
                *ngh = ngh1;
            }
            else {
                *ngh = ngh2;
            }
            break;
            
        default:
            break;
    }
    return k;
}

int FixBase(int *fix, struct TracesVars *tv, Candidate *Cand, int from, int to) {
    int i, j, k, go, nfix;
    
    nfix = j = 0;
    go = TRUE;
    for (i = from; i < to; i++) {
        k = Cand->lab[Spine[i+1].tgtpos];
        if (go && (nfix < tv->nfix) && (fix[nfix] == k)) {
            j++;
        }
        else {
            fix[nfix] = k;
            if (go) go = FALSE;
        }
        nfix++;
    }
    tv->nfix = nfix;
    return j;
}

boolean FixedBase(int *fix, struct TracesVars *tv, Candidate *Cand, int from, int to) {
    int i, k, nfix;
    
    nfix = 0;
    for (i = from; i < to; i++) {
        k = Cand->lab[Spine[i+1].tgtpos];
        if (fix[nfix] != k) {
            return FALSE;
        }
        nfix++;
    }
    return TRUE;
}

int FreeList(Candidate *List, int cond) {
    Candidate *Temp;
    int conta = 0;
    int conta1 = 0;
    
    while (List) {
        if (List->do_it == cond) {
            conta1++;
        }
        conta++;
        Temp = List;
        if (List->lab) free(List->lab);
        if (List->invlab) free(List->invlab);
        List = List->next;
        free(Temp);
    }
    
    if (cond) {
        return conta1;
    }
    else {
        return conta;
    }
}

/* Check if the permutations in the list gens are automorphisms,
 * also set mark and refcount fields and initialise orbits. */
int given_gens(sparsegraph *g, permnode *gens, int *orbits, boolean digraph) {
    int i, m, n, norbs;
    permnode *pn;
    
    n = g->nv;
    for (i = 0; i < n; ++i) orbits[i] = i;
    memcpy(IDENTITY_PERM, orbits, n*sizeof(int));
    norbs = n;
    
    if (!gens) return norbs;
    
    m = SETWORDSNEEDED(n);
    pn = gens;
    do {
        if (!isautom_sg((graph*)g, pn->p, digraph, m, n)) {
            fprintf(ERRFILE, "Input permutation is not an automorphism\n");
            exit(1);
        }
        norbs = orbjoin(orbits, pn->p, n);
        pn->mark = 1;
        pn->refcount = 0;
        pn = pn->next;
    } while (pn != gens);
    
    return norbs;
}

void grouporderplus(sparsegraph *sg_orig, Candidate *Cand, Partition *Part, permnode **ring,
                    double *grpsize1, int *grpsize2, int n, TracesVars *tv, TracesInfo *ti) {
    
    int i, i1, j, j0, j2, k, k1, k2, w, w1, w2, c, c1, c2, n1, n2;
    int prev, step, start, counts, StInd, CyInd, cycnum;
    int tmp, temp, halfsize, nghcell, numvertices;
    int arg, val;
    
    searchtrie *TrieNode;
    int NSFCInd, ind;
    boolean do_ngh = FALSE;
    
    numvertices = n;
    memcpy(CanonIndices, IDENTITY_PERM, n*sizeof(int));
    memset(TreeNodes, 0, n*sizeof(int));
    
    TrieNode = Spine[tv->maxtreelevel].liststart->stnode;
    if (TrieNode->father) {
        if (tv->options->verbosity >= 2) {
            LINE(32, "—")
            NEXTLINE
            fprintf(outfile, "group structure: ");
            while (TrieNode->father) {
                if (Factorials[TrieNode->level]) {
                    fprintf(outfile, "%d (%d!), ", TrieNode->name, TrieNode->index);
                    factorial(grpsize1, grpsize2, TrieNode->index);
                } else {
                    if (TrieNode->father->name) {
                        fprintf(outfile, "%d (%d), ", TrieNode->name, TrieNode->index);
                        MULTIPLY(tv->stats->grpsize1, tv->stats->grpsize2, TrieNode->index);
                    }
                    else {
                        temp = spinelementorbsize(tv->orbits, Spine[tv->maxtreelevel].liststart->lab+Spine[1].tgtcell, Spine[1].tgtsize, TrieNode->vtx);
                        MULTIPLY(*grpsize1, *grpsize2, temp);
                        fprintf(outfile, "%d (%d)\n",
                                TrieNode->name, temp);
                    }
                }
                TrieNode = TrieNode->father;
            }
        }
        else {
            while (TrieNode->father) {
                if (Factorials[TrieNode->level]) {
                    factorial(grpsize1, grpsize2, TrieNode->index);
                } else {
                    if (TrieNode->father->name) {
                        MULTIPLY(tv->stats->grpsize1, tv->stats->grpsize2, TrieNode->index);
                    }
                    else {
                        temp = spinelementorbsize(tv->orbits, Spine[tv->maxtreelevel].liststart->lab+Spine[1].tgtcell, Spine[1].tgtsize, TrieNode->vtx);
                        MULTIPLY(*grpsize1, *grpsize2, temp);
                    }
                }
                TrieNode = TrieNode->father;
            }
        }
    }
    
    if (Part->cells < n) {
        
        if (!ti->deg_one) {
            memcpy(tv->graph->e, sg_orig->e, tv->graph->elen*sizeof(int));
            for (i=0; i<n; i++) {
                TheGraph[i].e = tv->graph->e + sg_orig->v[i];
            }
        }
        NSFCInd = 0;
        
        /* Trees */
        if (tv->options->getcanon && tv->preprocessed) {
            for (i = 0; i < n; i += Part->cls[i]) {
                if (Part->cls[i] == 1) {
                    tmp = Cand->lab[i];
                    if ((TheGraph[tmp].d >= 0) && (TheGraph[tmp].d < sg_orig->d[tmp])) {
                        for (j=i; j<i+Part->cls[i]; j++) {
                            MakeCanTree(Cand->lab[j], sg_orig, n, Cand, Part, tv);
                        }
                    }
                }
            }
        }
        
        memset(SingNonSing, 0, n*sizeof(int));
        
        for (i = 0; i < n; i += Part->cls[i]) {
            if (Part->cls[i] > 1) {
                if (TheGraph[Cand->lab[i]].d > 2) {
                    for (j=i; j<i+Part->cls[i]; j++) {
                        SingNonSing[Cand->lab[j]] = 2;
                    }
                }
            }
            else {
                SingNonSing[Cand->lab[i]] = 1;
            }
        }
        
        for (i = 0; i < n; i += Part->cls[i]) {
            if (Part->cls[i] > 1) {
                if (TheGraph[Cand->lab[i]].d > 2) NonSingDegPlus1(Cand, Part, i, tv);
                NSFCells[NSFCInd++] = i;
            }
            else {
                NonSingDegPlus2(Cand, Part, i, tv);
                numvertices--;
            }
        }
        
        /* Degree 2 and at least one nghb with deg > 2 */
        SETMARK(StackMarkers, tv->stackmark)
        for (ind = 0; ind < NSFCInd; ind++) {
            i = NSFCells[ind];
            SETMARK(Markers, tv->mark)
            if (Part->cls[i] > 1) {
                tmp = Cand->lab[i];
                if ((TheGraph[tmp].d == 2) && ((TheGraph[TheGraph[tmp].e[0]].d > 2) || ((TheGraph[TheGraph[tmp].e[1]].d > 2)))) {
                    n1 = TheGraph[tmp].e[0];
                    n2 = TheGraph[tmp].e[1];
                    if (TheGraph[n1].d > 2) {
                        if (TheGraph[n2].d > 2) {
                            if (Cand->invlab[n1] < Cand->invlab[n2]) {
                                start = n1;
                            }
                            else {
                                start = n2;
                            }
                        }
                        else {
                            start = n1;
                        }
                    }
                    else {
                        start = n2;
                    }
                    counts = 0;
                    StInd = 0;
                    for (j=i; j<i+Part->cls[i]; j++) {
                        step = Cand->lab[j];
                        if (Markers[step] != tv->mark) {
                            prev = start;
                            counts++;
                            do {
                                Markers[step] = tv->mark;
                                PERMSTACK[StInd++] = step;
                                if (TheGraph[step].e[0] != prev) {
                                    prev = step;
                                    step = TheGraph[step].e[0];
                                }
                                else {
                                    prev = step;
                                    step = TheGraph[step].e[1];
                                }
                            } while (TheGraph[step].d == 2);
                            
                            if (TheGraph[step].d == 1) {
                                PERMSTACK[StInd++] = step;
                            }
                        }
                    }
                    
                    if (counts == Part->cls[i]) {
                        factorial(grpsize1, grpsize2, Part->cls[i]);
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        for (k=0; k<StInd/counts; k++) {
                            i1 = PERMSTACK[k];
                            for (j0=0; j0<counts-1; j0++) {
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                            }
                            SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                        }
                        SPECIALGENERATORS
                        if (counts > 2) {
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                        }
                    }
                    else {
                        factorial2(grpsize1, grpsize2, Part->cls[i]);
                        for (j=0; j<counts; j++) {
                            j0 = j*(StInd/counts);
                            k1 = (j+1)*(StInd/counts);
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=j0, i1=k1-1; k<k1; k++, i1--) {
                                SETPAIRSAUTANDTREE(PERMSTACK[k], PERMSTACK[i1])
                            }
                            SPECIALGENERATORS
                        }
                        if (counts > 1) {
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<counts-1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                        }
                        if (counts > 2) {
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                        }
                    }
                    for (j=0; j<StInd; j++) {
                        Place(PERMSTACK[j], Cand, Part);
                        if ((TheGraph[PERMSTACK[j]].d >= 0) && (TheGraph[PERMSTACK[j]].d < sg_orig->d[PERMSTACK[j]])) {
                            MakeCanTree(PERMSTACK[j], sg_orig, n, Cand, Part, tv);
                        }
                    }
                }
            }
        }
        
        /* Degree 2 and at least one nghb with == 1 */
        for (ind = 0; ind < NSFCInd; ind++) {
            SETMARK(Markers, tv->mark)
            i = NSFCells[ind];
            if (Part->cls[i] > 1) {
                tmp = Cand->lab[i];
                if ((TheGraph[tmp].d == 2) && ((TheGraph[TheGraph[tmp].e[0]].d == 1) || ((TheGraph[TheGraph[tmp].e[1]].d == 1)))) {
                    counts = 0;
                    StInd = 0;
                    for (j=i; j<i+Part->cls[i]; j++) {
                        step = Cand->lab[j];
                        if (Markers[step] != tv->mark) {
                            n1 = TheGraph[step].e[0];
                            n2 = TheGraph[step].e[1];
                            if (TheGraph[n1].d == 1) {
                                if (TheGraph[n2].d == 1) {
                                    if (Cand->invlab[n1] < Cand->invlab[n2]) {
                                        start = n1;
                                    }
                                    else {
                                        start = n2;
                                    }
                                }
                                else {
                                    start = n1;
                                }
                            }
                            else {
                                start = n2;
                            }
                            PERMSTACK[StInd++] = start;
                            prev = start;
                            counts++;
                            
                            do {
                                Markers[step] = tv->mark;
                                PERMSTACK[StInd++] = step;
                                if (TheGraph[step].e[0] != prev) {
                                    prev = step;
                                    step = TheGraph[step].e[0];
                                } else {
                                    prev = step;
                                    step = TheGraph[step].e[1];
                                }
                            } while (TheGraph[step].d == 2);
                            PERMSTACK[StInd++] = step;
                        }
                    }
                    
                    if (counts == Part->cls[i]) {
                        if (Part->inv[Cand->invlab[PERMSTACK[0]]] != Part->inv[Cand->invlab[PERMSTACK[StInd/counts-1]]]) {
                            factorial(grpsize1, grpsize2, Part->cls[i]);
                        }
                        else {
                            factorial2(grpsize1, grpsize2, 2*Part->cls[i]);
                            for (j=0; j<counts; j++) {
                                j0 = j*(StInd/counts);
                                k1 = (j+1)*(StInd/counts);
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                for (k=j0, i1=k1-1; k<k1; k++, i1--) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[k], PERMSTACK[i1])
                                }
                                SPECIALGENERATORS
                            }
                        }
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        for (k=0; k<StInd/counts; k++) {
                            i1 = PERMSTACK[k];
                            for (j0=0; j0<counts-1; j0++) {
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                            }
                            SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                        }
                        SPECIALGENERATORS
                        if (counts > 2) {
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                        }
                    }
                    else {
                        factorial2(grpsize1, grpsize2, Part->cls[i]);
                        for (j=0; j<counts; j++) {
                            j0 = j*(StInd/counts);
                            k1 = (j+1)*(StInd/counts);
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=j0, i1=k1-1; k<k1; k++, i1--) {
                                SETPAIRSAUTANDTREE(PERMSTACK[k], PERMSTACK[i1])
                            }
                            SPECIALGENERATORS
                        }
                        if (counts > 1) {
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<counts-1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                        }
                        if (counts > 2) {
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                        }
                    }
                    for (j=0; j<StInd; j++) {
                        Place(PERMSTACK[j], Cand, Part);
                        if ((TheGraph[PERMSTACK[j]].d >= 0) && (TheGraph[PERMSTACK[j]].d < sg_orig->d[PERMSTACK[j]])) {
                            MakeCanTree(PERMSTACK[j], sg_orig, n, Cand, Part, tv);
                        }
                    }
                }
            }
        }
        
        /* Cycles */
        for (ind = 0; ind < NSFCInd; ind++) {
            i = NSFCells[ind];
            SETMARK(Markers, tv->mark)
            if (Part->cls[i] > 1) {
                tmp = Cand->lab[i];
                if (TheGraph[tmp].d == 2) {
                    CyInd = StInd = cycnum = 0;
                    for (j=i; j<i+Part->cls[i]; j++) {
                        start = Cand->lab[j];
                        if (Markers[start] != tv->mark) {
                            counts = 1;
                            CYCLES[StInd] = start;
                            CYCOLR[StInd++] = Part->inv[Cand->invlab[start]];
                            Markers[start] = tv->mark;
                            k = Cand->invlab[TheGraph[start].e[0]];
                            k1 = Cand->invlab[TheGraph[start].e[1]];
                            if (Part->inv[k] < Part->inv[k1]) {
                                step = TheGraph[start].e[0];
                            }
                            else {
                                step = TheGraph[start].e[1];
                            }
                            prev = start;
                            do {
                                counts++;
                                Markers[step] = tv->mark;
                                CYCLES[StInd] = step;
                                CYCOLR[StInd++] = Part->inv[Cand->invlab[step]];
                                
                                if (TheGraph[step].e[0] != prev) {
                                    prev = step;
                                    step = TheGraph[step].e[0];
                                }
                                else {
                                    prev = step;
                                    step = TheGraph[step].e[1];
                                }
                            } while (step != start);
                            CYLGTH[CyInd++] = counts;
                            cycnum++;
                        }
                    }
                    
                    CYCPOS[0] = 0;
                    for (j=1; j<CyInd; j++) {
                        CYCPOS[j] = CYCPOS[j-1]+CYLGTH[j-1];
                    }
                    memcpy(WorkArray, CYLGTH, CyInd*sizeof(int));
                    sort2ints(WorkArray, CYCPOS, CyInd);
                    
                    k = 0;
                    for (i1=0; i1<CyInd; i1++) {
                        k1 = CYCOLR[k];
                        k2 = CYCOLR[k+1];
                        for (j=1; j<=CYLGTH[i1]/2; j++) {
                            w1 = CYCOLR[j+k];
                            w2 = CYCOLR[j+1+k];
                            if ((w1 == k1) && (w2 == k2)) {
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                for (w=0; w<CYLGTH[i1]; w++) {
                                    if (CYCOLR[w+k] == CYCOLR[((w+j) % CYLGTH[i1]) + k]) {
                                        SETPAIRSAUTANDTREE(CYCLES[w+k], CYCLES[((w+j) % CYLGTH[i1]) + k])
                                    }
                                    else {
                                        break;
                                    }
                                }
                                if (w == CYLGTH[i1]) { SPECIALGENERATORS }
                                if (w == CYLGTH[i1]) {
                                    MULTIPLY(*grpsize1, *grpsize2, CYLGTH[i1]/j);
                                    break;
                                }
                            }
                        }
                        
                        if (Part->cls[k1] >= Part->cls[k2]) {
                            for (j=CYLGTH[i1]-1; j>0; j--) {
                                w1 = CYCOLR[j % CYLGTH[i1] + k];
                                w2 = CYCOLR[(j-1) % CYLGTH[i1] + k];
                                if ((w1 == k1) && (w2 == k2)) {
                                    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                    for (w=0; w<CYLGTH[i1]; w++) {
                                        SETPAIRSAUTANDTREE(CYCLES[w+k], CYCLES[((j-w+(w>j)*CYLGTH[i1]) % CYLGTH[i1]) + k])
                                    }
                                    SPECIALGENERATORS
                                    MULTIPLY(*grpsize1, *grpsize2, 2);
                                    break;
                                }
                            }
                        }
                        else {
                            j=CYLGTH[i1]-1;
                            w2 = CYCOLR[j % CYLGTH[i1] + k];
                            if (w2 == k2) {
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                for (w=1; w<CYLGTH[i1]; w++) {
                                    SETPAIRSAUTANDTREE(CYCLES[w+k], CYCLES[CYLGTH[i1]-w+k])
                                }
                                SPECIALGENERATORS
                                MULTIPLY(*grpsize1, *grpsize2, 2);
                            }
                        }
                        k += CYLGTH[i1];
                    }
                    k = 0;
                    for (i1=0; i1<CyInd; i1++) {
                        if (CYLGTH[i1] > 0) {
                            CYMULT[0] = k;
                            k1 = k;
                            counts = 1;
                            for (j0=i1+1; j0<CyInd; j0++) {
                                k1 += abs(CYLGTH[j0]);
                                if (CYLGTH[j0] == CYLGTH[i1]) {
                                    CYMULT[counts++] = k1;
                                    CYLGTH[j0] = -CYLGTH[j0];
                                }
                            }
                            if (counts > 1) {
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                for (j0=0; j0<CYLGTH[i1]; j0++) {
                                    for (j2 = 0; j2<counts-1; j2++) {
                                        SETPAIRSAUTANDTREE(CYCLES[CYMULT[j2]+j0], CYCLES[CYMULT[j2+1]+j0])
                                    }
                                    SETPAIRSAUTANDTREE(CYCLES[CYMULT[j2]+j0], CYCLES[CYMULT[0]+j0])
                                }
                                SPECIALGENERATORS
                                if (counts > 2) {
                                    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                    for (j0=0; j0<CYLGTH[i1]; j0++) {
                                        SETPAIRSAUTANDTREE(CYCLES[CYMULT[1]+j0], CYCLES[CYMULT[0]+j0])
                                        if (tv->build_autom) {
                                            SETPAIRSAUT(CYCLES[CYMULT[0]+j0], CYCLES[CYMULT[1]+j0])
                                        }
                                        MakeTree(CYCLES[CYMULT[0]+j0], CYCLES[CYMULT[1]+j0], sg_orig, n, tv, FALSE);
                                    }
                                    SPECIALGENERATORS
                                }
                                factorial(grpsize1, grpsize2, counts);
                            }
                        }
                        k += abs(CYLGTH[i1]);
                        CYLGTH[i1] = -CYLGTH[i1];
                    }
                    
                    for (c1=0; c1<CyInd; c1++) {
                        c = CYCPOS[c1]+WorkArray[c1];
                        for (c2=CYCPOS[c1]; c2<c; c2++) {
                            Place(CYCLES[c2], Cand, Part);
                            if ((TheGraph[CYCLES[c2]].d >= 0) && (TheGraph[CYCLES[c2]].d < sg_orig->d[CYCLES[c2]])) {
                                MakeCanTree(CYCLES[c2], sg_orig, n, Cand, Part, tv);
                            }
                        }
                    }
                }
            }
        }
        
        /* Degree 1, and nghb too */
        SETMARK(Markers, tv->mark)
        for (ind = 0; ind < NSFCInd; ind++) {
            i = NSFCells[ind];
            if (Part->cls[i] > 1) {
                tmp = Cand->lab[i];
                if ((TheGraph[tmp].d == 1) && (TheGraph[TheGraph[tmp].e[0]].d == 1) && (i == Part->inv[Cand->invlab[TheGraph[tmp].e[0]]])) {
                    factorial2(grpsize1, grpsize2, Part->cls[i]);
                    /* the cell has size two */
                    if (Part->cls[i] == 2) {
                        val = Cand->lab[i+1];
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        arg = tmp;
                        SETPAIRSAUTANDTREE(arg, val)
                        SETPAIRSAUTANDTREE(val, arg)
                        SPECIALGENERATORS
                    }
                    else {
                        /* the cell has size greater than two */
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        SETMARK(Markers, tv->mark)
                        halfsize = Part->cls[i]/2;
                        i1 = 0;
                        for (j=i; j<i+Part->cls[i]; j++) {
                            if (Markers[Cand->lab[j]] != tv->mark) {
                                Markers[TheGraph[Cand->lab[j]].e[0]] = tv->mark;
                                PERMSTACK[i1] = Cand->lab[j];
                                PERMSTACK[i1+halfsize] = TheGraph[Cand->lab[j]].e[0];
                                i1++;
                            }
                        }
                        temp = PERMSTACK[0];
                        for (j=0; j<Part->cls[i]-1; j++) {
                            SETPAIRSAUTANDTREE(PERMSTACK[j], PERMSTACK[j+1])
                        }
                        SETPAIRSAUTANDTREE(PERMSTACK[j], temp)
                        SPECIALGENERATORS
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        memmove(PERMSTACK+halfsize, PERMSTACK+halfsize+1, (halfsize-1)*sizeof(int));
                        temp = PERMSTACK[1];
                        for (j=1; j<Part->cls[i]-2; j++) {
                            SETPAIRSAUTANDTREE(PERMSTACK[j], PERMSTACK[j+1])
                        }
                        SETPAIRSAUTANDTREE(PERMSTACK[j], temp)
                        SPECIALGENERATORS
                    }
                    
                    SETMARK(Markers, tv->mark)
                    for (j=i; j<i+Part->cls[i]; j++) {
                        temp = Cand->lab[j];
                        if (Markers[temp] != tv->mark) {
                            if ((TheGraph[temp].d >= 0) && (TheGraph[temp].d < sg_orig->d[temp])) {
                                MakeCanTree(temp, sg_orig, n, Cand, Part, tv);
                            }
                            tmp = Cand->lab[j+1];
                            Markers[TheGraph[temp].e[0]] = tv->mark;
                            i1 = Cand->invlab[TheGraph[temp].e[0]];
                            Cand->lab[j+1] = TheGraph[temp].e[0];
                            
                            if ((TheGraph[TheGraph[temp].e[0]].d >= 0) && (TheGraph[TheGraph[temp].e[0]].d < sg_orig->d[TheGraph[temp].e[0]])) {
                                MakeCanTree(TheGraph[temp].e[0], sg_orig, n, Cand, Part, tv);
                            }
                            Cand->invlab[TheGraph[temp].e[0]] = j+1;
                            Cand->lab[i1] = tmp;
                            Cand->invlab[tmp] = i1;
                        }
                    }
                }
            }
        }
        
        /* Degree 0 */
        for (ind = 0; ind < NSFCInd; ind++) {
            i = NSFCells[ind];
            if (Part->cls[i] > 1) {
                tmp = Cand->lab[i];
                if ((TheGraph[0].e != NULL) && (TheGraph[tmp].d != 0) && (TheGraph[tmp].d != numvertices-1))
                    nghcell = Part->inv[Cand->invlab[TheGraph[tmp].e[0]]]; else nghcell = i;
                if ((TheGraph[tmp].d == 0) ||
                    ((TheGraph[tmp].d == numvertices-1) && (TheGraph[tmp].d > 2)) ||
                    ((TheGraph[tmp].d == 1) && (TheGraph[TheGraph[tmp].e[0]].d == 1) && (i < nghcell))) {
                    do_ngh = FALSE;
                    if ((TheGraph[tmp].d == 1) && (TheGraph[TheGraph[tmp].e[0]].d == 1) && (i != nghcell)) {
                        do_ngh = TRUE;
                    }
                    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                    for (j=i; j<i+Part->cls[i]-1; j++) {
                        arg = Cand->lab[j];
                        val = Cand->lab[j+1];
                        SETPAIRSAUTANDTREE(arg, val)
                        if (do_ngh) {
                            SETPAIRSAUTANDTREE(TheGraph[arg].e[0], TheGraph[val].e[0])
                            
                        }
                    }
                    arg = Cand->lab[j];
                    val = tmp;
                    SETPAIRSAUTANDTREE(arg, val)
                    if (do_ngh) {
                        SETPAIRSAUTANDTREE(TheGraph[arg].e[0], TheGraph[val].e[0])
                    }
                    SPECIALGENERATORS
                    if (Part->cls[i] > 2) {
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        arg = tmp;
                        val = Cand->lab[i+1];
                        SETPAIRSAUTANDTREE(arg, val)
                        if (do_ngh) {
                            SETPAIRSAUTANDTREE(TheGraph[arg].e[0], TheGraph[val].e[0])
                        }
                        arg = Cand->lab[i+1];
                        val = tmp;
                        SETPAIRSAUTANDTREE(arg, val)
                        if (do_ngh) {
                            SETPAIRSAUTANDTREE(TheGraph[arg].e[0], TheGraph[val].e[0])
                        }
                        SPECIALGENERATORS
                    }
                    factorial(grpsize1, grpsize2, Part->cls[i]);
                    if (do_ngh) {
                        for (j=i; j<i+Part->cls[i]; j++) {
                            temp = TheGraph[Cand->lab[j]].e[0];
                            Cand->lab[nghcell] = temp;
                            Cand->invlab[temp] = nghcell;
                            nghcell++;
                        }
                    }
                    
                    k = i+Part->cls[i];
                    for (j=i; j<k; j++) {
                        Place(Cand->lab[j], Cand, Part);
                        if ((TheGraph[Cand->lab[j]].d >= 0) && (TheGraph[Cand->lab[j]].d < sg_orig->d[Cand->lab[j]])) {
                            MakeCanTree(Cand->lab[j], sg_orig, n, Cand, Part, tv);
                            if (do_ngh) {
                                MakeCanTree(TheGraph[Cand->lab[j]].e[0], sg_orig, n, Cand, Part, tv);
                            }
                        }
                    }
                }
            }
        }
        
    }
    
    /* Orbit Count */
    SETMARK(Markers, tv->mark)
    i1=0;
    for (c1=0; c1<n; c1++) {
        if (Markers[tv->orbits[c1]] != tv->mark) {
            i1++;
            Markers[tv->orbits[c1]] = tv->mark;
        }
    }
    tv->stats->numorbits = i1;
    return;
}

void Individualize(Partition *NextPart, Candidate *NextCand, int K, int Tc, int Cl, int Pos) {
    int i, j;
    
    NextCand->do_it = TRUE;
    if (NextPart->cls[Tc] > 1) {
        NextPart->cells = Cl+1;
        NextPart->active = 1;
        NextPart->cls[Tc]--;
        NextPart->cls[Pos] = 1;
    }
    NextPart->inv[Pos] = Pos;
    
    j = NextCand->lab[Pos];
    i = NextCand->invlab[K];
    NextCand->lab[Pos] = K;
    NextCand->invlab[K] = Pos;
    NextCand->lab[i] = j;
    NextCand->invlab[j] = i;
    return;
}

void Initialize_Traces_Variables(TracesVars *tv, TracesOptions *options_arg,
                                 TracesStats *stats_arg, int *orbits_arg,
                                 sparsegraph *g_arg, sparsegraph *canong_arg,
                                 int n) {
    tv->augmented_cells = n;
    tv->canlist = 0;
    tv->compstage = 0;
    tv->expathlength = n;
    tv->finalnumcells = n;
    tv->firstpathlength = 0;
    tv->group_level = 0;
    tv->lev_of_lastauto = 0;
    tv->linelgth = 0;
    tv->name = 0;
    tv->maxspineorblevel = 0;
    tv->nfix = 0;
    tv->options = options_arg;
    tv->orbits = orbits_arg;
    tv->permInd = 0;
    tv->maxdeg = 0;
    tv->mindeg = n;
    
    if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc)
        tv->build_autom = TRUE;
    else
        tv->build_autom = FALSE;
    
    tv->specialgens = 0;
    tv->stats = stats_arg;
    tv->treedepth = 0;
    tv->gotonode = NULL;
    tv->input_graph = tv->graph = g_arg;
    tv->cangraph = canong_arg;
    tv->mark = tv->stackmark = tv->treemark = tv->autmark = tv->markcell1 = tv->markcell2 = NAUTY_INFINITY-1;
    tv->conta0 = tv->conta1 = tv->conta2 = tv->conta3 = tv->conta4 = tv->conta5 = tv->conta6 = tv->conta7 = tv->contatc = 0;
    
    if (tv->options->strategy == 0) {
        tv->steps = n;
        tv->strategy = 0;
    }
    else {
        tv->strategy = 1;
        tv->steps = tv->options->strategy;
        if (tv->steps > n) {
            tv->steps = n;
        }
    }
}

void Initialize_Traces_Statistics (TracesStats *stats_arg, int n) {
    stats_arg->grpsize1 = 1;
    stats_arg->grpsize2 = 0;
    stats_arg->numorbits = n;
    stats_arg->treedepth= 0;
    stats_arg->numgenerators = 0;
    stats_arg->numnodes = 1;
    stats_arg->errstatus = 0;
    stats_arg->interrupted = 0;
    stats_arg->canupdates = 0;
    stats_arg->peaknodes = 0;
}

void Initialize_Traces_Time_Variables (TracesVars *tv) {
    tv->autchk = 0;
    tv->expaths = 0;
    tv->schreier1 = 0;
    tv->schreier2 = 0;
    tv->schreier3 = 0;
}

boolean isautom_sg_pair(graph *g, int *p, boolean digraph, int m, int n, struct TracesVars *tv) {
    int *d, *e;
    size_t *v;
    int i, k, pi, di;
    size_t vi, vpi, j;
    
    SG_VDE(g, v, d, e);
    
    for (k = 0; k < tv->permInd; ++k)
    {
        i = PrmPairs[k].arg;
        pi = p[i];
        di = d[i];
        if (d[pi] != di) return FALSE;
        
        vi = v[i];
        vpi = v[pi];
        SETMARK(AutMarkers, tv->autmark)
        for (j = 0; j < di; ++j) AutMarkers[p[e[vi+j]]] = tv->autmark;
        for (j = 0; j < di; ++j) if (AutMarkers[e[vpi+j]] != tv->autmark) {
            return FALSE;
        }
    }
    
    return TRUE;
}

boolean lookup(searchtrie *t) {
    searchtrie *TreeNode;
    
    TreeNode = t;
    while (TreeNode->level >= 1) {
        if (TreeNode->goes_to) {
            return FALSE;
        }
        TreeNode = TreeNode->father;
    }
    return TRUE;
}

void MakeCanTree(int v1, sparsegraph *sg_orig, int n, Candidate *Cand, Partition *Part, struct TracesVars* tv) {
    int ind, vtx, ngh, trind, deg0, deg1;
    size_t j1;
    int *sge1;
    
    trind = 1;
    ind = 0;
    TreeStack[0] = v1;
    SETMARK(TreeMarkers, tv->treemark);
    
    while (ind < trind) {
        vtx = TreeStack[ind++];
        if (TreeNodes[vtx]) {
            return;
        }
        
        if (TheGraph[vtx].d == -1) {
            Place(vtx, Cand, Part);
            TreeNodes[vtx] = TRUE;
        }
        
        TreeMarkers[vtx] = tv->treemark;
        deg0 = max(TheGraph[vtx].d, 0);
        deg1 = sg_orig->d[vtx];
        sge1 = TheGraph[vtx].e;
        
        for (j1 = deg0; j1 < deg1; j1++) {
            ngh = sge1[j1];
            if ((TheGraph[ngh].d == -1) && (TreeMarkers[ngh] != tv->treemark)) {
                TreeStack[trind++] = ngh;
            }
        }
    }
    return;
}

void MakeDiscrete(Partition *Part, int cell) {
    int i, k;
    
    Part->cells += (Part->cls[cell] - 1);
    k = cell + Part->cls[cell];
    
    for (i = cell; i < k; i++) {
        Part->cls[i] = 1;
        Part->inv[i] = i;
    }
}

void MakeTree(int v1, int v2, sparsegraph *sg, int n, struct TracesVars* tv, boolean forceautom) {
    int ind, vtx1, vtx2, ngh1, ngh2, trind, deg0, deg1;
    size_t j1;
    int *sge1, *sge2;
    boolean build_autom;
    
    if (v1 == v2) return;
    build_autom = tv->build_autom || forceautom;
    trind = 2;
    ind = 0;
    TreeStack[0] = v1;
    TreeStack[1] = v2;
    SETMARK(TreeMarkers, tv->treemark);
    
    while (ind < trind) {
        vtx1 = TreeStack[ind++];
        vtx2 = TreeStack[ind++];
        
        TreeMarkers[vtx1] = tv->treemark;
        TreeMarkers[vtx2] = tv->treemark;
        
        deg0 = max(TheGraph[vtx1].d, 0);
        deg1 = sg->d[vtx1];
        sge1 = TheGraph[vtx1].e;
        sge2 = TheGraph[vtx2].e;
        for (j1 = deg0; j1 < deg1; j1++) {
            ngh1 = sge1[j1];
            ngh2 = sge2[j1];
            if ((TreeMarkers[ngh1] != tv->treemark) && (ngh1 != ngh2)) {
                TreeStack[trind++] = ngh1;
                TreeStack[trind++] = ngh2;
                if (ngh1 != ngh2) {
                    if (build_autom) {
                        AUTPERM[ngh1] = ngh2;
                        PrmPairs[tv->permInd].arg = ngh1;
                        PrmPairs[tv->permInd].val = ngh2;
                        tv->permInd++;
                    }
                    orbjoin_sp_pair(tv->orbits, OrbList, n,
                                    ngh1, ngh2, &tv->stats->numorbits);
                }
            }
        }
    }
    return;
}

int max(int u, int v) {
    if (u > v) {
        return u;
    }
    else {
        return v;
    }
}

int min(int u, int v) {
    if (u < v) {
        return u;
    }
    else {
        return v;
    }
}

int NextNeighbour(int vtx, Candidate *Cand, Partition *Part, int* Markers, int mark, int *ngh, int n) {
    int *e_vtx;
    int i, deg;
    int cell1;
    
    deg = TheGraph[vtx].d;
    e_vtx = TheGraph[vtx].e;
    
    if (deg == n-1) {
        return 0;
    }
    
    for (i=0; i<deg; i++) {
        if (Markers[e_vtx[i]] != mark) {
            cell1 = Part->inv[Cand->invlab[e_vtx[i]]];
            if (Part->cls[cell1] > 1) {
                *ngh = e_vtx[i];
                break;
            }
        }
    }
    if (i<deg) return 1; else return 0;
}

int NonSingDeg(int vtx, Candidate *Cand, Partition *Part) {
    int *e_vtx;
    int i, deg, retdeg;
    
    retdeg = TheGraph[vtx].d;
    deg = retdeg;
    e_vtx = TheGraph[vtx].e;
    for (i=0; i<deg; i++) {
        if (Part->cls[Part->inv[Cand->invlab[e_vtx[i]]]] == 1) {
            retdeg--;
        }
    }
    return retdeg;
}

int NonSingDegPlus1(Candidate *Cand, Partition *Part, int cell, TracesVars *tv) {
    
    int *e_vtx;
    int vtx, sing;
    int i, j, deg, retdeg, n, singcount;
    
    n = tv->input_graph->nv;
    singcount = 0;
    
    SETMARK(StackMarkers, tv->stackmark)
    
    for (j=cell; j<cell+Part->cls[cell]; j++) {
        vtx = Cand->lab[j];
        deg = TheGraph[vtx].d;
        retdeg = 0;
        e_vtx = TheGraph[vtx].e;
        
        for (i=0; i<deg; i++) {
            if (SingNonSing[e_vtx[i]] != 1) {
                e_vtx[retdeg++] = e_vtx[i];
            }
            else {
                if (StackMarkers[e_vtx[i]] != tv->stackmark) {
                    sing = e_vtx[i];
                    WorkArray2[singcount] = Part->inv[Cand->invlab[sing]];
                    WorkArray[singcount++] = sing;
                    
                    StackMarkers[e_vtx[i]] = tv->stackmark;
                }
            }
        }
        if (j == cell) {
            sort2ints(WorkArray2, WorkArray, singcount);
        }
        if (deg != retdeg) {
            memcpy(e_vtx+retdeg, WorkArray, singcount*sizeof(int));
            TheGraph[vtx].d = retdeg;
        }
    }
    return retdeg;
}

void NonSingDegPlus2(Candidate *Cand, Partition *Part, int cell, TracesVars *tv) {
    
    int *e_sing;
    int sing;
    int k, deg1, singdeg, singcount;
    
    singcount = 0;
    
    sing = Cand->lab[cell];
    singdeg = 0;
    deg1 = TheGraph[sing].d;
    e_sing = TheGraph[sing].e;
    
    for (k=0; k<deg1; k++) {
        if (SingNonSing[e_sing[k]] != 2) {
            e_sing[singdeg++] = e_sing[k];
        }
    }
    TheGraph[sing].d = singdeg;
}

void orbjoin_sp_pair(int *orbits, int *list, int n, int u, int v, int *numorbs) {
    int j1, j2, k1, k2;
    
    j1 = orbits[u];
    while (orbits[j1] != j1) j1 = orbits[j1];
    j2 = orbits[v];
    while (orbits[j2] != j2) j2 = orbits[j2];
    
    if (j1 != j2) {
        k1 = j1;
        k2 = j2;
        if (k1 < k2) {
            (*numorbs)--;
            while (list[j2] != k2) {
                orbits[j2] = k1;
                j2 = list[j2];
            }
            orbits[j2] = k1;
            k1 = list[k1];
            list[j2] = k1;
            list[j1] = k2;
        }
        else if (k1 > k2) {
            (*numorbs)--;
            while (list[j1] != k1) {
                orbits[j1] = k2;
                j1 = list[j1];
            }
            orbits[j1] = k2;
            k2 = list[k2];
            list[j1] = k2;
            list[j2] = k1;
        }
    }
    return;
}

void orbjoin_sp_perm(int *orbits, int *map, int *list, int n, int *numorbs) {
    int i, j1, j2, k1, k2;
    
    for (i = 0; i < n; ++i)
        if (map[i] != i)
        {
            j1 = orbits[i];
            while (orbits[j1] != j1) j1 = orbits[j1];
            j2 = orbits[map[i]];
            while (orbits[j2] != j2) j2 = orbits[j2];
            k1 = j1;
            k2 = j2;
            if (k1 < k2) {
                (*numorbs)--;
                while (OrbList[j2] != k2) {
                    orbits[j2] = k1;
                    j2 = OrbList[j2];
                }
                orbits[j2] = k1;
                k1 = OrbList[k1];
                OrbList[j2] = k1;
                OrbList[j1] = k2;
            }
            else if (k1 > k2) {
                (*numorbs)--;
                while (OrbList[j1] != k1) {
                    orbits[j1] = k2;
                    j1 = OrbList[j1];
                }
                orbits[j1] = k2;
                k2 = OrbList[k2];
                OrbList[j1] = k2;
                OrbList[j2] = k1;
            }
        }
}

struct Partition *NewPartition(int n) {
    struct Partition *P;
    
    P = malloc(sizeof(*(P)));
    if (P == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
    P->cls = malloc(n*sizeof(int));
    if (P->cls == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
    P->inv = malloc(n*sizeof(int));
    if (P->inv == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
    P->code = -1;
    P->cells = 0;
    return P;
}

void NewPartSpine(int Lev, int n) {
    
    if (Lev > 3) {
        Spine[Lev].part = malloc(sizeof(*(Spine[Lev].part)));
        if (Spine[Lev].part == NULL) {
            fprintf(ERRFILE, "\nError, memory not allocated.\n");
            exit(1);
        }
        Spine[Lev].part->cls = Spine[Lev-3].part->cls;
        Spine[Lev].part->inv = Spine[Lev-3].part->inv;
        Spine[Lev-3].part->cls = Spine[Lev-3].part->inv = NULL;
        Spine[Lev].part->code = -1;
        Spine[Lev].part->cells = 0;
    }
    else {
        Spine[Lev].part = NewPartition(n);
    }
}

void Place(int vtx, Candidate *Cand, Partition *Part) {
    int vtxto, vtxpos;
    
    vtxpos = Cand->invlab[vtx];
    vtxto = CanonIndices[Part->inv[vtxpos]]++;
    if (Cand->lab[vtxpos] != Cand->lab[vtxto]) {
        Cand->lab[vtxpos] = Cand->lab[vtxto];
        Cand->lab[vtxto] = vtx;
        Cand->invlab[Cand->lab[vtxpos]] = vtxpos;
        Cand->invlab[Cand->lab[vtxto]] = vtxto;
    }
    if (Part->cls[vtxto] > 1) {
        Part->cls[vtxto+1] = Part->cls[vtxto]-1;
        Part->cls[vtxto] = 1;
    }
}

boolean Prefix(Candidate *Cand1, Candidate *Cand2, int k) {
    int i;
    
    for (i=1; i<=k; i++) {
        if (Cand1->lab[Spine[k].tgtpos] != Cand2->lab[Spine[k].tgtpos]) {
            break;
        }
    }
    return (i>k);
}

int Preprocess(sparsegraph *sg,
               permnode **ring,
               Candidate *Cand,
               int n,
               Partition *Part,
               struct TracesVars* tv) {
    
    int i, j, j0, k, curr_cell, ind, ind0, ind1, ind2;
    int *sge;
    int HitClsInd, labi, nghb, value, SplInd, SplCntInd, sc, iend, CStackInd, SingInd, newcell, TraceInd;
    
#define SETPAIRSAUTANDTREE_PREPROC(arg, val) \
if (tv->build_autom) SETPAIRSAUT(arg, val) \
if (arg != val) \
orbjoin_sp_pair(tv->orbits, OrbList, n, arg, val, &tv->stats->numorbits); \
MakeTree(arg, val, sg, n, tv, FALSE);
    
    CStackInd = 0;
    for (i = 0; i < n; i += Part->cls[i]) {
        if (TheGraph[Cand->lab[i]].d == 1) {
            CStack[CStackInd++] = i;
        }
    }
    
    TraceInd = Part->cells;
    
    if (CStackInd > 0) {
        ind = 0;
        while (ind < CStackInd) {
            
            if (tv->mark > (NAUTY_INFINITY-2)) {
                memset(Markers, 0, n*sizeof(int));
                memset(MarkHitVtx, 0, n*sizeof(int));
                tv->mark = 0;
            }
            tv->mark++;
            
            curr_cell = CStack[ind++];
            ind2 = curr_cell+Part->cls[curr_cell];
            HitClsInd = 0;
            for (i = curr_cell; i < ind2; i++) {
                labi = Cand->lab[i];
                nghb = *(TheGraph[labi].e);
                
                if (TheGraph[nghb].d != 1) {
                    TheGraph[labi].d = -1;
                    TheGraph[labi].one = TRUE;
                }
                
                if (MarkHitVtx[nghb] == tv->mark) {
                    NghCounts[nghb]++;
                }
                else {
                    value = Part->inv[Cand->invlab[nghb]];
                    MarkHitVtx[nghb] = tv->mark;
                    NghCounts[nghb] = 1;
                    if (Markers[value] != tv->mark) {
                        HitCls[HitClsInd++] = value;
                        Markers[value] = tv->mark;
                        HitVtx[value] = nghb;
                        ElmHitCll[value] = 1;
                    }
                    else {
                        HitVtx[value+ElmHitCll[value]++] = nghb;
                    }
                }
            }
            tv->mark++;
            
            sort_Split_Array(HitCls,HitClsInd);
            SplInd = 0;
            SplCls[0] = n;
            for (j = 0; j < HitClsInd; j++) {
                ind1 = HitCls[j];
                if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < Part->cls[ind1])) {
                    SplCls[SplInd++] = ind1;
                }
                else {
                    ind2 = ind1+Part->cls[ind1];
                    value = NghCounts[Cand->lab[ind1++]];
                    for (i = ind1; i < ind2; i++) {
                        if (NghCounts[Cand->lab[i]] != value) {
                            SplCls[SplInd++] = HitCls[j];
                            break;
                        }
                    }
                    if (i == ind2) {
                        ind1 = HitCls[j];
                        if (TheGraph[Cand->lab[ind1]].d != 1) {
                            for (i = ind1; i < ind2; i++) {
                                value = Cand->lab[i];
                                Edge_Delete(value, NghCounts[value], Cand, tv);
                                sge = TheGraph[value].e+TheGraph[value].d;
                                if (NghCounts[value]>1) {
                                    factorial(&(tv->stats->grpsize1), &(tv->stats->grpsize2), NghCounts[value]);
                                    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                    for (j0=0; j0<NghCounts[value]-1; j0++) {
                                        SETPAIRSAUTANDTREE_PREPROC(sge[j0], sge[j0+1])
                                    }
                                    SETPAIRSAUTANDTREE_PREPROC(sge[j0], sge[0])
                                    SPECIALGENERATORS
                                    if (NghCounts[value] > 2) {
                                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                        SETPAIRSAUTANDTREE_PREPROC(sge[0], sge[1])
                                        if (tv->build_autom) {
                                            SETPAIRSAUT(sge[1], sge[0])
                                        }
                                        MakeTree(sge[1], sge[0], sg, n, tv, FALSE);
                                        SPECIALGENERATORS
                                    }
                                }
                            }
                            if (TheGraph[Cand->lab[ind1]].d == 1) {
                                CStack[CStackInd++] = ind1;
                            }
                        }
                    }
                }
            }
            
            if (SplInd) {
                
                for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
                    ind0 = SplCls[sc];
                    ind1 = ind0 + Part->cls[ind0];
                    SplCntInd = 0;
                    if (ElmHitCll[ind0] < Part->cls[ind0]) {
                        SplCnt[SplCntInd++] = 0;
                        SplPos[0] = Part->cls[ind0] - ElmHitCll[ind0];
                    }
                    
                    /* According to the numbers of neighbors of C into the current cell */
                    /* compute how many vertices in C will be placed into the same new cell */
                    iend = ind0 + ElmHitCll[ind0];
                    for (i = ind0; i < iend; i++) {
                        value = NghCounts[HitVtx[i]];
                        if (Markers[value] != tv->mark) {
                            Markers[value] = tv->mark;
                            SplCnt[SplCntInd++] = value;
                            SplPos[value] = 1;
                        }
                        else {
                            SplPos[value]++;
                        }
                    }
                    tv->mark++;
                    
                    /* Sort the values deriving from the previous step */
                    sort_Split_Array(SplCnt, SplCntInd);
                    
                    Part->cells += SplCntInd-1;
                    
                    /* Split the cell C and update the information for sizes of new cells */
                    /* Put the new cells into the stack */
                    i = ind0;
                    for (k = 0; k < SplCntInd; k++) {
                        value = SplPos[SplCnt[k]];
                        Part->cls[i] = value;
                        SplPos[SplCnt[k]] = i;
                        i += value;
                        if (i < ind1) {
                            TheTrace[TraceInd++] = i;
                        }
                    }
                    
                    /* Permute elements of the cell C */
                    iend = ind0 + ElmHitCll[ind0];
                    
                    for (i = ind0; i < iend; i++) {
                        value = HitVtx[i];
                        Edge_Delete(value, NghCounts[value], Cand, tv);
                        sge = TheGraph[value].e+TheGraph[value].d;
                        if (NghCounts[value] > 1) {
                            factorial(&(tv->stats->grpsize1), &(tv->stats->grpsize2), NghCounts[value]);
                            
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (j0=0; j0<NghCounts[value]-1; j0++) {
                                SETPAIRSAUTANDTREE_PREPROC(sge[j0], sge[j0+1])
                            }
                            SETPAIRSAUTANDTREE_PREPROC(sge[j0], sge[0])
                            SPECIALGENERATORS
                            if (NghCounts[value] > 2) {
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                SETPAIRSAUTANDTREE_PREPROC(sge[0], sge[1])
                                if (tv->build_autom) {
                                    SETPAIRSAUT(sge[1], sge[0])
                                }
                                MakeTree(sge[1], sge[0], sg, n, tv, FALSE);
                                SPECIALGENERATORS
                            }
                        }
                        
                        j = SplPos[NghCounts[value]]++;         /* where HitVtx[i] goes */
                        k = Cand->invlab[value];				/* where HitVtx[i] is in lab */
                        Cand->lab[k] = Cand->lab[j];
                        Cand->lab[j] = value;
                        Cand->invlab[value] = j;
                        Cand->invlab[Cand->lab[k]] = k;
                        NghCounts[value] = 0;
                    }
                    
                    /* Reconstruct the cell C and update the inverse partition */
                    newcell = ind1 - ElmHitCll[ind0];
                    i = newcell;
                    ind2 = newcell+Part->cls[newcell]-1;
                    do {
                        Part->inv[i] = newcell;
                        if (i == ind2) {
                            newcell = i+1;
                            if (newcell < n) ind2 = newcell+Part->cls[newcell]-1;
                        }
                    }
                    while (++i < ind1);
                    
                    for (i = ind0, k = 0; k < SplCntInd; i+=Part->cls[i], k++) {
                        if ((k > 0) || (SplCnt[0] > 0)) {
                            if (TheGraph[Cand->lab[i]].d == 1) {
                                CStack[CStackInd++] = i;
                            }
                        }
                        if (Part->cls[i] == 1) {
                            Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                            SingInd++;
                        }
                    }
                    
                }
            }
        }
        return 1;
    }
    else {
        return 0;
    }
}

int Preprocess_refine(sparsegraph *sg,
                      permnode **ring,
                      Candidate *Cand,
                      int n,
                      Partition *Part,
                      struct TracesVars* tv) {
    
    int i, j, j0, k, curr_cell, ind, ind0, ind1, ind2;
    int *sge;
    int HitClsInd, labi, nghb, value, SplInd, SplCntInd, sc, iend, CStackInd, SingInd, newcell, TraceInd;
    
#define SETPAIRSAUTANDTREE_PREPROC_REFINE(arg, val) \
MakeTree(arg, val, sg, n, tv, FALSE);
    
    CStackInd = 0;
    for (i = 0; i < n; i += Part->cls[i]) {
        if (TheGraph[Cand->lab[i]].d == 1) {
            CStack[CStackInd++] = i;
        }
    }
    
    TraceInd = Part->cells;
    
    if (CStackInd > 0) {
        ind = 0;
        while (ind < CStackInd) {
            
            if (tv->mark > (NAUTY_INFINITY-2)) {
                memset(Markers, 0, n*sizeof(int));
                memset(MarkHitVtx, 0, n*sizeof(int));
                tv->mark = 0;
            }
            tv->mark++;
            
            curr_cell = CStack[ind++];
            ind2 = curr_cell+Part->cls[curr_cell];
            HitClsInd = 0;
            for (i = curr_cell; i < ind2; i++) {
                labi = Cand->lab[i];
                nghb = *(TheGraph[labi].e);
                
                if (TheGraph[nghb].d != 1) {
                    TheGraph[labi].d = -1;
                    TheGraph[labi].one = TRUE;
                }
                
                if (MarkHitVtx[nghb] == tv->mark) {
                    NghCounts[nghb]++;
                }
                else {
                    value = Part->inv[Cand->invlab[nghb]];
                    MarkHitVtx[nghb] = tv->mark;
                    NghCounts[nghb] = 1;
                    if (Markers[value] != tv->mark) {
                        HitCls[HitClsInd++] = value;
                        Markers[value] = tv->mark;
                        HitVtx[value] = nghb;
                        ElmHitCll[value] = 1;
                    }
                    else {
                        HitVtx[value+ElmHitCll[value]++] = nghb;
                    }
                }
            }
            
            tv->mark++;
            
            sort_Split_Array(HitCls,HitClsInd);
            
            SplInd = 0;
            SplCls[0] = n;
            for (j = 0; j < HitClsInd; j++) {
                ind1 = HitCls[j];
                if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < Part->cls[ind1])) {
                    SplCls[SplInd++] = ind1;
                }
                else {
                    ind2 = ind1+Part->cls[ind1];
                    value = NghCounts[Cand->lab[ind1++]];
                    for (i = ind1; i < ind2; i++) {
                        if (NghCounts[Cand->lab[i]] != value) {
                            SplCls[SplInd++] = HitCls[j];
                            break;
                        }
                    }
                    if (i == ind2) {
                        ind1 = HitCls[j];
                        if (TheGraph[Cand->lab[ind1]].d != 1) {
                            for (i = ind1; i < ind2; i++) {
                                value = Cand->lab[i];
                                Edge_Delete(value, NghCounts[value], Cand, tv);
                                sge = TheGraph[value].e+TheGraph[value].d;
                                if (NghCounts[value]>1) {
                                    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                    for (j0=0; j0<NghCounts[value]-1; j0++) {
                                        SETPAIRSAUTANDTREE_PREPROC_REFINE(sge[j0], sge[j0+1])
                                    }
                                    SETPAIRSAUTANDTREE_PREPROC_REFINE(sge[j0], sge[0])
                                    if (NghCounts[value] > 2) {
                                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                        SETPAIRSAUTANDTREE_PREPROC_REFINE(sge[0], sge[1])
                                        MakeTree(sge[1], sge[0], sg, n, tv, FALSE);
                                    }
                                }
                            }
                            if (TheGraph[Cand->lab[ind1]].d == 1) {
                                CStack[CStackInd++] = ind1;
                            }
                        }
                    }
                }
            }
            
            if (SplInd) {
                
                for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
                    ind0 = SplCls[sc];
                    ind1 = ind0 + Part->cls[ind0];
                    SplCntInd = 0;
                    if (ElmHitCll[ind0] < Part->cls[ind0]) {
                        SplCnt[SplCntInd++] = 0;
                        SplPos[0] = Part->cls[ind0] - ElmHitCll[ind0];
                    }
                    
                    /* According to the numbers of neighbors of C into the current cell */
                    /* compute how many vertices in C will be placed into the same new cell */
                    iend = ind0 + ElmHitCll[ind0];
                    for (i = ind0; i < iend; i++) {
                        value = NghCounts[HitVtx[i]];
                        if (Markers[value] != tv->mark) {
                            Markers[value] = tv->mark;
                            SplCnt[SplCntInd++] = value;
                            SplPos[value] = 1;
                        }
                        else {
                            SplPos[value]++;
                        }
                    }
                    tv->mark++;
                    
                    /* Sort the values deriving from the previous step */
                    sort_Split_Array(SplCnt, SplCntInd);
                    
                    Part->cells += SplCntInd-1;
                    
                    /* Split the cell C and update the information for sizes of new cells */
                    /* Put the new cells into the stack */
                    i = ind0;
                    for (k = 0; k < SplCntInd; k++) {
                        value = SplPos[SplCnt[k]];
                        Part->cls[i] = value;
                        SplPos[SplCnt[k]] = i;
                        i += value;
                        if (i < ind1) {
                            TheTrace[TraceInd++] = i;
                        }
                    }
                    
                    /* Permute elements of the cell C */
                    iend = ind0 + ElmHitCll[ind0];
                    
                    for (i = ind0; i < iend; i++) {
                        value = HitVtx[i];
                        Edge_Delete(value, NghCounts[value], Cand, tv);
                        sge = TheGraph[value].e+TheGraph[value].d;
                        if (NghCounts[value] > 1) {
                            
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (j0=0; j0<NghCounts[value]-1; j0++) {
                                SETPAIRSAUTANDTREE_PREPROC_REFINE(sge[j0], sge[j0+1])
                            }
                            SETPAIRSAUTANDTREE_PREPROC_REFINE(sge[j0], sge[0])
                            if (NghCounts[value] > 2) {
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                SETPAIRSAUTANDTREE_PREPROC_REFINE(sge[0], sge[1])
                                MakeTree(sge[1], sge[0], sg, n, tv, FALSE);
                            }
                        }
                        
                        j = SplPos[NghCounts[value]]++;         /* where HitVtx[i] goes */
                        k = Cand->invlab[value];				/* where HitVtx[i] is in lab */
                        Cand->lab[k] = Cand->lab[j];
                        Cand->lab[j] = value;
                        Cand->invlab[value] = j;
                        Cand->invlab[Cand->lab[k]] = k;
                        NghCounts[value] = 0;
                    }
                    
                    /* Reconstruct the cell C and update the inverse partition */
                    newcell = ind1 - ElmHitCll[ind0];
                    i = newcell;
                    ind2 = newcell+Part->cls[newcell]-1;
                    do {
                        Part->inv[i] = newcell;
                        if (i == ind2) {
                            newcell = i+1;
                            if (newcell < n) ind2 = newcell+Part->cls[newcell]-1;
                        }
                    }
                    while (++i < ind1);
                    
                    for (i = ind0, k = 0; k < SplCntInd; i+=Part->cls[i], k++) {
                        if ((k > 0) || (SplCnt[0] > 0)) {
                            if (TheGraph[Cand->lab[i]].d == 1) {
                                CStack[CStackInd++] = i;
                            }
                        }
                        if (Part->cls[i] == 1) {
                            Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                            SingInd++;
                        }
                    }
                }
            }
        }
        return 1;
    }
    else {
        return 0;
    }
}

void PrintPartition(int *v, int *cls, int n, int l, int line) {
    int i, j, k;
    
    k=0;
    fprintf(outfile, "[ ");
    for (i=0; i<n; i+=cls[i]) {
        if ((cls[i]<=0) || i>=n) {
            printf("WRONG");
            break;
        }
        for (j=i; j<i+cls[i]; j++) {
            fprintf(outfile, "%d ", v[j]+l);
            if (k++ > 50) {
                fprintf(outfile,"\n");
                k=0;
            }
        }
        if ((i+cls[i])<n) fprintf(outfile, "| ");
    }
    fprintf(outfile, "] at line %d\n", line);
    return;
}

void PrintVect(int *v, int z, int n, int l) {
    int i;
    printf("[");
    for (i = z; i<n; i++)
        printf(" %2d", v[i]+l);
    printf(" ]\n");
    return;
}

void PrintWeightedGraph1(sparsegraph *g_arg, int n, char msg[30]) {
    int i, j;
    int *ngh1, *wgh1;
    
    printf("%s\n",msg);
    for (i=0; i<n; i++) {
        ngh1 = g_arg->e+g_arg->v[i];
        wgh1 = g_arg->w+g_arg->v[i];
        printf("%2d: ",i+labelorg);
        for (j=0; j<g_arg->d[i]; j++) {
            printf("%2d ",ngh1[j]);
            printf("(%d) ",wgh1[j]);
        }
        printf("\n");
    }
    printf("\n");
}

void PrintWeightedGraph2(int n, char msg[30]) {
    int i, j;
    int *ngh1, *wgh1;
    
    printf("%s\n",msg);
    for (i=0; i<n; i++) {
        ngh1 = TheGraph[i].e;
        printf("%2d: ",i+labelorg);
        for (j=0; j<TheGraph[i].d; j++) {
            printf("%2d ",ngh1[j]+labelorg);
        }
        printf(";\n");
    }
    printf("\n");
}

void PrintBlissGraph(int n) {
    int i, j;
    int *ngh1, *wgh1;
    
    fprintf(outfile,"p edge %d\n",n);
    for (i=0; i<n; i++) {
        ngh1 = TheGraph[i].e;
        for (j=0; j<TheGraph[i].d; j++) {
            if (i < ngh1[j]) {
                fprintf(outfile, "e %d %d\n",i+labelorg,ngh1[j]+labelorg);
            }
        }
    }
    printf("\n");
}

void putgraphplus_sg(FILE *f, sparsegraph *sg, int linelength)
{
    int i,n,curlen,slen;
    int *d,*e;
    size_t *v,j;
    char s[60];
    
    n = sg->nv;
    SG_VDE(sg,v,d,e);
    
    for (i = 0; i < n; ++i)
    {
        fprintf(f,"%3d : ",i+labelorg);
        curlen = 7;
        
        for (j = v[i]; j < v[i]+d[i]; ++j)
        {
            if (sg->w) {
                if (sg->w[j] != 1) {
                    slen = itos(sg->w[j],s);
                    if (linelength > 0 && curlen + slen + 1 > linelength)
                    {
                        putstring(f,"\n  ");
                        curlen = 2;
                    }
                    PUTC(' ',f);
                    PUTC('w',f);
                    putstring(f,s);
                    curlen += slen + 3;
                }
            }
            
            slen = itos(e[j]+labelorg,s);
            if (linelength > 0 && curlen + slen + 1 > linelength)
            {
                putstring(f,"\n  ");
                curlen = 2;
            }
            PUTC(' ',f);
            putstring(f,s);
            curlen += slen + 1;
        }
        putstring(f,";\n");
    }
}

void quickSort(int *arr, int elements) {
    
#define MAX_LEVELS 300
    
    int piv, beg[MAX_LEVELS], end[MAX_LEVELS], i = 0, L, R, swap;
    int k, value;
    
    beg[0] = 0;
    end[0] = elements;
    while (i>= 0) {
        L = beg[i];
        R = end[i]-1;
        if (L<R-8) {
            piv = arr[(L+R)/2];
            arr[(L+R)/2] = arr[L];
            arr[L] = piv;
            while (L<R) {
                while (arr[R]>= piv && L<R) R--;
                if (L<R) arr[L++] = arr[R];
                while (arr[L]<= piv && L<R) L++;
                if (L<R) arr[R--] = arr[L];
            }
            arr[L] = piv;
            beg[i+1] = L+1;
            end[i+1] = end[i]; end[i++] = L;
            if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
                swap = beg[i];
                beg[i] = beg[i-1];
                beg[i-1] = swap;
                swap = end[i];
                end[i] = end[i-1];
                end[i-1] = swap;
            }
        }
        else {
            i--;
        }
    }
    for (k = 1; k < elements; ++k) {
        value = arr[k];
        i = k - 1;
        while ((i >= 0) && (value < arr[i])) {
            arr[i + 1] = arr[i];
            --i;
        }
        arr[i + 1] = value;
    }
}

void RemoveFromLevel(int from, int to, int strategy, boolean reinit) {
    int i;
    
    for (i=from; i<=to; i++) {
        if (Spine[i].listend) {
            (Spine[i].listend)->next = GarbList;
            GarbList = Spine[i].liststart;
            Spine[i].liststart = Spine[i].listend = NULL;
        }
        if (strategy == 0 || reinit) {
            Spine[i].listcounter = 0;
            if (i>from) {
                Spine[i].thetracexists = FALSE;
                Spine[i].part->code = -1;
            }
        }
    }
}

searchtrie *searchtrie_make(Candidate *CurrCand, Candidate *NextCand, int n, struct TracesVars *tv) {
    
    searchtrie *st;
    if (tv->strienext == n) {
        tv->strienext = 0;
        tv->strielist->next = malloc(sizeof(struct trielist));
        if (tv->strielist->next == NULL) {
            fprintf(ERRFILE, "\nError, memory not allocated.\n");
            exit(1);
        }
        tv->strielist->next->prev = tv->strielist;
        tv->strielist = tv->strielist->next;
        tv->strielist->next = NULL;
        tv->strielist->triearray = malloc(n*sizeof(searchtrie));
        if (tv->strielist->triearray == NULL) {
            fprintf(ERRFILE, "\nError, memory not allocated.\n");
            exit(1);
        }
    }
    st = &(tv->strielist->triearray[tv->strienext]);
    st->father = CurrCand->stnode;
    st->name = NextCand->name;
    st->index = tv->newindex+1;
    st->vtx = NextCand->vertex;
    st->level = tv->tolevel;
    st->first_child = st->next_sibling = st->last_child = st->goes_to = NULL;
    if (st->father) {
        if (st->father->first_child) {
            st->father->last_child->next_sibling = st;
            st->father->last_child = st;
        }
        else {
            st->father->first_child = st->father->last_child = st;
        }
    }
    NextCand->stnode = st;
    if (tv->newgotonode) {
        tv->newgotonode->goes_to = st;
    }
    if (tv->gotonode) {
        st->goes_to = tv->gotonode;
        tv->gotonode = NULL;
    }
    tv->strienext++;
    return st;
}

trielist *searchtrie_new(int n, struct TracesVars *tv) {
    
    tv->strielist = malloc(sizeof(struct trielist));
    if (tv->strielist == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
    tv->strielist->prev = tv->strielist->next = NULL;
    tv->strielist->triearray = malloc(n*sizeof(searchtrie));
    if (tv->strielist->triearray == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
    tv->strielist->triearray[0].father = tv->strielist->triearray[0].first_child = NULL;
    tv->strielist->triearray[0].next_sibling = tv->strielist->triearray[0].last_child = NULL;
    tv->strielist->triearray[0].goes_to = NULL;
    tv->strielist->triearray[0].index = 1;
    tv->strielist->triearray[0].name = tv->strielist->triearray[0].level = 0;
    tv->strielist->triearray[0].vtx = n;
    
    tv->strienext = 1;
    return tv->strielist;
}

int Select_from_CStack(int *cls, int CStackInd) {
    int j, k;
    
    j = CStackInd;
    k = CStackInd;
    while (--j > 0) {
        if (cls[CStack[j]] < cls[CStack[k]]) {
            k = j;
        }
        if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
            break;
        }
    }
    return k;
}

boolean SelectNextLevel(int n, struct TracesVars *tv, struct TracesInfo *ti) {
    int i, j, val;
    Candidate *FirstCand;
    boolean orbitcell;
    
    switch (tv->compstage) {
        case 2:
            tv->nextlevel = tv->maxtreelevel;
            while (tv->nextlevel >=0) {
                if (Spine[tv->nextlevel].liststart) {
                    break;
                }
                tv->nextlevel--;
            }
            if (tv->nextlevel < 0) {
                return FALSE;
            }
            break;
        default:
            switch (tv->strategy) {
                case 0:
                    tv->nextlevel = tv->fromlevel;
                    while (!Spine[tv->nextlevel].liststart) {
                        (tv->nextlevel)++;
                    }
                    PRINTF2("SelectNextLevel 1?: finalnumcells: %d; ", tv->finalnumcells);
                    PRINTF2("Spine[tv->nextlevel].part->cells: %d; ", Spine[tv->nextlevel].part->cells);
                    PRINTF2("tv->maxtreelevel: %d; ", tv->maxtreelevel);
                    PRINTF2("tv->nextlevel: %d\n", tv->nextlevel);
                    if ((Spine[tv->nextlevel].part->cells == tv->finalnumcells) || (tv->nextlevel > tv->maxtreelevel)) {
                        return FALSE;
                    } else {
                        /* Check the whole group */
                        if ((tv->group_level < tv->tolevel) && !ti->first_matching && ti->thegrouphaschanged) {
                            
                            FirstCand = Spine[tv->nextlevel].liststart;
                            val = tv->orbits[FirstCand->lab[Spine[1].tgtcell]];
                            for (i=Spine[1].tgtcell; i<Spine[1].tgtend; i++) {
                                if (tv->orbits[FirstCand->lab[i]] != val) {
                                    break;
                                }
                            }
                            if (i<Spine[1].tgtend) {
                                orbitcell = FALSE;
                            } else {
                                orbitcell = TRUE;
                            }
                            if (orbitcell) {
                                FirstCand = Spine[tv->nextlevel].liststart;
                                FixBase(fix, tv, FirstCand, 0, tv->firstpathlength);
                                if (tv->options->verbosity >= 2) tv->schreier1 -= CPUTIME;
                                getorbitsmin(fix, tv->nfix, gpB, &gensB, &tv->currorbit, NULL, n, n, TRUE);
                                if (tv->options->verbosity >= 2) tv->schreier1 += CPUTIME;
                                for (j=2; j<=tv->firstpathlength; j++) {
                                    tv->currorbit = findcurrorbits(gpB, j-1);
                                    val = tv->currorbit[FirstCand->lab[Spine[j].tgtcell]];
                                    for (i=Spine[j].tgtcell; i<Spine[j].tgtend; i++) {
                                        if (tv->currorbit[FirstCand->lab[i]] != val) {
                                            break;
                                        }
                                    }
                                    if (i<Spine[j].tgtend) {
                                        break;
                                    }
                                }
                                tv->group_level = j-1;
                                if (tv->group_level >= tv->tolevel) {
                                    ti->thegrouphaschanged = FALSE;
                                }
                            }
                            /* End check the whole group */
                        }
                        
                    }
                    break;
                case 1:
                    tv->nextlevel = tv->maxtreelevel;
                    PRINTF2("SelectNextLevel 2?: finalnumcells: %d; ", tv->finalnumcells);
                    PRINTF2("Spine[tv->nextlevel].part->cells: %d; ", Spine[tv->nextlevel].part->cells);
                    if (Spine[tv->nextlevel].part->cells == tv->finalnumcells) {
                        (tv->nextlevel)--;
                    }
                    while (tv->nextlevel >= 0) {
                        if (Spine[tv->nextlevel].liststart) {
                            break;
                        }
                        tv->nextlevel--;
                    }
                    if (tv->nextlevel < 0) {
                        return FALSE;
                    }
                    break;
                default:
                    break;
            }
            break;
    }
    return TRUE;
}

void SetAutom(int q, int n, struct TracesVars *tv) {
    int i;
    
    for (i=0; i<q; i++) {
        AUTPERM[PrmPairs[i].arg] = PrmPairs[i].val;
    }
    return;
}

void ResetAutom(int q, int n, struct TracesVars *tv) {
    int i;
    
    if (n/q < 256) {
        memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
    }
    else {
        for (i=0; i<q; i++) {
            AUTPERM[PrmPairs[i].arg] = PrmPairs[i].arg;
        }
    }
    tv->permInd = 0;
    return;
}

void sort_Split_Array(int *Array, int Ind){
    int i, k, value;
    
    switch (Ind) {
        case 0:
        case 1:
            break;
        case 2:
            if (Array[0] > Array[1]) {
                value = Array[0];
                Array[0] = Array[1];
                Array[1] = value;
            }
            break;
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
            for (k = 1; k < Ind; ++k) {
                value = Array[k];
                i = k - 1;
                while ((i >= 0) && (value < Array[i])) {
                    Array[i + 1] = Array[i];
                    --i;
                }
                Array[i + 1] = value;
            }
            break;
        default:
            quickSort(Array, Ind);
            break;
    }
}

int spinelementorbsize(int *orbits, int *lab, int size, int elem) {
    int i, j, val;
    
    j = 0;
    val = orbits[elem];
    for (i = 0; i < size; ++i) {
        if (orbits[lab[i]] == val) ++j;
    }
    return j;
}

boolean TargetCell(Candidate *TargCand, Partition *Part, int n, struct TracesVars* tv, int Lv) {
    int TCell = -1, TCSize = 1;
    int i;
    
    if (Part->cells == n) {
        tv->finalnumcells = n;
        return FALSE;
    }
    if (tv->maxdeg <=2) {
        return FALSE;
    }
    
    if (Lv < tv->tcellevel) {
        tv->tcell = Spine[Lv+1].tgtcell;
        return TRUE;
    }
    else {
        if (Part->cls[0] == n) {
            tv->tcell = 0;
            return TRUE;
        }
        while (TCell < 0) {
            for (i = Spine[Lv].tgtcell; i < Spine[Lv].tgtend; i += Part->cls[i]) {
                if (Part->cls[i] > TCSize) {
                    if (NonSingDeg(TargCand->lab[i], TargCand, Part) > 2) {
                        TCSize = Part->cls[i];
                        TCell = i;
                    }
                }
            }
            Lv--;
            if ((Lv < 0) && (TCell < 0)) return FALSE;
        }
        tv->tcell = TCell;
        return TRUE;
    }
}

int TargetCellExpPath(Candidate *TargCand, Partition *Part, struct TracesVars* tv) {
    int Lv, n;
    
    n = tv->input_graph->nv;
    if (Part->cells == n) {
        return 0;
    }
    
    Lv = tv->tolevel_tl+1;
    SpineTL_tl = Spine+Lv;
    if (tv->tolevel_tl < tv->tcellevel) {
        tv->tcellexpath = Part->inv[SpineTL_tl->tgtcell];
        tv->tolevel_tl++;
        if (Part->cls[tv->tcellexpath] == 1) {
            if (tv->options->verbosity >= 2) {
                if (tv->tolevel_tl-tv->tolevel == 6) {
                    fprintf(outfile, "... ");
                }
            }
            return TargetCellExpPath(TargCand, Part, tv);
        } else {
            return 1+((Spine[tv->tolevel_tl].tgtcell >= Spine[tv->tolevel_tl-1].tgtcell) && (Spine[tv->tolevel_tl].tgtend <= Spine[tv->tolevel_tl-1].tgtend));
        }
    }
    else {
        if (TargetCellFirstPath(TargCand, Part, tv)) {
            return 1+((Spine[tv->tolevel_tl].tgtcell >= Spine[tv->tolevel_tl-1].tgtcell) && (Spine[tv->tolevel_tl].tgtend <= Spine[tv->tolevel_tl-1].tgtend));
        }
        else {
            return 0;
        }
    }
}

boolean TargetCellFirstPath(Candidate *TargCand, Partition *Part, struct TracesVars* tv) {
    int n, TCell, TCSize, TCell1, TCSize1;
    int Lv, i, Lev, vtx, vtx_d;
    int loopstart, loopend;
    boolean divided;
    
    n = tv->input_graph->nv;
    
    if (Part->cells == n) {
        return 0;
    }
    
    Lev = tv->tolevel_tl;
    Lv = tv->tolevel_tl;
    
    TCell = TCell1 = -1;
    TCSize = TCSize1 = 1;
    
    while (TCell < 0) {
        
        if (tv->compstage == 2) {
            loopstart = Spine[Lv].tgtcell;
            divided = FALSE;
        } else {
            loopstart = Part->inv[Spine[Lv].tgtcell];
            divided = FALSE;
            
            if (Lv == tv->lastlev) {
                loopstart = Part->inv[tv->lastcell];
                divided = TRUE;
            }
        }
        
        i = loopstart;
        loopend = Spine[Lv].tgtend;
        while (i < loopend) {
            if (Part->cls[i] > TCSize) {
                vtx = TargCand->lab[i];
                vtx_d = TheGraph[vtx].d;
                if (vtx_d > 2) {
                    if (NonSingDeg(vtx, TargCand, Part) > 2) {
                        TCSize = Part->cls[i];
                        TCell = i;
                        if (TCSize == WorkArray[Lv]) {
                            break;
                        }
                    }
                }
            }
            i += Part->cls[i];
            if (divided && (i == loopend)) {
                i = loopstart = Part->inv[Spine[Lv].tgtcell];
                loopend = tv->lastcell;
                divided = FALSE;
                TCSize1 = TCSize;
                TCell1 = TCell;
                TCell = -1;
                TCSize = 1;
            }
        }
        
        if (TCSize1 > TCSize) {
            TCell = TCell1;
            TCSize = TCSize1;
        }
        
        if (TCell < 0) {
            if (Lv == 0) {
                tv->finalnumcells = min(Part->cells,tv->finalnumcells);    /* 160712 */
                return FALSE;
            } else {
                Lv = Spine[Lv].tgtfrom;
            }
        }
    }
    
    tv->tcellexpath = tv->lastcell = TCell;
    tv->tolevel_tl++;
    
    Spine[tv->tolevel_tl].tgtfrom = tv->lastlev = Lv;
    Spine[tv->tolevel_tl].tgtcell = tv->tcellexpath;
    Spine[tv->tolevel_tl].tgtsize = WorkArray[Lv] = TCSize;
    Spine[tv->tolevel_tl].tgtend = Spine[tv->tolevel_tl].tgtcell + TCSize;
    Spine[tv->tolevel_tl].tgtpos = Spine[tv->tolevel_tl].tgtend - 1;
    tv->tcellevel = tv->tolevel_tl;
    
    if (Lv != Lev) {
        BreakSteps[Lev] = ++tv->brkstpcount;
        if (Spine[tv->tolevel].liststart) {
            if (!Spine[tv->tolevel].liststart->firstsingcode) {
                Spine[tv->tolevel].liststart->firstsingcode = Spine[tv->tolevel].liststart->pathsingcode;
            }
        }
    }
    return TRUE;
}

void traces_freedyn(void) {
    /* Free the static dynamic memory used by Traces */
#if !MAXN
    DYNFREE(AUTPERM, AUTPERM_sz);
    DYNFREE(BreakSteps, BreakSteps_sz);
    DYNFREE(CStack, CStack_sz);
    DYNFREE(CurrOrbSize, CurrOrbSize_sz);
    DYNFREE(CurrRefCells, CurrRefCells_sz);
    DYNFREE(Diff, Diff_sz);
    DYNFREE(Factorials, Factorials_sz);
    DYNFREE(fix, fix_sz);
    DYNFREE(IDENTITY_PERM, IDENTITY_PERM_sz);
    DYNFREE(Markers, Markers_sz);
    DYNFREE(TreeMarkers, TreeMarkers_sz);
    DYNFREE(AutMarkers, AutMarkers_sz);
    DYNFREE(MarkHitVtx, MarkHitVtx_sz);
    DYNFREE(MultRefCells, MultRefCells_sz);
    DYNFREE(NghCounts, NghCounts_sz);
    DYNFREE(OrbSize, OrbSize_sz);
    DYNFREE(OrbList, OrbList_sz);
    DYNFREE(PrmPairs, PrmPairs_sz);
    DYNFREE(TempOrbList, TempOrbList_sz);
    DYNFREE(RefCells, RefCells_sz);
    DYNFREE(RefPath, RefPath_sz);
    DYNFREE(Singletons, Singletons_sz);
    DYNFREE(SplCls, SplCls_sz);
    DYNFREE(SplCnt, SplCnt_sz);
    DYNFREE(SplPos, SplPos_sz);
    DYNFREE(StackMarkers, StackMarkers_sz);
    DYNFREE(TheTrace, TheTrace_sz);
    DYNFREE(TheTraceCC, TheTraceCC_sz);
    DYNFREE(TheTraceSplNum, TheTraceSplNum_sz);
    DYNFREE(TheTraceSteps, TheTraceSteps_sz);
    DYNFREE(TEMPLAB, TEMPLAB_sz);
    DYNFREE(TEMPINVLAB, TEMPINVLAB_sz);
    DYNFREE(WeightsSeq, WeightsSeq_sz);
    DYNFREE(WorkArray, WorkArray_sz);
    DYNFREE(WorkArray0, WorkArray0_sz);
    DYNFREE(WorkArray1, WorkArray1_sz);
    DYNFREE(WorkArray2, WorkArray2_sz);
    DYNFREE(WorkArray3, WorkArray3_sz);
    DYNFREE(WorkArray4, WorkArray4_sz);
    DYNFREE(WorkArray5, WorkArray5_sz);
    DYNFREE(WorkArray6, WorkArray6_sz);
    DYNFREE(WorkArray7, WorkArray7_sz);
    DYNFREE(Neighbs1, Neighbs1_sz);
    DYNFREE(Neighbs2, Neighbs2_sz);
    DYNFREE(TreeStack, TreeStack_sz);
    DYNFREE(Spine, Spine_sz);
    DYNFREE(TrieArray, TrieArray_sz);
    DYNFREE(TheGraph, TheGraph_sz);
    DYNFREE(EPCodes, EPCodes_sz);
#endif
}

boolean TreeFyTwo(int From, Candidate *Cand1, Candidate *Cand2, Partition *Part, int n,
                  struct TracesVars* tv, struct TracesInfo *ti) {
    int i, i1, i2, j1, j2, k;
    int vtx1, vtx2, ngh1, ngh2, arg, val;
    int *tgtc1, *tgtc2, *adj1, *adj2;
    int iend;
    
    SETMARK(Markers, tv->mark)
    
    i2=0;
    
    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
    i1 = Spine[From].tgtsize;
    tgtc1 = Cand1->lab+Spine[From].tgtcell;
    tgtc2 = Cand2->lab+Spine[From].tgtcell;
    for (i=0; i<i1; i++) {
        arg = tgtc1[i];
        val = tgtc2[i];
        if ((Markers[arg] != tv->mark) && (Markers[val] != tv->mark)) {
            SETPAIRSAUT(arg, val)
            SETPAIRSAUT(val, arg)
            Markers[arg] = tv->mark;
            Markers[val] = tv->mark;
        } else {
            return FALSE;    /*  160715  */
        }
    }
    
    while (i2 < tv->permInd) {
        vtx1 = PrmPairs[i2].arg;
        vtx2 = PrmPairs[i2++].val;
        adj1 = TheGraph[vtx1].e;
        adj2 = TheGraph[vtx2].e;
        iend = TheGraph[vtx1].d;
        j1 = j2 = 0;
        for (k=0; k < iend; k++) {
            ngh1 = adj1[k];
            if (Markers[ngh1] != tv->mark) {
                Neighbs1[j1++] = Cand1->invlab[ngh1];
            }
            ngh2 = adj2[k];
            if (Markers[ngh2] != tv->mark) {
                Neighbs2[j2++] = Cand2->invlab[ngh2];
            }
        }
        
        k = tv->permInd;
        if (j1 == j2) {
            quickSort(Neighbs1, j1);
            quickSort(Neighbs2, j2);
            for (i=0; i<j1; i++) {
                arg = Cand1->lab[Neighbs1[i]];
                val = Cand2->lab[Neighbs2[i]];
                if ((Markers[arg] != tv->mark) && (Markers[val] != tv->mark)) {
                    SETPAIRSAUT(arg, val)
                    SETPAIRSAUT(val, arg)
                    Markers[arg] = tv->mark;
                    Markers[val] = tv->mark;
                }
            }
        }
    }
    return TRUE;
}

/*****************************************************************************
 *                                                                            *
 *  trie_class(t,c) classifies vertices according to weights of their edges   *
 *                                                                            *
 *****************************************************************************/

void  trie_class(trie *t, int *count) {
    
    if (t->first_child == NULL) {
        WeightsSeq[t->value] = *count;
        if (t->next_sibling == NULL) (*count)++;
        return;
    }
    else {
        t = t->first_child;
        while (t) {
            trie_class(t,count);
            t = t->next_sibling;
        }
    }
}

int trie_classify(int n, TracesVars *tv) {
    
    int i, j, ord;
    int *ngh1, *wgh1;
    
    trieroot = trie_new(n, tv);
    ord = 0;
    
    for (i=0; i<n; i++) {
        ngh1 = TheGraph[i].e;
        wgh1 = TheGraph[i].w;
        sort2ints(wgh1, ngh1, TheGraph[i].d);
        
        trieref = trieroot;
        for (j=0; j<TheGraph[i].d; j++) {
            trieref = trie_make(trieref, wgh1[j], n, tv);
        }
        trieref = trie_make(trieref, n, n, tv);
        trie_make(trieref, i, n, tv);
    }
    trie_class(trieroot,&ord);
    
    for (i=0; i<=tv->triepos; i++) {
        free(TrieArray[i]);
    }
    trieroot = NULL;
    return ord-1;
}

struct trie *trie_comp(trie *t, int value) {
    
    if (t->first_child) {
        t = t->first_child;
        while (t) {
            if  (value != t->value) {
                t = t->next_sibling;
            }
            else {
                break;
            }
        }
        return t;
    }
    else {
        return NULL;
    }
}

void  trie_dump(trie *t) {
    if (t->first_child == NULL) return;
    else {
        printf("( ");
        t = t->first_child;
        while (t) {
            printf("%d ",t->value);
            trie_dump(t);
            t = t->next_sibling;
        }
        printf(") ");
    }
}

/*****************************************************************************
 *                                                                            *
 *  trie_make(t,v,n,tv) places the value v into the trie t                    *
 *                                                                            *
 *****************************************************************************/

struct trie *trie_make(trie *t, int value, int n, struct TracesVars* tv) {
    trie *t1;
    
    t1 = t;
    if (tv->trienext == n) {
        tv->trienext = 0;
        tv->triepos++;
        TrieArray[tv->triepos] = malloc(n*sizeof(trie));
        if (TrieArray[tv->triepos] == NULL) {
            fprintf(ERRFILE, "\nError, memory not allocated.\n");
            exit(1);
        }
    }
    if (t->first_child) {
        t = t->first_child;
        if (value < t->value) {
            t1->first_child = &TrieArray[tv->triepos][tv->trienext++];
            t1->first_child->next_sibling = t;
            t1->first_child->first_child = NULL;
            t = t1->first_child;
            t->value = value;
            return t;
        }
        while (value > t->value) {
            t1 = t;
            if (t->next_sibling) {
                t = t->next_sibling;
            }
            else break;
        }
        if (value == t->value) {
            return t;
        }
        t1->next_sibling = &TrieArray[tv->triepos][tv->trienext++];
        t1->next_sibling->first_child = t1->next_sibling->next_sibling = NULL;
        if (t != t1) {
            t1->next_sibling->next_sibling = t;
        }
        t = t1->next_sibling;
    }
    else {
        t->first_child = &TrieArray[tv->triepos][tv->trienext++];
        t = t->first_child;
        t->first_child = t->next_sibling = NULL;
    }
    t->value = value;
    return t;
}

struct trie *trie_new(int n, struct TracesVars* tv) {
    
    TrieArray[0] = malloc(n*sizeof(trie));
    if (TrieArray[0] == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
    TrieArray[0][0].first_child = TrieArray[0][0].next_sibling = NULL;
    tv->triepos = 0;
    tv->trienext = 1;
    return TrieArray[0];
}

boolean VerifyCand(Candidate *Cand, int n, int line) {
    int i, k;
    
    for (i=0; i<n; i++) {
        k=Cand->lab[i];
        if (Cand->invlab[k] != i) {
            printf("Cand->invlab wrong at %d (vtx: %d), line %d\n", i, k, line);
            PrintVect(Cand->lab, 0, n, 0);
            PrintVect(Cand->invlab, 0, n, 0);
            return FALSE;
        }
    }
    return TRUE;
}

boolean VerifyId(int *p, int n) {
    int i, r;
    
    r = TRUE;
    for (i=0; i<n; i++) {
        if (p[i] != i) {
            printf("p[%d] = %d\n", i, p[i]);
            r = FALSE;
        }
    }
    return r;
}

boolean VerifyPart(Partition *Part, int start, int end) {
    int i,j;
    
    for (i=start; i<end; i+=Part->cls[i]) {
        if (Part->cls[i] == 0 || i>=end) {
            printf("WRONG cls\n");
            return FALSE;
        }
        for (j=0; j<Part->cls[i]; j++) {
            if (Part->inv[i+j] != i) {
                printf("WRONG inv\n");
                return FALSE;
            }
        }
    }
    printf("OK\n");
    return TRUE;
}

int VerifyPerm(int *perm, int n,int where) {
    int i;
    
    memset(Markers, 0, n*sizeof(int));
    
    for (i=0; i<n; i++) {
        if ((perm[i] >= n) || (Markers[perm[i]])) {
            fprintf(stderr,"wrong permutation @ %d\n",where);
            PrintVect(perm,0,i+1,labelorg);
        }
        Markers[perm[i]] = TRUE;
    }
    return TRUE;
}

/*****************************************************************************
 *                                                                            *
 *  WeightCodes(n) transforms the weight w(u,v) of an edge into a new weight  *
 *  W(u,v) such that W(a,b) = W(c,d) iff w(a,b) = w(c,d) and w(b,a) = w(d,c). *
 *                                                                            *
 *****************************************************************************/

void WeightCodes(int n) {
    int i,j,aux;
    int sumdegs;
    int deg, vtx1, vtx2, *ngh1, *ngh2, *wgh1, *wgh2, ord;
    
    sumdegs = 0;
    for (i=0; i<n; i++) {
        sumdegs += TheGraph[i].d;
    }
    
    DYNALLSTAT(int, VArray, VArray_sz);
    DYNALLSTAT(weightwhere, WArray, WArray_sz);
    DYNALLSTAT(grph_strct, TheAuxGraph, TheAuxGraph_sz);
    
    DYNALLOC1(int, VArray, VArray_sz, sumdegs, "WeightCodes");
    DYNALLOC1(weightwhere, WArray, WArray_sz, sumdegs, "WeightCodes");
    DYNALLOC1(grph_strct, TheAuxGraph, TheAuxGraph_sz, n, "WeightCodes");
    
    memcpy(TheAuxGraph,TheGraph,n*sizeof(grph_strct));
    
    ord = 0;
    for (vtx1=0; vtx1<n; vtx1++) {
        ngh1 = (TheAuxGraph[vtx1].e)++;
        wgh1 = TheAuxGraph[vtx1].w;
        deg = TheAuxGraph[vtx1].d;
        for (i=0; i<deg; i++) {
            vtx2 = ngh1[i];
            ngh2 = (TheAuxGraph[vtx2].e)++;
            wgh2 = (TheAuxGraph[vtx2].w)++;
            (TheAuxGraph[vtx1].d)--;
            (TheAuxGraph[vtx2].d)--;
            VArray[ord] = wgh1[i];
            WArray[ord].weight = wgh2[0];
            WArray[ord++].ref = (TheAuxGraph[vtx1].w)++;
            VArray[ord] = wgh2[0];
            WArray[ord].weight = wgh1[i];
            WArray[ord++].ref = wgh2;
        }
    }
    
    sortweights(VArray,WArray,ord);
    
    /* swap VArray and WArray.weight */
    for (i=0; i<sumdegs; i++) {
        aux = VArray[i];
        VArray[i] = WArray[i].weight;
        WArray[i].weight = aux;
    }
    
    i = j = 0;
    do {
        if (WArray[i].weight == WArray[j].weight) {
            j++;
        } else {
            sortweights(VArray+i,WArray+i,j-i);
            i = j;
        }
    } while (j<sumdegs);
    sortweights(VArray+i,WArray+i,j-i);
    
    /* weight class */
    ord = 0;
    *(WArray[0].ref) = 0;
    for (i=1; i<sumdegs; i++) {
        if ((WArray[i].weight != WArray[i-1].weight) || (VArray[i] != VArray[i-1])) {
            ord++;
        }
        *(WArray[i].ref) = ord;
    }
    
    DYNFREE(VArray, VArray_sz);
    DYNFREE(WArray, WArray_sz);
    DYNFREE(TheAuxGraph, TheAuxGraph_sz);
    
}



boolean TargetCellSmall(Candidate *TargCand, Partition *Part, int n, struct TracesVars* tv, int Lv) {
    int TCell = -1, TCSize = n;
    int i;
    
    if (tv->maxdeg <=2) {
        return FALSE;
    }
    
    if (Lv < tv->tcellevel) {
        tv->tcell = Spine[Lv+1].tgtcell;
        return TRUE;
    }
    else {
        if (Part->cls[0] == n) {
            tv->tcell = 0;
            return TRUE;
        }
        while (TCell < 0) {
            for (i = Spine[Lv].tgtcell; i < Spine[Lv].tgtend; i += Part->cls[i]) {
                if (Part->cls[i] < TCSize) {
                    if (NonSingDeg(TargCand->lab[i], TargCand, Part) > 2) {
                        TCSize = Part->cls[i];
                        TCell = i;
                    }
                }
            }
            Lv--;
            if ((Lv < 0) && (TCell < 0)) return FALSE;
        }
        tv->tcell = TCell;
        return TRUE;
    }
}

int TargetCellExpPathSmall(Candidate *TargCand, Partition *Part, struct TracesVars* tv) {
    int Lv, n;
    
    n = tv->input_graph->nv;
    if (Part->cells == n) {
        return 0;
    }
    
    Lv = tv->tolevel_tl+1;
    SpineTL_tl = Spine+Lv;
    if (tv->tolevel_tl < tv->tcellevel) {
        tv->tcellexpath = Part->inv[SpineTL_tl->tgtcell];
        tv->tolevel_tl++;
        if (Part->cls[tv->tcellexpath] == 1) {
            if (tv->options->verbosity >= 2) {
                if (tv->tolevel_tl-tv->tolevel == 6) {
                    fprintf(outfile, "... ");
                }
            }
            return TargetCellExpPath(TargCand, Part, tv);
        } else {
            return 1+((Spine[tv->tolevel_tl].tgtcell >= Spine[tv->tolevel_tl-1].tgtcell) && (Spine[tv->tolevel_tl].tgtend <= Spine[tv->tolevel_tl-1].tgtend));
        }
    }
    else {
        if (TargetCellFirstPath(TargCand, Part, tv)) {
            return 1+((Spine[tv->tolevel_tl].tgtcell >= Spine[tv->tolevel_tl-1].tgtcell) && (Spine[tv->tolevel_tl].tgtend <= Spine[tv->tolevel_tl-1].tgtend));
        }
        else {
            return 0;
        }
    }
}

boolean TargetCellFirstPathSmall(Candidate *TargCand, Partition *Part, struct TracesVars* tv) {
    int n, TCell, TCSize, TCell1, TCSize1;
    int Lv, i, Lev, vtx, vtx_d;
    int loopstart, loopend;
    boolean divided;
    
    n = tv->input_graph->nv;
    if (Part->cells == n) {
        return 0;
    }
    Lev = tv->tolevel_tl;
    Lv = tv->tolevel_tl;
    TCell = TCell1 = -1;
    TCSize = TCSize1 = n;
    
    while (TCell < 0) {
        loopstart = Part->inv[Spine[Lv].tgtcell];
        divided = FALSE;
        
        if (Lv == tv->lastlev) {
            loopstart = Part->inv[tv->lastcell];
            divided = TRUE;
        }
        
        i = loopstart;
        loopend = Spine[Lv].tgtend;
        
        while (i < loopend) {
            if ((Part->cls[i] > 1) && (Part->cls[i] < TCSize)) {
                vtx = TargCand->lab[i];
                vtx_d = TheGraph[vtx].d;
                
                if (vtx_d > 2) {
                    if (NonSingDeg(vtx, TargCand, Part) > 2) {
                        TCSize = Part->cls[i];
                        TCell = i;
                        if (TCSize == WorkArray[Lv]) {
                            break;
                        }
                    }
                }
            }
            i += Part->cls[i];
            if (divided && (i == loopend)) {
                i = loopstart = Spine[Lv].tgtcell;
                loopend = tv->lastcell;
                divided = FALSE;
                TCSize1 = TCSize;
                TCell1 = TCell;
                TCell = -1;
                TCSize = n;
            }
        }
        
        if (TCSize1 < TCSize) {
            TCell = TCell1;
            TCSize = TCSize1;
        }
        
        if (TCell < 0) {
            if (Lv == 0) {
                tv->finalnumcells = min(Part->cells,tv->finalnumcells);    /* 160712 */
                return FALSE;
            } else {
                Lv = Spine[Lv].tgtfrom;
            }
        }
    }
    tv->tcellexpath = tv->lastcell = TCell;
    tv->tolevel_tl++;
    
    Spine[tv->tolevel_tl].tgtfrom = tv->lastlev = Lv;
    Spine[tv->tolevel_tl].tgtcell = tv->tcellexpath;
    Spine[tv->tolevel_tl].tgtsize = WorkArray[Lv] = TCSize;
    Spine[tv->tolevel_tl].tgtend = Spine[tv->tolevel_tl].tgtcell + TCSize;
    Spine[tv->tolevel_tl].tgtpos = Spine[tv->tolevel_tl].tgtend - 1;
    tv->tcellevel = tv->tolevel_tl;
    
    if (Lv != Lev) {
        BreakSteps[Lev] = ++tv->brkstpcount;
        if (Spine[tv->tolevel].liststart) {
            if (!Spine[tv->tolevel].liststart->firstsingcode) {
                Spine[tv->tolevel].liststart->firstsingcode = Spine[tv->tolevel].liststart->pathsingcode;
            }
        }
    }
    return TRUE;
}
