/* This is dretodot.c version 1.0 (Feb 2016), which is included in     *
 * nauty distribution (version 2.6).                                   *
 * See the file COPYRIGHT for the details of the software license.     *
 *                                                                     *
 * gcc -03 -o dretodot dretodot.c nauty.a                              *
 * nauty.a is produced by: make all                                    */

#define USAGE "dretodot [-S#:#ixF#o#m#n#r#:#[r#]d#g] [infile.dre [outfile.dot [outfile.dre]]]"

#define HELPTEXT \
" Read graphs and initial coloring in dreadnaut format.\n\
 Write graphs in dot format to outfile.dot.\n\
 If outfile.dre is given, write the input graph and the partition,\n\
 as modified by the -F and -i options, to outfile.dre. outfile.dre\n\
 is allowed to be the same file as infile.dre.\n\
    -V    Set max number of vertices (default 1000).\n\
    -E    Set max number of edges (default 5000).\n\
    -v    Set verbose mode (default NO).\n\
   -S#:# Set maximum width and height of the drawing, in inches\n\
         (default 10 x 6.18).\n\
   -i    Refine the partition before drawing (default NO).\n\
   -x    Draw the orbit partition, computed by Traces. (default NO).\n\
   -F#   Individualize vertex # (and refine the partition).\n\
   -o#   Label vertices starting at # (default 0). This can be\n\
         overridden in the input.\n\
   -m#   Set the drawing model (see http://www.graphviz.org): \n\
         0 (or any value different from 1,...,5)=dot (default 0),\n\
         1=neato, 2=fdp, 3=sfdp, 4=twopi, 5=circo.\n\
   -n#   Scale the size of vertices in the drawing (#=0,1,2; default 1).\n\
   -r#:# (-r#) Set the vertices to be drawn at the topmost level\n\
         in a hierarchical (dot model) drawing (default none).\n\
         Any sequence of -r#:# (r#) options is allowed.\n\
   -d#   Draw the graph induced by vertices at topmost level\n\
         and by vertices at distance # from them; example:\n\
         ./dretodot -n2 -r1 -r12:17 -d2 MyGraph.dre Outfile.dot.\n\
   -g    Highlight the induced subgraph into the whole graph.\n"

/*************************************************************************/

#include <math.h>
#include "gtools.h"  /* which includes nauty.h and stdio.h */
#include "traces.h"

#define SORT_OF_SORT 2
#define SORT_TYPE1 int
#define SORT_TYPE2 int
#define SORT_NAME sort2ints
#include "sorttemplates.c"

#define WRITEDREFILE if (argnum == 3) { \
fprintf(drefile, "n=%d $=%d g\n", n, labelorg); \
putgraph_sg(drefile, &g, 0); \
if (numcells > 1) { \
fprintf(drefile, "f = "); \
putptn(drefile,lab,ptn,0,0,n); }}


extern int labelorg;

/**************************************************************************/
typedef struct NodeShape {
    char color[6];
    int labcol;
} NodeShape;

typedef struct list {
    int vtx;
    struct list *next;
} list;

struct list *NewListelem(int el)
{
	struct list *L;
	
    L = malloc(sizeof(list));
    if (L == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
	L->vtx = el;
	L->next = NULL;
	return L;
}

#define FONT_NAME "Arial Narrow"
#if MAXN
int lab[MAXN];
int invlab[MAXN];
int ptn[MAXN];
NodeShape NShape[MAXN];
int CurrVertices[MAXN];
int DistStack[MAXN];
int Ranks[MAXN];
int orbits[MAXN];
#else
DYNALLSTAT(int, lab, lab_sz);
DYNALLSTAT(int, invlab, invlab_sz);
DYNALLSTAT(int, ptn, ptn_sz);
DYNALLSTAT(NodeShape, NShape, NShape_sz);
DYNALLSTAT(int, CurrVertices, CurrVertices_sz);
DYNALLSTAT(int, DistStack, DistStack_sz);
DYNALLSTAT(int, Ranks, Ranks_sz);
DYNALLSTAT(int, orbits, orbits_sz);
#endif

DEFAULTOPTIONS_TRACES(traces_opts);
TracesStats traces_stats;

void CreateRandColors(int n)
{
	int i=0, ind0, ind1, ind2;
    char aux[6];
    
    ran_init(n);
	for (i = 0; i<n; i++)
	{
		strcpy(NShape[i].color, "\0");
		strcpy(aux, "\0");

        ind1 = KRAN(4096);
        
	    NShape[i].labcol=(ind1/128) % 2;
	    for (ind0 = 0; ind0<3; ind0++)
		{
            ind2 = ind1 % 16;
            ind1 = (ind1-ind2) / 16;
            switch (ind2)
            {
				case 0:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "00");
                    strcat(NShape[i].color, aux);
                    break;
				case 1:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "11");
                    strcat(NShape[i].color, aux);
                    break;
				case 2:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "22");
                    strcat(NShape[i].color, aux);
                    break;
				case 3:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "33");
                    strcat(NShape[i].color, aux);
                    break;
				case 4:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "44");
                    strcat(NShape[i].color, aux);
                    break;
				case 5:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "55");
                    strcat(NShape[i].color, aux);
                    break;
				case 6:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "66");
                    strcat(NShape[i].color, aux);
                    break;
				case 7:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "77");
                    strcat(NShape[i].color, aux);
                    break;
				case 8:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "88");
                    strcat(NShape[i].color, aux);
                    break;
				case 9:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "99");
                    strcat(NShape[i].color, aux);
                    break;
				case 10:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "AA");
                    strcat(NShape[i].color, aux);
                    break;
				case 11:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "BB");
                    strcat(NShape[i].color, aux);
                    break;
				case 12:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "CC");
                    strcat(NShape[i].color, aux);
                    break;
				case 13:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "DD");
                    strcat(NShape[i].color, aux);
                    break;
				case 14:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "EE");
                    strcat(NShape[i].color, aux);
                    break;
				case 15:
                    strcpy(aux, NShape[i].color);
                    strcpy(NShape[i].color, "FF");
                    strcat(NShape[i].color, aux);
                    break;
            }
		}
	}
}

double ComputeFontsize(int vtx) {
    int ndigits;
    double fsize;
    
    if (vtx == 0) {
        fsize = 13.0;
    }
    else {
        ndigits=log10(vtx)+1;
        switch (ndigits) {
            case 1:
                fsize = 13.0;
                break;
            case 2:
                fsize = 13.0;
                break;
            case 3:
                fsize = 11;
                break;
            case 4:
                fsize = 9.5;
                break;
            case 5:
                fsize = 8;
                break;
            case 6:
                fsize = 6.5;
                break;
            default:
                fsize = 5.0;
                break;
        }
    }
    return fsize;
}

int
main(int argc, char *argv[])
{
	int m, n, c, a, i, j, k, end, n_count;
    size_t j1, e_count;
	int argnum, initorg, cell, numCol, modcode, vtx, refcode, MaxV, MaxE;
	char *arg, sw;
	boolean badargs, make_ranks, too_big;
	boolean rswitch, dswitch, mswitch, nswitch, oswitch, iswitch, xswitch, Sswitch, gswitch, fswitch, vswitch, rswitch1, Eswitch, Vswitch;
	char *infilename, *outfilename, *drefilename;
	FILE *infile, *outfile, *drefile;
	nauty_counter nin;
	char s[10];
	char model[10];
    int numcells, flind, indivtx;
    long minil, maxil;
    double hsize, vsize, fsize;
    sparsegraph g, g1;
    list *liststart, *listend;
    int nodescale, distance, StInd, StIndDist, RnkInd, RnkIndDist;
    
	HELP; PUTVERSION;
    
    liststart = listend = NULL;
	rswitch = dswitch = mswitch = nswitch = oswitch = iswitch = xswitch = Sswitch = gswitch = fswitch = vswitch = rswitch1 = Vswitch = Eswitch = FALSE;
	infilename = outfilename = drefilename = NULL;
	initorg = minil = maxil = flind = 0;
    hsize = 10.00;
    vsize = 6.18;
    modcode = 0;
	n = -1;
    nodescale = 1;
    distance = 0;
    
	argnum = 0;
	badargs = make_ranks = FALSE;
    MaxV = 1000;
    MaxE = 5000;
    
	for (j = 1; !badargs && j < argc; ++j)
	{
	    arg = argv[j];
	    if (arg[0] == '-' && arg[1] != '\0')
	    {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                SWBOOLEAN('i',iswitch)
                else SWBOOLEAN('x',xswitch)
                else SWBOOLEAN('g',gswitch)
                else SWBOOLEAN('v',vswitch)
                else SWINT('V', Vswitch, MaxV, ">E dretodot -V")
                else SWINT('E', Eswitch, MaxE, ">E dretodot -E")
                else SWINT('d', dswitch, distance, ">E dretodot -d")
                else SWINT('F', fswitch, indivtx, ">E dretodot -F")
                else SWINT('m', mswitch, modcode, ">E dretodot -m")
                else SWINT('n', nswitch, nodescale, ">E dretodot -n")
                else SWINT('o', oswitch, initorg, ">E dretodot -o")
                else SWRANGE('r',":-",rswitch,minil,maxil,">E dretodot -r")
                else SWREALRANGE('S',":-",Sswitch,hsize,vsize,">E dretodot -S")
                else badargs = TRUE;
            }
            if (rswitch) {
                rswitch1 = TRUE;
                rswitch = FALSE;
                for (i=minil; i<=maxil; i++) {
                    if (liststart == NULL) {
                        liststart = NewListelem(i);
                        listend = liststart;
                    } else {
                        listend->next = NewListelem(i);
                        listend = listend->next;
                    }
                }
            }
	    }
	    else
	    {
            ++argnum;
            if      (argnum == 1) infilename = arg;
	        else if (argnum == 2) outfilename = arg;
	        else if (argnum == 3) drefilename = arg;
            else                  badargs = TRUE;
	    }
	}
    
	if (labelorg < 0) gt_abort(">E dretodot: negative origin forbidden\n");
    
	if (badargs || argnum > 3)
	{
	    fprintf(stderr, ">E Usage: %s\n", USAGE);
	    GETHELP;
	    exit(1);
	}
    
	if (!infilename || infilename[0] == '-')
	{
	    infilename = "stdin";
	    infile = stdin;
	}
	else if ((infile = fopen(infilename, "r")) == NULL)
	{
	    fprintf(stderr, "Can't open input file %s\n", infilename);
	    gt_abort(NULL);
	}
    
	if (!outfilename || outfilename[0] == '-')
	{
	    outfilename = "stdout";
	    outfile = stdout;
	}
	else if ((outfile = fopen(outfilename, "w")) == NULL)
	{
	    fprintf(stderr, "Can't open output file %s\n", outfilename);
	    gt_abort(NULL);
	}

	labelorg = initorg;
	nin = 0;
    
	while (fscanf(infile, "%1s", s) == 1)
	{
	    if (s[0] == 'n')
	    {
            if (fscanf(infile, "%1s", s) == 1 && s[0] != '=')
                ungetc(s[0], infile);
            if (fscanf(infile, "%d", &n) != 1)
            {
                fprintf(stderr, ">E dretodot: invalid n=# command\n");
                gt_abort(NULL);
            }
            if (n <= 0)
                gt_abort(">E dretodot: n can't be <= 0\n");
	    }
	    else if (s[0] == '"')
	    {
            while ((c = getc(infile)) != '"' && c != EOF) {}
	    }
	    else if (s[0] == '!')
	    {
            while ((c = getc(infile)) != '\n' && c != EOF) {}
	    }
	    else if (s[0] == '$')
	    {
            if ((s[0] = getc(infile)) == '$')
                labelorg = initorg;
            else
            {
                if (s[0] != '=') ungetc(s[0], infile);
                if (fscanf(infile, "%d", &labelorg) != 1)
                    gt_abort(">E dretodot: invalid $=# command\n");
                if (labelorg < 0)
                    gt_abort(">E dretodot: must have labelorg >= 0\n");
            }
        }
	    else if (s[0] == 'g' || (s[0] >= '0' && s[0] <= '9')
                 || s[0] == ';')
	    {
            if (n < 0)
                gt_abort(">E dretodot: g command before n is defined\n");
            if (s[0] != 'g') ungetc(s[0], infile);
            m = (n + WORDSIZE - 1) / WORDSIZE;
#if MAXN
            if (n > MAXN || m > MAXM)
                gt_abort(">E n or m too big\n");
#else
            DYNALLOC2(int, lab, lab_sz, n, m, "dretodot");
            DYNALLOC2(int, invlab, invlab_sz, n, m, "dretodot");
            DYNALLOC2(int, ptn, ptn_sz, n, m, "dretodot");
            DYNALLOC2(NodeShape, NShape, NShape_sz, n, m, "dretodot");
            DYNALLOC2(int, orbits, orbits_sz, n, m, "dretodot");
#endif
            
            for (i=0; i<n; i++) {
                lab[i] = i;
                ptn[i] = 1;
            }
            ptn[n-1] = 0;
            numcells = 1;
            ++nin;
            SG_INIT(g);
            SG_INIT(g1);
            readgraph_sg(infile, &g, FALSE, FALSE, 0, n);
	    }
        else if (s[0] == 'f') {
            readptn(infile, lab, ptn, &numcells, FALSE, n);
            if (numcells > 1) {
                traces_opts.defaultptn = FALSE;
            }
        }
        else
	    {
            fprintf(stderr, ">E dretodot: invalid command \"%c\"\n", s[0]);
            gt_abort(NULL);
	    }
	}

    fclose(infile);
    
    if ((nodescale != 0) && (nodescale != 1) && (nodescale != 2)) {
        nodescale = 1;
    }

    if (!drefilename || drefilename[0] == '-')
	{
	    drefilename = "stdout";
	    drefile = stdout;
	}
	else if ((drefile = fopen(drefilename, "w")) == NULL)
	{
	    fprintf(stderr, "Can't open output file %s\n", drefilename);
	    gt_abort(NULL);
	}
    
    if (!dswitch) distance = n;
    
    if (iswitch && !xswitch && !fswitch) {
        if (vswitch) fprintf(stderr, ">Z  partition refinement...");
        refine_tr(&g,lab,ptn,&numcells,&refcode,&traces_opts);
        if (vswitch) fprintf(stderr, "done.\n");
        WRITEDREFILE
    }
    
    if (fswitch) {
        indivtx -= labelorg;
        refine_tr(&g,lab,ptn,&numcells,&refcode,&traces_opts);
        cell = 0;
        for (i=0; i<n; i++) {
            if (lab[i] == indivtx) break;
            if (!ptn[i]) {
                cell = i+1;
            }
        }
        lab[i] = lab[cell];
        lab[cell] = indivtx;
        ptn[cell] = 0;
        numcells++;
        traces_opts.defaultptn = FALSE;
        refine_tr(&g,lab,ptn,&numcells,&refcode,&traces_opts);
        if (!xswitch) WRITEDREFILE
    }
    
    if (xswitch) {
        if (vswitch) fprintf(stderr, ">Z  running Traces...");
        Traces(&g,lab,ptn,orbits,&traces_opts,&traces_stats,&g1);
        if (vswitch) fprintf(stderr, "done.\n");
        WRITEDREFILE
        for (i=0; i<n; i++) {
            lab[i] = i;
        }
        sort2ints(orbits,lab,n);
        memcpy(ptn,orbits,n*sizeof(int));
    }
    
    strcpy(model, "\0");
    switch (modcode) {
        case 1:
            strcpy(model, "neato");
            break;
        case 2:
            strcpy(model, "fdp");
            break;
        case 3:
            strcpy(model, "sfdp");
            break;
        case 4:
            strcpy(model, "twopi");
            break;
        case 5:
            strcpy(model, "circo");
            break;
        default:
            strcpy(model, "dot");
            modcode = 0;
            break;
    }
    fprintf(outfile, "graph G\n{\n");
    
    if (vswitch) fprintf(stderr, ">Z  selected drawing model \"%s\".\n",model);
    
#if !MAXN
    DYNALLOC2(int, CurrVertices, CurrVertices_sz, n, m, "dretodot");
    DYNALLOC2(int, DistStack, DistStack_sz, n, m, "dretodot");
    DYNALLOC2(int, Ranks, Ranks_sz, n, m, "dretodot");
#endif

    if (vswitch) fprintf(stderr, ">Z  analyzing graph: step 1");

    if (!liststart) {
        for (i=0; i<n; i++) {
            CurrVertices[i] = 1;
            DistStack[i] = i;
            StInd = n;
        }
    } else {
        memset(CurrVertices, 0, n*sizeof(int));
        StInd = 0;
        listend = liststart;
        while (listend) {
            if (listend->vtx-labelorg < 0 || listend->vtx-labelorg >= n) {
                fprintf(stderr, ">W dretodot: invalid vertex %d\n", listend->vtx);
            } else {
                DistStack[StInd++] = listend->vtx-labelorg;
                CurrVertices[listend->vtx-labelorg] = 1;
            }
            listend = listend->next;
        }
    }
    
    if (vswitch) fprintf(stderr, ", 2");

    cell = 0;
    for (i=0; i<n; i++) {
        invlab[lab[i]] = i;
        if (!xswitch) {
            if (ptn[i] == 0) {
                ptn[i] = cell++;
            } else {
                ptn[i] = cell;
            }
        }
    }

    if (vswitch) fprintf(stderr, ", 3");

    k = RnkInd = 0;
    StIndDist = StInd;
    while (k<StInd) {
        end = StInd;
        for (j=k; j<end; j++) {
            vtx = DistStack[j];
                for (j1 = g.v[vtx]; j1<g.v[vtx]+g.d[vtx]; j1++) {
                    if (!CurrVertices[g.e[j1]]) {
                        if (RnkInd<distance) {
                            CurrVertices[g.e[j1]] = 1;
                        } else {
                            CurrVertices[g.e[j1]] = 2;
                        }
                        DistStack[StInd++] = g.e[j1];
                    }
                }
        }
        if (RnkInd<distance) {
            StIndDist = StInd;
            RnkIndDist = RnkInd;
        }
        Ranks[RnkInd++] = k;
        k=end;
    }

    if (gswitch) {
        n_count = n;
        e_count = g.nde;
    } else {
        n_count = 0;
        e_count = 0;
        for (i=0; i<n; i++) {
            if (CurrVertices[i] == 1) {
                n_count++;
                for (j1 = g.v[i]; j1<g.v[i]+g.d[i]; j1++) {
                    if ((i < g.e[j1]) && (CurrVertices[g.e[j1]] == 1)) {
                        e_count++;
                    }
                }
            }
        }
    }
    
    too_big = FALSE;
    if (n_count > MaxV) {
        fprintf(stderr, ">E Too many vertices (%d, max: %d; use -V# (at your own risk))\n", n_count, MaxV);
        too_big = TRUE;
    }
    if (e_count > MaxE) {
        fprintf(stderr, ">E Too many edges (%lu, max: %d; use -E# (at your own risk))\n", e_count, MaxE);
        too_big = TRUE;
    }
    
    if (too_big) {
        exit(1);
    }
    if (vswitch) fprintf(stderr, ", 4");

    if (Sswitch) {
        if (gswitch) {
            hsize *= log10(n);
            vsize *= log10(n);
        } else {
            if (StIndDist) {
                hsize *= log10(StIndDist);
                vsize *= log10(StIndDist);
            }
        }
    }

    if (vswitch) fprintf(stderr, ".\n");

    if (vswitch) fprintf(stderr, ">Z  creating colors...");
    CreateRandColors(n);
    if (vswitch) fprintf(stderr, "done\n");
    
    if (vswitch) fprintf(stderr, ">Z  drawing size %.2f x %.2f\n", hsize, vsize);
    
    /* Preamble */
	fprintf(outfile, "graph [center=\"true\", size=\"%.2f, %.2f\", ratio=\"fill\", ranksep=\"0.25\", nodesep=\"0.40\",\n", hsize, vsize);
    fprintf(outfile, "       outputorder=\"edgesfirst\", overlap=\"scale\", layout=\"%s\"];\n",model);
	fprintf(outfile, "node  [shape=\"circle\", width=\"%.2f\", height=\"%.2f\", fixedsize=\"true\",\n",
            0.5*nodescale, 0.5*nodescale);
    fprintf(outfile, "       style=\"filled\", color=\"black\",\n");
    fprintf(outfile, "       fontsize=\"13\", fontname=\"%s\"];\n", FONT_NAME);
    
    /* Nodes */
    if (vswitch) fprintf(stderr, ">Z  drawing vertices...");

    fprintf(outfile, "node  [color=\"black\", fontcolor=\"black\"]\n");

    a = 1;
    if (dswitch) end = StIndDist; else end = StInd;
    make_ranks = ((modcode == 0) || (modcode > 5)) && rswitch1;

    if (make_ranks) fprintf(outfile, "{ rank=\"source\";\n");
    for (i=0; i<end; i++) {
        vtx = DistStack[i];
        fsize = ComputeFontsize(vtx) * (double)nodescale * 1.666;
        if (make_ranks && (i==Ranks[a])) {
            fprintf(outfile, "}\n{ rank=\"same\";\n");
            a++;
        }
        numCol = ptn[invlab[vtx]];
        if (fsize == 13.0) {
            if (NShape[numCol].labcol == 0)
                fprintf(outfile, "%6d [fontcolor=\"white\", fillcolor=\"#%s\"]\n", vtx + labelorg,
                        NShape[numCol].color);
            else
                fprintf(outfile, "%6d [fillcolor=\"#%s\"]\n", vtx + labelorg,
                        NShape[numCol].color);
        } else {
            if (NShape[numCol].labcol == 0)
                fprintf(outfile, "%6d [fontsize=\"%.2f\", fontcolor=\"white\", fillcolor=\"#%s\"]\n", vtx + labelorg,
                        fsize, NShape[numCol].color);
            else
                fprintf(outfile, "%6d [fontsize=\"%.2f\", fillcolor=\"#%s\"]\n", vtx + labelorg,
                        fsize, NShape[numCol].color);
        }
    }
    if (make_ranks) fprintf(outfile, "}\n");
    
    if (gswitch) {
        fprintf(outfile, "node  [fontsize=\"%.2f\", fontcolor=\"gray\", color=\"gray\", fillcolor=\"gray97\"]\n", fsize);
        for (i=0; i<n; i++) {
            if (CurrVertices[i] != 1) {
                fsize = ComputeFontsize(i) * (double)nodescale * 1.666;
                fprintf(outfile, "%6d\n", i+labelorg);
            }
        }
    }

    if (vswitch) fprintf(stderr, "done\n");

    /* Edges */
    if (vswitch) fprintf(stderr, ">Z  drawing edges...");

    if (gswitch) {
        fprintf(outfile, "edge  [penwidth=\"0.4\", color=\"gray\", weight=\"10\"];\n");
        for (vtx = 0; vtx < n; vtx++) {
            for (j1 = g.v[vtx]; j1 < g.v[vtx] + g.d[vtx]; ++j1) {
                if (CurrVertices[vtx] != 1 || CurrVertices[g.e[j1]] != 1) {
                    if (g.e[j1] > vtx)
                        fprintf(outfile, "%d -- %d;\n", vtx + labelorg, g.e[j1] + labelorg);
                }
            }
        }
    }

    fprintf(outfile, "edge  [penwidth=\"0.8\", color=\"black\", weight=\"10\"];\n");
	for (i = 0; i < StIndDist; ++i) {
        vtx = DistStack[i];
		for (j1 = g.v[vtx]; j1 < g.v[vtx] + g.d[vtx]; ++j1) {
            if (CurrVertices[g.e[j1]] == 1) {
                if (g.e[j1] > vtx)
                    fprintf(outfile, "%d -- %d;\n", vtx + labelorg, g.e[j1] + labelorg);
            }
        }
    }
    
    if (gswitch) {
        RnkIndDist = RnkInd;
    }
    if (RnkIndDist>1) {
        fprintf(outfile, "%d ", DistStack[0] + labelorg);
        for (i=1; i<RnkIndDist; i++) {
            fprintf(outfile, "-- %d ", DistStack[Ranks[i]] + labelorg);
        }
        fprintf(outfile, "[style=\"invis\"];\n");
    }

    fprintf(outfile, "}\n");

    if (vswitch) fprintf(stderr, "done\n>Z  " COUNTER_FMT
            " graph drawn from %s to %s\n",
            nin, infilename, outfilename);
	exit(0);
}
