/*****************************************************************************
*                                                                            *
*   Quarticgen by Brendan D. McKay and S. Narjess Afzaly                     *
*   This program generates all quartic graphs up to a given order            *
*                                                                            *
*   Parameters - <r> means read-only, <w> means write-only, <wr> means both: *
*          nmax  <r>  # vertices.  This must be at least 5 and               *
*                  at most MAXN.                                             *
*                                                                            *
*   numread     <w> # irreducible graphs with n <= nmax vertices             *
*   numwritten  <w> # quartic graphs of order nmax                           *
*   diffnum     <w> # all non-root (not irreducible ones) nodes in the       *
*   search tree (both accepted and rejected ones)                            *
*   numaac      <w> # accepted non-root nodes in the search tree             *
*                                                                            *
 ****************************************************************************/
//Just the nodes above or on the SL will be generated. The children of nodes oh the SL won't be investigated or generated.
#define MAXN 28
#define LEV1 2         //if (nmax<= threshold) then (SplitLevel = nmax-LEV1 ), otherwise (SplitLevel = LEV2).
#define LEV2 15
#define threshold 16
#define WORDSIZE 64
#include "gtools.h"
#include "quarticirred28.h"


#define MAXE (2*MAXN)
#define MAXP ((MAXE*MAXE))

static void (*outproc)(FILE*, graph*, int);
#ifdef OUTPROC
extern void OUTPROC(FILE*, graph*, int);
#endif

static DEFAULTOPTIONS_GRAPH(options);
static int perm[MAXN];
static setword workspace[50];
static statsblk(stats);





// badan -c#: (-c1 --> only connected ones ) (-c2 --> only biconnected ones)
#define USAGE "genquarticg [-ugs -h -c -l] n [res/mod] [file]"


#define HELPTEXT \
"  generate all non-isomorphic quartic graphs of a given order \n\
\n\
n     : the number of the vertices\n\
file  : the name of the output file (default stdout)\n\
-u    : do not output any graphs, just generate and count them\n\
-g    : use graph6 format for output (default)\n\
-s    : use sparse6 format for output\n\
-h      write a header (only with -g or -s). \n\
-c    : only write connected graphs\n\
-C    : only write biconnected graphs\n\
res/mod : only generate subset res out of subsets 0..mod-1\n\
-l    : canonically label output graphs.\n" //??????? is it a good feature to add?



static FILE *outfile, *msgfile;             /* file for output graphs */
static char *outfilename;
boolean     nooutput;                       /* presence of -u */
boolean     graph6;                         /* presence of -g */
boolean     sparse6;                        /* presence of -s */
boolean     header;                         /* presence of -h */
static int  connec;                          /* 1 for -c, 2 for -C, 0 for neither */
boolean     connec1;                        /* presence of -c */
boolean     connec2;                        /* presence of -C */
boolean     canonise;                       /* presence of -l */  //????????


typedef struct
{
	int first;
    int sec;
    setword fn;
    setword sn;
    setword end;
    int intersect;
    int cond;
} edgestruct;

typedef struct
{
	int first;
    int sec;
    int multp;
} pairstruct;

typedef struct
{
	int base;
	int first1;
	int sec1;
	int first2;
	int sec2;
} dovistruct;

typedef   enum
{ accept,reject,undef } CHOISE;

int         *pCNT, *pnumpair, *pdoviorbit, *pepairorbit;
int         (*pantidovi)[MAXN], (*pantipair)[MAXE], (*pantiedge)[MAXN];
int         count[MAXN];
long        numread;
nauty_counter numwritten;
pairstruct  *pepair;
dovistruct  *pdovi;
edgestruct  *pedge;
setword     active;
boolean     goodret;
static int  nmax, m, code, numcells, mod, res, splitlevel, splitcount;
//TMP static splitlevel,splitcount,mod,res;    ??????TMP?
static void    extend(int, graph *, edgestruct *, pairstruct *, int, int * , int *, setword *, int *, boolean );
static int     init_refinex( int *, int *, int *, set *, int);
static void     refinex( graph *, int *, int *, int , int *, int *, set *, boolean , int *, int , int );
static UPROC    userautom1();
static UPROC 	userautom2();
static UPROC 	userautom3();


/************************************************************************/

static void
writeg6x(FILE *f, graph *g, int n)
/* write graph g with n vertices to file f in graph6 format */
{
    writeg6(f, g, 1, n);
}

/************************************************************************/

static void
writes6x(FILE *f, graph *g, int n)
/* write graph g with n vertices to file f in graph6 format */
{
    writes6(f, g, 1, n);
}

/***********************************************************************/

static void
nullwrite(FILE *f, graph *g, int n)
/* don't write graph g to file f */
{
}

/***********************************************************************/

static boolean
isconnected(graph *g, int n)
/* test if g is connected with vertices 0,1.., n-1 */
{
    setword seen,expanded,toexpand,allbits;
    int i;

    allbits = ALLMASK(n);

    expanded = bit[n-1];
    seen = expanded | g[n-1];

    while (seen != allbits && (toexpand = (seen & ~expanded))) /* not == */
    {
        i = FIRSTBITNZ(toexpand);
        expanded |= bit[i];
        seen |= g[i];
    }

    return  seen == allbits;
}

/**********************************************************************/

static boolean
isbiconnected(graph *g, int n)
/* test if g is biconnected */
{
    int sp,v,w;
    setword sw;
    setword visited;
    int numvis,num[MAXN],lp[MAXN],stack[MAXN];

    if (n <= 2) return FALSE;

    visited = bit[0];
    stack[0] = 0;
    num[0] = 0;
    lp[0] = 0;
    numvis = 1;
    sp = 0;
    v = 0;

    for (;;)
    {
        if ((sw = g[v] & ~visited))           /* not "==" */
        {
            w = v;
            v = FIRSTBITNZ(sw);       /* visit next child */
            stack[++sp] = v;
            visited |= bit[v];
            lp[v] = num[v] = numvis++;
            sw = g[v] & visited & ~bit[w];
            while (sw)
            {
                w = FIRSTBITNZ(sw);
                sw &= ~bit[w];
                if (num[w] < lp[v])  lp[v] = num[w];
            }
        }
        else
        {
            w = v;                  /* back up to parent */
            if (sp <= 1)          return numvis == n;
            v = stack[--sp];
            if (lp[w] >= num[v])  return FALSE;
            if (lp[w] < lp[v])    lp[v] = lp[w];
        }
    }
}


/****************************************************************************
 *    extend recursively extends quartic graph until number of vertices      *
 *    of the graph reaches nmax after which it is written to outputfile      *
 ****************************************************************************/

static void
extend(int n, graph *g, edgestruct *edge, pairstruct *epair, int numpair,
   int *epairorbit, int *multar, setword *zar, int *col00w, boolean connectflag)
{
    int   vm1, vm2, vm3, vm4, vt1, vt2, vt3, vt4, c, b, mcol1, mcol,
          tcol, got_one, i, j, k, numpair1, numdovi, maxdovi, i1, j1, i2, j2,
          numrival, temp, mult, multm, mult2, rely, numedge, dcol, dcolp, e1, e2;
    int   firsttime[MAXN], firsttimey[MAXN], mult2i[MAXN], multar1[MAXN],
          col00[MAXN], col00w1[MAXN], rival[MAXN], tagriv[MAXN], doviorbit[3*MAXN],
          epairorbit1[MAXP], l[MAXN], lab[MAXN], ptn[MAXN], orbits[MAXN];
    int 		antipair[MAXE][MAXE], antiedge[MAXN][MAXN], antidovi[MAXN][MAXN];
    setword     x, y, z, x1, bitj1, bitj2, biti2;
    setword     yi[MAXN], zar1[MAXN];
    CHOISE	colours; 
    pairstruct 	epair1[MAXP];
    edgestruct 	edge1[MAXE];
    graph       gi[MAXN];
    dovistruct 	dovi[3*MAXN];
    dovistruct  dovimax;
    register    setword    gi1, gi2, gj1, gj2;
    boolean     conf;

    /////////////////////////////////////////////////////////////////////////
    if( n == splitlevel )
    {
        if (splitcount-- != 0) return;
        splitcount = mod - 1;
    }
    /////////////////////////////////////////////////////////////////////////
    for( c = 1; c < numpair; c++)
    {
        if( epairorbit[c] == c  )
        {

            int eqn[MAXN],neqn;
            setword zval[MAXN],zz,eqcol0;
            int dcol0,dcol1;

            e1 = epair[c].first;
            e2 = epair[c].sec;

            vm1 = edge[e1].first;
            vm2 = edge[e1].sec;
            vm3 = edge[e2].first;
            vm4 = edge[e2].sec;

            g[vm1] ^= bit[n] | bit[vm2];
            g[vm2] ^= bit[n] | bit[vm1];
            g[vm3] ^= bit[n] | bit[vm4];
            g[vm4] ^= bit[n] | bit[vm3];

            g[n] = bit[vm1] | bit[vm2] | bit[vm3] | bit[vm4];

            multar1[n] = multm = epair[c].multp;

            got_one = 0;
            dcol0 = ( (g[vm1]&g[vm2])!= bit[n]) + ((g[vm3]&g[vm4])!=bit[n]);


            zar1[n] = z = g[vm1] ^ g[vm2] ^ g[vm3] ^ g[vm4];
            col00w1[n] = col00[n] = -POPCOUNT(z) + ( (4-multm) << 6 );

            zval[n] = z; eqcol0 = bit[n]; neqn = 0;
            setword xn;
            xn = g[n]|bit[n];


            b = n-1;		
            while( (!got_one)  &&   (b >= 0) )
            {					


               if(  !((g[b]|bit[b]) & xn )   &&  (multar[b] < 9 ) )
                {
                    mult = multar1[b] = multar[b];
                    if( mult )
                    {
                        z = zar1[b] = zar[b];
                        col00[b] = col00w1[b] = col00w[b];
                    }

                }
                else
                {
                    x = g[b];
                    bitj2 = x&(-x); x ^= bitj2; j2 = FIRSTBITNZ(bitj2);
                    biti2 = x&(-x); x ^= biti2; i2 = FIRSTBITNZ(biti2);
                    bitj1 = x&(-x); x ^= bitj1; j1 = FIRSTBITNZ(bitj1);
                    i1 = FIRSTBITNZ(x);
                    mult = 0;
                    if(   (!(g[i1] & bitj1))  &&   (!(g[i2] & bitj2))  )
                        mult++;
                    if(   (!(g[i1] & bitj2))  &&   (!(g[i2] & bitj1))  )
                        mult++;
                    if(   (!(g[i1] & biti2))  &&   (!(g[j1] & bitj2))  )
                        mult++;
                    multar1[b] = mult;


                    if( mult )
                    {

                        z = zar1[b] = g[i1] ^ g[i2] ^ g[j1] ^ g[j2];
                        col00[b] = col00w1[b] = -POPCOUNT(z) + ( (4-mult) << 6 );
                    }
                }


                if( mult )
                {
                    if( col00[b] > col00[n] )
                            got_one = 1;
                    else if ( col00[b] == col00[n] )
                    {
                            eqcol0 |= bit[b];
                            zval[b] = z;
                            eqn[neqn++] = b;
                    }
                }
                else
                        col00w1[b] = col00[b] = 0;	
                b--; 	
            }




            if (!got_one && neqn > 0)
            {
                int col00x,col00y,v;
                setword yz;

                yz = eqcol0 & zval[n];
                col00x = MAXN-POPCOUNT(yz);
                col00[n] += col00x;
                while (--neqn >= 0)
                {
                    v = eqn[neqn];
                    yz = eqcol0 & zval[v];
                    col00y = MAXN-POPCOUNT(yz);
                    if (col00y > col00x)
                    {
                        got_one = 1;
                        break;
                    }
                    else col00[v] += col00y;
                }
            }

            if (  !got_one	) 	
            {
                numcells = init_refinex( col00, lab, ptn, &active, n);		

                refinex(g, lab, ptn, 0, &numcells, count, &active, TRUE, &code, 1, n+1);	

                if( code == -1 )
                    got_one = 1;

                if (	  !got_one	)
                    mcol1 =  ((col00[vm1] + col00[vm2] ) * (col00[vm3] + col00[vm4]));
            }


            numdovi = 0;
            j = n;				
            while(  (!got_one)       &&       (  (j==n)  || (ptn[j])  )     )
            {			  			
                b = lab[j];
                if( b != n)
                {				
                    x = g[b];
                    TAKEBIT(i1,x);				
                    TAKEBIT(j1,x);
                    TAKEBIT(i2,x);
                    TAKEBIT(j2,x);  			
                    for(i=0; i < 3 && !got_one; i++)
                    {
                        if(i == 0)
                        {
                            vt1 = i1;
                            vt2 = j1;
                            vt3 = i2;
                            vt4 = j2;
                        }

                        if(i == 1)
                        {
                            vt1 = i1;
                            vt2 = i2;
                            vt3 = j1;
                            vt4 = j2;
                        }

                        if(i == 2)
                        {
                            vt1 = i1;
                            vt2 = j2;
                            vt3 = j1;
                            vt4 = i2;
                        }

                        if( vt3 > vt4 )
                        {
                            temp = vt3;
                            vt3 = vt4;
                            vt4 = temp;
                        }						

                        if( (!(g[vt1] & bit[vt2]))  &&   (!(g[vt3] & bit[vt4]))   &&  (!got_one) )
                        {	
                            tcol =  ((col00[vt1] + col00[vt2] ) * (col00[vt3] + col00[vt4]));
                            dcol1 = ((g[vt1]&g[vt2])!=bit[b]) + ((g[vt3]&g[vt4])!=bit[b]);
                            if( dcol1 < dcol0 || (dcol1 == dcol0 && tcol < mcol1 ))
                                got_one = 1;
                            else if( dcol1 == dcol0 && tcol == mcol1 )
                            {

                                    dovi[numdovi].base = b;
                                    dovi[numdovi].first1 = vt1;
                                    dovi[numdovi].sec1 = vt2;
                                    dovi[numdovi].first2 = vt3;
                                    dovi[numdovi].sec2 = vt4;
                                    antidovi[b][vt2] =numdovi;
                                    numdovi++;				
                            }	
                        }
                    }							
                }

                j--;
                if(j<0)
                    break; 	
            }



            if( (!(g[vm1] & bit[vm3]))  &&   (!(g[vm2] & bit[vm4]))   &&  (!got_one) )
            {   	
                if( vm2 < vm4 )
                {
                    tcol =  ((col00[vm1] + col00[vm3] ) * (col00[vm2] + col00[vm4]));
                    dcol1 = ((g[vm1]&g[vm3])!=bit[n]) + ((g[vm2]&g[vm4])!=bit[n]);
                    if (dcol1 < dcol0 || (dcol1 == dcol0 &&  tcol < mcol1 ))
                        got_one = 1;
                    else if( dcol1 == dcol0 && tcol == mcol1 )
                    {
                            dovi[numdovi].base = n;
                            dovi[numdovi].first1 = vm1;
                            dovi[numdovi].sec1 = vm3;
                            dovi[numdovi].first2 = vm2;
                            dovi[numdovi].sec2 = vm4;
                            antidovi[n][vm3] =numdovi;
                            numdovi++;  			
                    }

                }
                else
                {
                    tcol = ((col00[vm1] + col00[vm3] ) * (col00[vm4] + col00[vm2]));

                    dcol1 = ((g[vm1]&g[vm3])!=bit[n]) + ((g[vm2]&g[vm4])!=bit[n]);
                    if (dcol1 < dcol0 || (dcol1 == dcol0 &&  tcol < mcol1 ))
                        got_one = 1; 				
                    else if( dcol1 == dcol0 && tcol == mcol1 )
                    {
                            dovi[numdovi].base = n;
                            dovi[numdovi].first1 = vm1;
                            dovi[numdovi].sec1 = vm3;
                            dovi[numdovi].first2 = vm4;
                            dovi[numdovi].sec2 = vm2;
                            antidovi[n][vm3] =numdovi;
                            numdovi++;								
                    }
                }	
            }

            if( (!(g[vm1] & bit[vm4]))  &&   (!(g[vm2] & bit[vm3]))   &&  (!got_one) )
            {   	
                if( vm2 < vm3 )
                {
                    tcol = ((col00[vm1] + col00[vm4] ) * (col00[vm2] + col00[vm3]));
                    dcol1 = ((g[vm1]&g[vm4])!=bit[n]) + ((g[vm2]&g[vm3])!=bit[n]);
                    if (dcol1 < dcol0 || (dcol1 == dcol0 && tcol < mcol1 ))
                        got_one = 1;
                    if( dcol1 == dcol0 && tcol == mcol1 )
                    {
                            dovi[numdovi].base = n;
                            dovi[numdovi].first1 = vm1;
                            dovi[numdovi].sec1 = vm4;
                            dovi[numdovi].first2 = vm2;
                            dovi[numdovi].sec2 = vm3;
                            antidovi[n][vm4] =numdovi;
                            numdovi++;		
                    }

                }
                else
                {	
                    tcol =  ((col00[vm1] + col00[vm4] ) * (col00[vm3] + col00[vm2]));
                    dcol1 = ((g[vm1]&g[vm4])!=bit[n]) + ((g[vm2]&g[vm3])!=bit[n]);
                    if (dcol1 < dcol0 || (dcol1 == dcol0 && tcol < mcol1 ))
                        got_one = 1;
                    if( dcol1 == dcol0 && tcol == mcol1 )
                    {
                            dovi[numdovi].base = n;
                            dovi[numdovi].first1 = vm1;
                            dovi[numdovi].sec1 = vm4;
                            dovi[numdovi].first2 = vm3;
                            dovi[numdovi].sec2 = vm2;
                            antidovi[n][vm4] =numdovi;
                            numdovi++;       	
                    }
                }
            }					


            if( got_one == 1 )
                colours = reject;
            if(   (!got_one) 	&&	 (!numdovi)   	)
                colours = accept;
            if( (!got_one) && (numdovi) )
            {
                colours = undef;
                dovi[numdovi].base = n;
                dovi[numdovi].first1 = vm1;
                dovi[numdovi].sec1 = vm2;
                dovi[numdovi].first2 = vm3;
                dovi[numdovi].sec2 = vm4;
                antidovi[n][vm2] =numdovi;
                numdovi++;	                              	
            }				

/*****************************************************************************/
            if( n+1 == nmax )
            {
                if( connec1 )
                {
                    conf = TRUE;
                    if( !connectflag )
                        conf = isconnected(g, n+1);
                }
                if(  !connec  ||  (connec1 && conf ) || (connec2  &&  isbiconnected(g, n+1))  )
                {
                    if( colours == accept )
                    {
                        (*outproc)( outfile, g, n+1);
                        ++numwritten;
                    }

                    if (colours == undef)
                    {  	
                        for (i = 0; i < numdovi; i++)
                            doviorbit[i] = i;

                        pdoviorbit = doviorbit;
                        pdovi = dovi;
                        pantidovi = antidovi;
                        pCNT = &numdovi;					

                        options.userautomproc = userautom2;
                        options.writeautoms = FALSE;
                        options.writemarkers = FALSE;
                        options.getcanon = TRUE;
                        options.defaultptn = FALSE;   ///// 1.badan in FALSE SHAVAD, 2.COULR AVALI TARIF SHAD, 3.lab o ptn initialize shavand

                        nauty( g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50, m, n+1, gi);				

                        for (i = 0; i < n+1; i++)
                            l[lab[i]] = i;

                        dovimax.base = -1;
                        dovimax.first1= -1;
                        dovimax.sec1 = -1;
                        dovimax.first2 = -1;
                        dovimax.sec2 = -1;

                        for (i = 0; i < numdovi; i++)
                        {
                            b = l[dovi[i].base];
                            if( b >= dovimax.base )
                            {
                                vt1 = l[dovi[i].first1];
                                vt1 = l[dovi[i].first1];
                                vt2 = l[dovi[i].sec1];
                                vt3 = l[dovi[i].first2];
                                vt4 = l[dovi[i].sec2];


                                if( vt1 > vt2 )
                                {
                                    temp = vt1;
                                    vt1 = vt2;
                                    vt2 = temp;
                                }
                                if( vt3 > vt4 )
                                {
                                    temp = vt3;
                                    vt3 = vt4;
                                    vt4 = temp;
                                }
                                if( vt1 > vt3 )
                                {
                                    temp = vt1;
                                    vt1 = vt3;
                                    vt3 = temp;

                                    temp = vt2;
                                    vt2 = vt4;
                                    vt4 = temp;
                                }

                                if( (b > dovimax.base)  || ( vt1 > dovimax.first1)  ||
                                   (  vt1 == dovimax.first1  &&  vt2 > dovimax.sec1)  ||
                                   (   vt1 == dovimax.first1  &&  vt2 == dovimax.sec1  &&   vt3 > dovimax.first2) ||
                                   ( vt1 == dovimax.first1  &&  vt2 == dovimax.sec1  &&   vt3 == dovimax.first2  &&  vt4 > dovimax.sec2)    )
                                {
                                    dovimax.base = b;
                                    dovimax.first1 = vt1;
                                    dovimax.sec1 = vt2;
                                    dovimax.first2 = vt3;
                                    dovimax.sec2 = vt4;

                                    maxdovi = i;
                                }
                            }

                        }

                        if (doviorbit[maxdovi] == doviorbit[numdovi-1])
                        {
                            (*outproc)( outfile, g, n+1);
                            ++numwritten;
                        }

                    }
                }
            }
/*****************************************************************************/

            else
            {		        						                   								
                if( colours != reject)
                {
                    epairorbit1[0] = 0;
                        numpair1 = 1;

                        setword yn;
                        yn = g[n] | bit[n];

                        numedge = 0;
                        for (i = 0; i < n+1; i++)
                        {
                            x = g[i] & BITMASK(i);
                            while (x)
                            {
                                y = x&(-x);
                                x ^= y;
                                j = FIRSTBITNZ(y);

                                edge1[numedge].first = i;
                                edge1[numedge].sec = j;	
                                edge1[numedge].fn = (g[i] & (~bit[j]) );
                                edge1[numedge].sn = (g[j] & (~bit[i]) );
                                edge1[numedge].intersect = ( (g[i] & g[j]) != 0 ) ;
                                edge1[numedge].cond = !( (bit[i]&yn) && (bit[j]&yn ) );
                                antiedge[i][j] = numedge;
                                numedge++;
                            }
                            firsttimey[i] = firsttime[i] = 1;
                        }

                        for (e1 = 0; e1 < numedge; e1++)
                        {						
                            i1 = edge1[e1].first;
                            j1 = edge1[e1].sec;
                            gi1 = edge1[e1].fn;
                            gj1 = edge1[e1].sn;
                            dcolp = edge1[e1].intersect;

                            for (e2 = e1+1; e2 < numedge; e2++)
                            {
                                antipair[e1][e2] = 0;
                                // while( (i1 == edge1[e2].first)    &&    (e2 < numedge) )
                                while( e2 < numedge   &&   (i1 == edge1[e2].first) )
                                {
                                    e2++;
                                    antipair[e1][e2] = 0;
                                }

                                if( e2 < numedge )
                                    if( ( j1 != edge1[e2].first )    &&    ( j1 != edge1[e2].sec )  )
                                    {
                                        i2 = edge1[e2].first;
                                        j2 = edge1[e2].sec;

                                        gi2 = edge1[e2].fn;
                                        gj2 = edge1[e2].sn;

                                        dcol = dcolp + edge1[e2].intersect;



                                        mult = 1;
                                        if(   (!(g[i1] & bit[j2]))  &&   (!(g[i2] & bit[j1]))  )
                                        {
                                            mult++;
                                            if(  ( ( ( gi1 & gj2 )!=0) + (( gi2 & gj1 )!=0) )  < dcol   )
                                                continue;
                                        }

                                        if(   (!(g[i1] & bit[i2]))  &&   (!(g[j1] & bit[j2]))  )
                                        {
                                            mult++;
                                            if(  ( ( ( gi1 & gi2 )!=0) + (( gj1 & gj2 )!=0) )  < dcol   )
                                                continue;
                                        }


                                        rely = got_one = 0;						
                                        if( mult > multm  )
                                        {
                                            if(  edge1[e1].cond   &&   edge1[e2].cond )
                                               continue;

                                            else for(i=1; i <= n/3    &&    rely < 1 &&    !got_one ; i++)
                                            {
                                                if( firsttimey[i] )
                                                {
                                                    yi[i] = (g[n-i]|bit[n-i]);
                                                    firsttimey[i] = 0;
                                                }
                                                if(  !( (bit[i1]&yi[i]) && (bit[j1]&yi[i]) )    &&   !( (bit[i2]&yi[i]) && (bit[j2]&yi[i]) )   )
                                                {
                                                    rely++;  							                                  								
                                                    if(  mult > multar1[n-i]  &&  multar1[n-i] != 0 )	
                                                        got_one = 1;
                                                }
                                            }
                                        }


                                        if( !got_one )
                                        {
                                            epair1[numpair1].first = e1;
                                            epair1[numpair1].sec = e2;
                                            epair1[numpair1].multp = mult;
                                            epairorbit1[numpair1] = numpair1;
                                            antipair[e1][e2] = numpair1;
                                            numpair1++;
                                        }

                                    }						

                            }	
                        }													                    								

                        pepair = epair1;
                        pedge = edge1;

                        pnumpair = &numpair1;
                        pantipair = antipair;
                        pantiedge = antiedge;
                        pepairorbit = epairorbit1;


                    if( colours == accept )
                    {
                        options.userautomproc = userautom1;			
                        options.writeautoms = FALSE;
                        options.writemarkers = FALSE;
                        options.getcanon = FALSE;
                        options.defaultptn = FALSE;
                        nauty(g, lab, ptn, NILSET, orbits, &options, &stats, workspace, 50, m, n+1, NILGRAPH); 	
                        if( connec1 )
                        {
                            conf = TRUE;
                            if( !connectflag )
                                conf = isconnected(g, n+1);
                        }
                        extend(n+1, g, edge1, epair1, numpair1, epairorbit1, multar1, zar1, col00w1, conf);
                    }
                    else
                    {												
                        for (i = 0; i < numdovi; i++)
                         	doviorbit[i] = i ;

                        pdoviorbit = doviorbit;
                        pdovi = dovi;
                        pantidovi = antidovi;
                        pCNT = &numdovi;


                        options.userautomproc = userautom3;
                        options.writeautoms = FALSE;
                        options.writemarkers = FALSE;
                        options.getcanon = TRUE;
                        options.defaultptn = FALSE;
                        nauty( g, lab, ptn, 0, orbits, &options, &stats, workspace, 50, m, n+1, gi);  	
                        for (i = 0; i < n+1; i++)
                            l[lab[i]]=i;


                        dovimax.base = -1;
                        dovimax.first1= -1;
                        dovimax.sec1 = -1;
                        dovimax.first2 = -1;
                        dovimax.sec2 = -1;

                        for (i = 0; i < numdovi; i++)
                        {
                            b   = l[dovi[i].base];
                            if( b >= dovimax.base )
                            {
                                vt1 = l[dovi[i].first1];
                                vt2 = l[dovi[i].sec1];
                                vt3 = l[dovi[i].first2];
                                vt4 = l[dovi[i].sec2];


                                if( vt1 > vt2 )
                                {
                                    temp = vt1;
                                    vt1 = vt2;
                                    vt2 = temp;
                                }
                                if( vt3 > vt4 )
                                {
                                    temp = vt3;
                                    vt3 = vt4;
                                    vt4 = temp;
                                }
                                if( vt1 > vt3 )
                                {
                                    temp = vt1;
                                    vt1 = vt3;
                                    vt3 = temp;

                                    temp = vt2;
                                    vt2 = vt4;
                                    vt4 = temp;
                                }

                                if( (b > dovimax.base)  ||
                                    (  vt1 > dovimax.first1 )  ||
                                    ( vt1 == dovimax.first1  &&  vt2 > dovimax.sec1)  ||
                                    ( vt1 == dovimax.first1  &&  vt2 == dovimax.sec1  &&   vt3 > dovimax.first2) ||
                                    ( vt1 == dovimax.first1  &&  vt2 == dovimax.sec1  &&   vt3 == dovimax.first2  &&  vt4 > dovimax.sec2)    )
                                {
                                    dovimax.base = b;
                                    dovimax.first1 = vt1;
                                    dovimax.sec1 = vt2;
                                    dovimax.first2 = vt3;
                                    dovimax.sec2 = vt4;

                                    maxdovi = i;
                                }
                            }

                        }

                        if (doviorbit[maxdovi] == doviorbit[numdovi-1])
                        {
                            if( connec1 )
                            {
                                conf = TRUE;
                                if( !connectflag )
                                    conf = isconnected(g, n+1);
                            }
                            extend(n+1, g, edge1, epair1, numpair1, epairorbit1, multar1, zar1, col00w1, conf);
                        }
                    }



                }  // end if col..= rej..

            }  // end else >>> (if namx ==n+1)

            g[vm1] ^= bit[n] | bit[vm2];
            g[vm2] ^= bit[n] | bit[vm1];
            g[vm3] ^= bit[n] | bit[vm4];
            g[vm4] ^= bit[n] | bit[vm3];
            g[n] = 0;

        }  // end if c..

    } //end for c..

}

/*****************************************************************************
 *                                                                            *
 *  userautom1(count,perm,orbits,numorbits,stabvertex,n) is a simple          *
 *  version of the procedure named by options.userautomproc.                  *
 *                                                                            *
 *****************************************************************************/

static void
userautom1(int count, int *perm, int *orbits, int numorbits, int stabvertex, int n)
{
    int  epairperm[MAXP];
    int  vn1, vn2, vn3, vn4, etmp, i, e1, e2;

    epairperm[0] = 0;
    for (i = 1; i < *pnumpair; i++)
    {
        e1 = pepair[i].first;
        e2 = pepair[i].sec;

        vn1 = perm[pedge[e1].first];
        vn2 = perm[pedge[e1].sec];
        vn3 = perm[pedge[e2].first];
        vn4 = perm[pedge[e2].sec];

        if( vn1 > vn2 )
        {
            etmp = vn1;
            vn1 = vn2;
            vn2 = etmp;
        }
        if( vn3 > vn4 )
        {
            etmp = vn3;
            vn3 = vn4;
            vn4 = etmp;
        }
        if( vn1 > vn3 )
        {
            etmp = vn1;
            vn1 = vn3;
            vn3 = etmp;

            etmp = vn2;
            vn2 = vn4;
            vn4 = etmp;
        }

        e1 = pantiedge[vn1][vn2];
        e2 = pantiedge[vn3][vn4];
        epairperm[i] = pantipair[e1][e2];  			
    }
    orbjoin(pepairorbit,epairperm,*pnumpair);
}

/*****************************************************************************
 *                                                                            *
 *  userautom2(count,perm,orbits,numorbits,stabvertex,n) is a simple          *
 *  version of the procedure named by options.userautomproc.                  *
 *                                                                            *
 *****************************************************************************/
static void
userautom2(int count, int *perm, int *orbits, int numorbits, int stabvertex, int n)
{
    int 	doviperm[3*MAXN];
    int             vn1, vn2, vn3, vn4, vnb, i, etmp;

    for (i = 0; i < *pCNT; i++)
    {
        vnb = perm[pdovi[i].base];
        vn1 = perm[pdovi[i].first1];
        vn2 = perm[pdovi[i].sec1];
        vn3 = perm[pdovi[i].first2];
        vn4 = perm[pdovi[i].sec2];


        if( vn1 > vn2 )
        {
            etmp = vn1;
            vn1 = vn2;
            vn2 = etmp;
        }
        if( vn3 > vn4 )
        {
            etmp = vn3;
            vn3 = vn4;
            vn4 = etmp;
        }
        if( vn1 > vn3 )
        {
            etmp = vn1;
            vn1 = vn3;
            vn3 = etmp;

            etmp = vn2;
            vn2 = vn4;
            vn4 = etmp;
        }			

        doviperm[i] = pantidovi[vnb][vn2]; 	
    } 						 		
    orbjoin(pdoviorbit,doviperm,*pCNT); 										
}

/*****************************************************************************
 *                                                                            *
 *  userautom3(count,perm,orbits,numorbits,stabvertex,n) is a simple          *
 *  version of the procedure named by options.userautomproc.                  *
 *                                                                            *
 *****************************************************************************/
static void
userautom3( int count, int *perm, int *orbits, int numorbits, int stabvertex, int n)
{

    int 	doviperm[3*MAXN], epairperm[MAXP];
    int             vn1, vn2, vn3, vn4, vnb, i, etmp, e1, e2;

    for (i = 0; i < *pCNT; i++)
    {
        vnb = perm[pdovi[i].base];
        vn1 = perm[pdovi[i].first1];
        vn2 = perm[pdovi[i].sec1];
        vn3 = perm[pdovi[i].first2];
        vn4 = perm[pdovi[i].sec2];

        if( vn1 > vn2 )
        {
            etmp = vn1;
            vn1 = vn2;
            vn2 = etmp;
        }
        if( vn3 > vn4 )
        {
            etmp = vn3;
            vn3 = vn4;
            vn4 = etmp;
        }
        if( vn1 > vn3 )
        {
            etmp = vn1;
            vn1 = vn3;
            vn3 = etmp;

            etmp = vn2;
            vn2 = vn4;
            vn4 = etmp;
        }

        doviperm[i] = pantidovi[vnb][vn2];  	

    }					 		
    orbjoin(pdoviorbit,doviperm,*pCNT);		

    epairperm[0] = 0;
    for (i = 1; i < *pnumpair; i++)
    {
        e1 = pepair[i].first;
        e2 = pepair[i].sec;

        vn1 = perm[pedge[e1].first];
        vn2 = perm[pedge[e1].sec];
        vn3 = perm[pedge[e2].first];
        vn4 = perm[pedge[e2].sec];


        if( vn1 > vn2 )
        {
            etmp = vn1;
            vn1 = vn2;
            vn2 = etmp;
        }
        if( vn3 > vn4 )
        {
            etmp = vn3;
            vn3 = vn4;
            vn4 = etmp;
        }
        if( vn1 > vn3 )
        {
            etmp = vn1;
            vn1 = vn3;
            vn3 = etmp;

            etmp = vn2;
            vn2 = vn4;
            vn4 = etmp;
        }

        e1 = pantiedge[vn1][vn2];
        e2 = pantiedge[vn3][vn4];
        epairperm [i] = pantipair[e1][e2];
    }
    orbjoin(pepairorbit,epairperm,*pnumpair);   		
}

/*****************************************************************************
 *                                                                            *
 *  init_refinex(clr, lb, p, active, n) initializes some parameters           *
 *  of the the procedure refinex                                              *
 *                                                                            *
 *****************************************************************************/

static int
init_refinex( int *clr, int *lb, int *p, set *active, int n)
{
    register int	i, j, ci, ncell;
	
	ncell = 1;
	*active = bit[0];
    for (i = 0; i < n; i++)
	{
        ci = clr[i];
        for (j = i-1; (j >= 0)  &&    (clr[lb[j]] > ci) ; j--)
            lb[j+1] = lb[j];
        lb[j+1] = i;
    }
	lb[n] = n;
    for (i = 0; i < n; i++)
	{
        if( clr[lb[i]] != clr[lb[i+1]] )
		{
			p[i] = 0;
			ncell++;
			*active |= bit[i+1];
		}
		else
			p[i] = 1;

    }
	p[n] = 0;
	return ncell;
}

/*****************************************************************************
 *                                                                            *
 *  refinex(g,lab,ptn,level,numcells,count,active,goodret,code,m,n) is a      *
 *  custom version of refine() which can exit quickly if required.            *
 *                                                                            *
 *  Only use at level==0.                                                     *
 *  goodret : whether to do an early return for code 1                        *
 *  code := -1 for n-1 not max, 0 for maybe, 1 for definite                   *
 *                                                                            *
 *****************************************************************************/
static void
refinex(graph *g, int *lab, int *ptn, int level, int *numcells, int *count,
           set *active, boolean goodret, int *code, int m, int n)
{
    int     i, c1, c2, labc1, split1, split2, cell1, cell2, cnt, bmin, bmax;
    int     workperm[MAXN], bucket[MAXN+2];
    setword x, lact, workset;
    set     *gptr;

    if (n == 1)
    {
        *code = 1;
        return;
    }

    *code = 0;
    lact = *active;

    split1 = -1;
    while (*numcells < n && lact)
    {
        TAKEBIT(split1,lact);

        for (split2 = split1; ptn[split2] > 0; ++split2) {}
        if (split1 == split2)       /* trivial splitting cell */
        {
            gptr = GRAPHROW(g,lab[split1],1);
            for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
            {
                for (cell2 = cell1; ptn[cell2] > 0; ++cell2) {}
                if (cell1 == cell2) continue;

                c1 = cell1;
                c2 = cell2;
                while (c1 <= c2)
                {
                    labc1 = lab[c1];
                    if (ISELEMENT1(gptr,labc1))
                        ++c1;
                    else
                    {
                        lab[c1] = lab[c2];
                        lab[c2] = labc1;
                        --c2;
                    }
                }
                if (c2 >= cell1 && c1 <= cell2)
                {
                    ptn[c2] = 0;
                    ++*numcells;
                    lact |= bit[c1];
                }
            }
        }

        else        /* nontrivial splitting cell */
        {
            workset = 0;
            for (i = split1; i <= split2; ++i) workset |= bit[lab[i]];

            for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
            {
                for (cell2 = cell1; ptn[cell2] > 0; ++cell2) {}
                if (cell1 == cell2) continue;
                i = cell1;
                if ((x = workset & g[lab[i]]) != 0) cnt = POPCOUNT(x);
                else                                cnt = 0;
                count[i] = bmin = bmax = cnt;
                bucket[cnt] = 1;
                while (++i <= cell2)
                {
                    if ((x = workset & g[lab[i]]) != 0)
                        cnt = POPCOUNT(x);
                    else
                        cnt = 0;

                    while (bmin > cnt) bucket[--bmin] = 0;
                    while (bmax < cnt) bucket[++bmax] = 0;
                    ++bucket[cnt];
                    count[i] = cnt;
                }
                if (bmin == bmax) continue;
                c1 = cell1;
                for (i = bmin; i <= bmax; ++i)
                    if (bucket[i])
                    {
                        c2 = c1 + bucket[i];
                        bucket[i] = c1;
                        if (c1 != cell1)
                        {
                            lact |= bit[c1];
                            ++*numcells;
                        }
                        if (c2 <= cell2) ptn[c2-1] = 0;
                        c1 = c2;
                    }
                for (i = cell1; i <= cell2; ++i)
                    workperm[bucket[count[i]]++] = lab[i];
                for (i = cell1; i <= cell2; ++i) lab[i] = workperm[i];
            }
        }

        if (ptn[n-2] == 0)
        {
            if (lab[n-1] == n-1)
            {
                *code = 1;
                if (goodret) return;
            }
            else
            {
                *code = -1;
                return;
            }
        }
        else
        {
            i = n - 1;
            while (TRUE)
            {
                if (lab[i] == n-1) break;
                --i;
                if (ptn[i] == 0)
                {
                    *code = -1;
                    return;
                }
            }
        }
    }
}


/*****************************************************************************/

int main(int argc, char *argv[])
{
    //argc = #entered arguments by the user + 1;
    int 	n, cntr, numpair, numedge, i, j, i1, i2, j1, j2, multm, e1, e2, argnum, sw; //char  sw??
    int         multar[MAXN], col00w[MAXN], epairorbit[MAXP], lab[MAXN], ptn[MAXN], orbits[MAXN];
    int 	antipair[MAXE][MAXE], antiedge[MAXN][MAXN];
    setword 	x, y;
    pairstruct 	epair[MAXP];
    edgestruct  edge[MAXE];	
    graph 	g[MAXN];
    double      timebefore, timeafter;
    setword     zar[MAXN];
    boolean     badargs, gotf, gotmr, quiet;
    char        *arg;
    char        msg[201];

    HELP; PUTVERSION;
    nauty_check(WORDSIZE,1,MAXN,NAUTYVERSIONID);   ////?????

    if( MAXN > WORDSIZE )
    {
        fprintf(stderr,"quarticgen: incompatible MAXN or WORDSIZE\n");
        exit(1);
    }

    badargs = FALSE;
    connec1 = connec2 = FALSE;
    nooutput = FALSE;
    canonise = FALSE;  //?????
    graph6 = FALSE;
    sparse6 = FALSE;
    header = FALSE;
    quiet = FALSE;

    gotf = FALSE;
    gotmr = FALSE;
    outfilename = NULL;

    argnum = 0;

    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                SWBOOLEAN('u', nooutput)
                else SWBOOLEAN('l', canonise)
                    else SWBOOLEAN('c',connec1)
                        else SWBOOLEAN('C',connec2)
                            else SWBOOLEAN('g', graph6)
                                else SWBOOLEAN('s', sparse6)
                                    else SWBOOLEAN('h', header)
                                        else SWBOOLEAN('q', quiet)
#ifdef PLUGIN_SWITCHES
                                PLUGIN_SWITCHES
#endif
                                else badargs = TRUE;
            }
        }
        else if (arg[0] == '-' && arg[1] == '\0')
            gotf = TRUE;
        else
        {
            if (argnum == 0)
            {
                if (sscanf(arg,"%d",&nmax) != 1)
                    badargs = TRUE;
                ++argnum;
            }
            else if (gotf)
                badargs = TRUE;
            else
            {
                if (!gotmr)
                {
                    if (sscanf(arg,"%d/%d",&res,&mod) == 2)
                    {
                        gotmr = TRUE;
                        continue;
                    }
                }
                if (!gotf)
                {
                    outfilename = arg;
                    gotf = TRUE;
                    continue;
                }
            }
        }
    }



    if( !badargs  &&  gotmr  &&  (res < 0 || res >= mod) )
    {
        fprintf(stderr,
                "E quarticgen: must have 0 <= res < mod\n");
        badargs = TRUE;
    }

    if      (connec2) connec = 2;
    else if (connec1) connec = 1;
    else              connec = 0;

    if (!argnum )
        badargs = TRUE;
    else if (nmax < 1 || nmax > MAXN )
    {
        fprintf(stderr,
                ">E quarticgen: must have n =1..%d \n",MAXN);
        badargs = TRUE;
    }

    if (!gotmr)
    {
        mod = 1;
        res = 0;
    }
   /* else if (argnum == 5 || argnum > 6)
        badargs = TRUE;*/ // argnum will never exceeds 1 here andeven in genbg it never exceeds 2!


    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if ( (graph6!=0) + (sparse6!=0) + (nooutput!=0) > 1)
        gt_abort(">E quarticgen: -ugs are incompatible\n");

    if ( nooutput && header)
        gt_abort(">E quarticgen: -u -h are incompatible\n");


#ifdef OUTPROC
    outproc = OUTPROC;
#else
    if (nooutput) outproc = nullwrite;
    else if (sparse6)  outproc = writes6x;
    else               outproc = writeg6x;
#endif


#ifdef PLUGIN_INIT
    PLUGIN_INIT
#endif

    if (nooutput)
        outfile = stdout;
    else if (!gotf || outfilename == NULL)
    {
        outfilename = "stdout";
        outfile = stdout;
    }
    else if ((outfile = fopen(outfilename, "w")) == NULL)
    {
        fprintf(stderr,
                ">E quarticgen: can't open %s for writing\n",outfilename);
        gt_abort(NULL);
    }

    if (!quiet)
    {
        msg[0] = '\0';
        if (strlen(argv[0]) > 75) fprintf(stderr,">A %s",argv[0]);
        else CATMSG1(">A %s",argv[0]);

        CATMSG1(" n = %d", nmax);
   // if (connec) CATMSG0(connec2 ? "C" : connec1 ? "c" : "",);
        if (connec2) CATMSG0(" C"); else if (connec1) CATMSG0(" c");
        if (mod > 1) CATMSG2(" class= %d/%d", res, mod);
        CATMSG0("\n");
        fputs(msg,stderr);
        fflush(stderr);
    }

    if (header)
    {
	if (SPARSE6)
            writeline(outfile,SPARSE6_HEADER);        // No enter after the header is ok?
	else
            writeline(outfile,GRAPH6_HEADER);
        fflush(outfile);
    }

 //   if (mod > 1 && nmax > 10)
    if( mod > 1 )
    {
        if( nmax <= threshold )
            splitlevel = nmax - LEV1;
        else
            splitlevel = LEV2;
        splitcount = res;
    }
    else
    {
        splitlevel = -1;
        mod = 1;
        res = 0;  // narjess for the sake of n>splitlevel (below)
    }

    timebefore = CPUTIME;
    m = 1;

    cntr = numread = numwritten = 0;
    for(cntr = 0; cntr < NUMIRRED; cntr++)
    {
        n = graphsize(irred[cntr]);
        //////////////////////////////////////////////////////////////////
        if( n == splitlevel  )
        {
            if (splitcount-- != 0) continue;
            if( n == nmax )
                splitcount = mod-1;
            else
                splitcount = 0;
        }
        if( n > splitlevel )
        {
            if ( gotmr && res ) continue;
        }
        /////////////////////////////////////////////////////////////////////////


        if( n < nmax )
        {
            numread++;
            stringtograph(irred[cntr],g,m);
            epairorbit[0] = 0;
            numpair = 1;
            numedge = 0;
            for (i = 0; i < n; i++)
            {
                multar[i] = 10;

                x = g[i] & BITMASK(i);
                while (x)
                {
                    TAKEBIT(j,x);
                    edge[numedge].first = i;
                    edge[numedge].sec = j;
                    antiedge[i][j] = numedge;
                    numedge++;

                }
            }
            for (e1 = 0; e1 < numedge; e1++)
            {
                i1 = edge[e1].first;
                j1 = edge[e1].sec;
                for (e2 = e1+1; e2 < numedge; e2++)
                {
                    i2 = edge[e2].first;
                    j2 = edge[e2].sec;
                    if( i1 != i2   &&   j1 != i2   &&  j1 != j2   )
                    {
                        epair[numpair].first = e1;
                        epair[numpair].sec = e2;

                        multm = 1;

                        if(   (!(g[i1] & bit[i2]))  &&   (!(g[j1] & bit[j2]))  )
                            multm++;

                        if(   (!(g[i1] & bit[j2]))  &&   (!(g[j1] & bit[i2]))  )
                            multm++;
                        epair[numpair].multp = multm;

                        epairorbit[numpair] = numpair;
                        antipair[e1][e2] = numpair;
                        numpair++;
                    }
                }
            }
            pepair = epair;
            pedge = edge;
            pantipair = antipair;
            pantiedge = antiedge;
            pnumpair = &numpair;
            pepairorbit = epairorbit;

            options.userautomproc = userautom1;
            options.writeautoms = FALSE;
            options.writemarkers = FALSE;
            options.getcanon = FALSE;
            options.defaultptn = TRUE;

            nauty(g, lab, ptn, NILSET, orbits, &options, &stats, workspace, 50, m, n, NILGRAPH);
            extend(n, g, edge, epair, numpair, epairorbit, multar, zar, col00w, isconnected(g, n));
        }
        if( n == nmax )
        {
            numread++;
            stringtograph(irred[cntr],g,m);

            if( connec2 && !isbiconnected(g, n) )
                continue;
            if( connec1 && !isconnected(g, n) )
                continue;
            (*outproc)( outfile, g, n);
            numwritten++;
        }
    }
    timeafter=CPUTIME;

    if (!quiet)
    {
	if (nooutput)
	    fprintf(stderr,">Z " COUNTER_FMT " graphs generated in %3.2f seconds\n",
		numwritten,timeafter-timebefore);
	else
	    fprintf(stderr,">Z " COUNTER_FMT " graphs written to %s in %3.2f seconds\n",
		numwritten,outfilename,timeafter-timebefore);
    }

    exit(0);
}
