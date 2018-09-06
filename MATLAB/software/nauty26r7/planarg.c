/* planarg : test for planarity and find embedding or obstruction */

#define USAGE "planarg [-v] [-nVq] [-p|-u] [infile [outfile]]"

#define HELPTEXT \
" For each input, write to output if planar.\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -v  Write non-planar graphs instead of planar graphs\n\
    -V  Write report on every input\n\
    -u  Don't write anything, just count\n\
    -p  Write in planar_code if planar (without -p, same format as input)\n\
    -k  Follow each non-planar output with an obstruction in sparse6\n\
        format (implies -v, incompatible with -p)\n\
    -n  Suppress checking of the result\n\
    -q  Suppress auxiliary information\n\
\n\
    This program permits multiple edges and loops\n"

/*************************************************************************/

#include "gtools.h" 
#include "planarity.h"

/*************************************************************************/

static void
write_planarcode(FILE *f, t_ver_sparse_rep *VR, t_adjl_sparse_rep *A,
                 t_embed_sparse_rep *ER, int n, int ne)
/* Write the embedding to f in planar_code */
{
    int bytes;
    size_t i,j,len,k,k0;
    unsigned int w;
    DYNALLSTAT(unsigned char,buff,buff_sz);
#define PUT1(x) buff[j++]=(x);
#define PUT2(x) w=(x); buff[j++]=(w>>8)&0xFF; buff[j++]=w&0xff;
#define PUT4(x) w=(x); buff[j++]=(w>>24)&0xFF; buff[j++]=(w>>16)&0xff; \
             buff[j++]=(w>>8)&0xFF; buff[j++]=w&0xff;

    if (n <= 255)        bytes = 1;
    else if (n <= 65535) bytes = 2;
    else                 bytes = 4;

    len = bytes * (1 + n + 2*(size_t)ne);
    if (bytes == 2)      len += 1;
    else if (bytes == 4) len += 3;

    DYNALLOC1(unsigned char,buff,buff_sz,len,"planarg");

    if (bytes == 1)
    {
        j = 0;
        PUT1(n);
        for (i = 0; i < n; ++i)
        {
            k = k0 = VR[i].first_edge;
            if (k != NIL)
                do
                {
                    PUT1(A[ER[k].in_adjl].end_vertex+1);
                    k = ER[k].next;
                } while (k != k0);
            PUT1(0);
        }
    }
    else if (bytes == 2)
    {
        j = 0;
        PUT1(0);
        PUT2(n);
        for (i = 0; i < n; ++i)
        {
            k = k0 = VR[i].first_edge;
            if (k != NIL)
                do
                {
                    PUT2(A[ER[k].in_adjl].end_vertex+1);
                    k = ER[k].next;
                } while (k != k0);
            PUT2(0);
        }
    }
    else
    {
        j = 0;
        PUT1(0);
        PUT2(0);
        PUT4(n);
        for (i = 0; i < n; ++i)
        {
            k = k0 = VR[i].first_edge;
            if (k != NIL)
                do
                {
                    PUT4(A[ER[k].in_adjl].end_vertex+1);
                    k = ER[k].next;
                } while (k != k0);
            PUT4(0);
        }
    }

    if (fwrite((void*)buff,1,len,f) != len)
    {
        fprintf(stderr,">E write_planarcode : error on writing\n");
        ABORT(">E write_planarcode");
    }
}

/*************************************************************************/

static boolean
isplanar(t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A, int e,
         int *c, t_ver_sparse_rep **VR, t_adjl_sparse_rep **AR,
         t_embed_sparse_rep  **ER, int *nbr_e_obs,
         boolean planarcheck, boolean nonplanarcheck)
/*
  The input graph is given as an adjacency list:
      V: array of vertices 
      n: size of graph
      A: adjacency list
      e: number of edges
      
  If the graph is planar the embedding is stored  in VR and ER;
  the embedding contains e edges (nbr_e_obs not used)
      
  If the graph is non planar the obstruction is returned in
  VR and AR together with the number of edges in nbr_e_obs.

  In all cases the number of components is return in c.

  planarcheck and nonplanarcheck determine if checking is done.
  The embedding and obstruction outputs are only made if the
  appropriate checking is done.
*/
{
    t_dlcl           **dfs_tree, **back_edges, **mult_edges;
    int              edge_pos, v, w;
    boolean          ans;
    t_ver_edge       *embed_graph;
 
    ans = sparseg_adjl_is_planar(V, n, A, c,
                                 &dfs_tree, &back_edges, &mult_edges,
                                 &embed_graph, &edge_pos, &v, &w);
    
    if (!ans && nonplanarcheck)
    {
        embedg_obstruction(V, A, dfs_tree, back_edges,
                           embed_graph, n, &edge_pos,
                           v, w, VR, AR, nbr_e_obs);
    }
    else if (planarcheck)
    {
        embedg_embedding(V, A, embed_graph, n, e, *c, edge_pos, mult_edges,
                         VR, ER);
    }
    
    sparseg_dlcl_delete(dfs_tree, n);
    sparseg_dlcl_delete(back_edges, n);
    sparseg_dlcl_delete(mult_edges, n);
    embedg_VES_delete(embed_graph, n);
 
    return ans;
}

/*************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    sparsegraph sg;
    boolean badargs;
    boolean verbose,nonplanar,quiet;
    boolean planarcode,nowrite,nocheck;
    int i,j,k,n,argnum,ne,loops;
    int codetype,outcode;
    t_ver_sparse_rep *VR;
    t_adjl_sparse_rep *AR;
    t_embed_sparse_rep *ER;
    int nbr_c,nbr_e_obs;
    nauty_counter nin,nout,nplan;
    char *arg,sw;
    double t0,tp,tnp,netotalp,netotalnp;
    DYNALLSTAT(t_ver_sparse_rep,V,V_sz);
    DYNALLSTAT(t_adjl_sparse_rep,A,A_sz);

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    quiet = nowrite = planarcode = FALSE;
    verbose = nonplanar = nocheck = FALSE;

    argnum = 0;
    badargs = FALSE;
    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                     SWBOOLEAN('v',nonplanar)
                else SWBOOLEAN('q',quiet)
                else SWBOOLEAN('V',verbose)
                else SWBOOLEAN('p',planarcode)
                else SWBOOLEAN('u',nowrite)
                else SWBOOLEAN('n',nocheck)
                else badargs = TRUE;
            }
        }
        else
        {
            ++argnum;
            if      (argnum == 1) infilename = arg;
            else if (argnum == 2) outfilename = arg;
            else                  badargs = TRUE;
        }
    }

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        exit(1);
    }

    if (planarcode && nonplanar)
    {
        fprintf(stderr,">E planarg: -p and -v are incompatible\n");
        exit(1);
    }

    if (!quiet)
    {
        fprintf(stderr,">A planarg");
        if (nonplanar||planarcode||nowrite||nocheck)
            fprintf(stderr," -");
        if (nonplanar) fprintf(stderr,"v");
        if (nowrite) fprintf(stderr,"u");
        if (planarcode) fprintf(stderr,"p");
        if (nocheck) fprintf(stderr,"n");
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,FALSE,1);
    if (!infile) exit(1);
    if (!infilename) infilename = "stdin";

    NODIGRAPHSYET(codetype);

    if (!nowrite)
    {
        if (!outfilename || outfilename[0] == '-')
        {
            outfilename = "stdout";
            outfile = stdout;
        }
        else if ((outfile = fopen(outfilename,"w")) == NULL)
        {
            fprintf(stderr,"Can't open output file %s\n",outfilename);
            gt_abort(NULL);
        }

        if (planarcode)            outcode = PLANARCODE;
        else if (codetype&SPARSE6) outcode = SPARSE6;
        else                       outcode = GRAPH6;

        if (outcode == PLANARCODE)
            writeline(outfile,PLANARCODE_HEADER);
        else if (codetype&HAS_HEADER)
        {
            if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
            else                    writeline(outfile,GRAPH6_HEADER);
        }
    }

    nin = nout = nplan = 0;
    netotalp = netotalnp = 0.0;
    SG_INIT(sg);

    tp = tnp = 0.0;
    while (TRUE)
    {
        if (read_sg_loops(infile,&sg,&loops) == NULL) break;
        n = sg.nv;
        ne = (sg.nde+loops)/2;
        ++nin;

        DYNALLOC1(t_ver_sparse_rep,V,V_sz,n,"planarg");
        DYNALLOC1(t_adjl_sparse_rep,A,A_sz,2*ne+1,"planarg");
        k = 0;
        for (i = 0; i < n; ++i)
            if (sg.d[i] == 0)
                V[i].first_edge = NIL;
            else
            {
                V[i].first_edge = k;
                for (j = sg.v[i]; j < sg.v[i]+sg.d[i]; ++j)
                {
                    A[k].end_vertex = sg.e[j];
                    A[k].next = k+1;
                    ++k;
                    if (A[k-1].end_vertex == i)  /* loops go in twice */
                    {
                        A[k].end_vertex = i;
                        A[k].next = k+1;
                        k++;
                    }
                }
                A[k-1].next = NIL;
            }

        if (k != 2*ne)
        {
            fprintf(stderr,
               ">E planarg: decoding error; nin=" COUNTER_FMT "\n",nin);
            fprintf(stderr,"n=%d nde=%lu ne=%d loops=%d\n",
                    n,(unsigned long)sg.nde,ne,loops);
            exit(1);
        }

        VR = NULL;
        AR = NULL;
        ER = NULL;
        t0 = CPUTIME;
        if (isplanar(V,n,A,ne,&nbr_c,&VR,&AR,&ER,&nbr_e_obs,
                     !nocheck||planarcode,!nocheck))
        {
            ++nplan;
            tp += CPUTIME - t0;
            netotalp += ne;
            if (!nowrite && !nonplanar)
            {
                if (planarcode)
                    write_planarcode(outfile,VR,A,ER,n,ne);
                else
                    writelast(outfile);
                ++nout;
            }
            if (verbose)
                fprintf(stderr,"graph " COUNTER_FMT ": n=%d ne=%d planar\n",
                        nin,n,ne);
        }
        else
        {
            tnp += CPUTIME - t0;
            netotalnp += ne;
            if (!nowrite && nonplanar)
            {
                writelast(outfile);
                ++nout;
            }
            if (verbose)
                fprintf(stderr,"graph " COUNTER_FMT ": n=%d ne=%d non-planar\n",
                        nin,n,ne);
        }
        FREES(VR);
        FREES(AR);
        FREES(ER);
    }


    if (!nowrite)
    {
        if (!quiet)
            fprintf(stderr,
            ">Z  " COUNTER_FMT " graphs read from %s, "
                  COUNTER_FMT " written to %s; %3.2f sec.\n",
                nin,infilename,nout,outfilename,tp+tnp);
    }
    else
    {
        fprintf(stderr,
                " " COUNTER_FMT " graphs input\n "
                    COUNTER_FMT " graphs planar\n",nin,nplan);
        fprintf(stderr," cpu = %3.2f sec. ",tp+tnp);
        if (netotalp)
            fprintf(stderr," planar:%.3f",1e5*tp/netotalp);
        if (netotalnp)
            fprintf(stderr," nonplanar:%.3f",1e5*tnp/netotalnp);
        fprintf(stderr,"\n");
    }

    return 0;
}
