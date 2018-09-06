/* planarity.c - code for planarity testing of undirected graphs.
 * Method of Boyer and Myrvold, programmed by Paulette Lieby.
 * The copyright of this program is owned by the Magma project.
 * Distributed with nauty by permission.
 ***************************************************************/
 
/*
 *  sparseg_adjl.c
 */
 
/*
  What:
  *****
  
  Implementing:

  Some high-level functions on the sparse graph as
  an adjacency list.
  In particular, testing if it is planar.
  


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

#include "planarity.h"

#define IF_DEB(x)    {}
#define IF_VERB(x)   {}



/* aproto: header embed_graph_protos.h */


#ifndef PLANAR_IN_MAGMA
#endif


boolean 
sparseg_adjl_plan_and_iso (t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A,
	int e, int *c, t_ver_sparse_rep **VR, t_adjl_sparse_rep **AR,
	t_embed_sparse_rep **ER, int *nbr_e_obs)
    /*
      the input graph is given as an adjacency list:
      V: array of vertices 
      n: size of graph
      A: adjacency list
      e: number of edges
      
      if the graph is planar the embedding is stored  in VR and ER;
      the embedding contains e edges
      (nbr_e_obs not used)
      
      if the graph is non planar the obstruction is returned in
      VR and AR together with the number of edges in nbr_e_obs

      in all cases is also returned the number of components (in c)
    */
{
    t_dlcl           **dfs_tree, **back_edges, **mult_edges;
    int              edge_pos, v, w;
    boolean          ans;
    t_ver_edge       *embed_graph;
 
    ans = sparseg_adjl_is_planar(V, n, A, c,
                                 &dfs_tree, &back_edges, &mult_edges,
                                 &embed_graph, &edge_pos, &v, &w);
    
    if (!ans)
    {
        embedg_obstruction(V, A, dfs_tree, back_edges,
                           embed_graph, n, &edge_pos,
                           v, w, VR, AR, nbr_e_obs);
    }
    else
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



int *
sparseg_adjl_footprint (t_ver_sparse_rep *V, int n,
	t_adjl_sparse_rep *A, int v)
    /*
      return v's footprint:
      an array fp of size n where fp[i] = index of (directed)
      edge [v, i] in A
    */
{
    /*
      note that we won't initialise the array:
      its subsequent usage doesn't require it
    */
    int        *fp, e;

    fp = (int *) mem_malloc(sizeof(int) * n);

    if (V[v].first_edge == NIL)
        /*
          do nothing
        */
        return fp;

    e = V[v].first_edge;
    while (e != NIL)
    {
        fp[A[e].end_vertex] = e;
        e = A[e].next;
    }

    return fp;
}


void 
sparseg_adjl_print (t_ver_sparse_rep *V, int n,
	t_adjl_sparse_rep *A, boolean user_level)
{
    int        v;

    for (v = 0; v < n; v++)
    {
        int     next;

        if (user_level)
            fprintf(stdout, "%d:\t", v + 1);
        else
            fprintf(stdout, "%d:\t", v);
        
        next = V[v].first_edge;
        while (next != NIL)
        {
            if (user_level)
                fprintf(stdout, "%d ", A[next].end_vertex + 1);
            else
                fprintf(stdout, "%d ", A[next].end_vertex);
            
            next = A[next].next;
        }
        fprintf(stdout, "\n");
    }
}




void 
sparseg_adjl_embed_print (t_ver_sparse_rep *V_e, int n,
	t_adjl_sparse_rep *A, t_embed_sparse_rep *E, boolean user_level)
    /*
      print the embedding given by E,
      edges are referred to by their index in A

      and V_e[v].first_edge is the index in E of the first edge
      (in the embedding's order) incident from v

      note that E is NOT indexed by the same vertices' array
      that indexes A (at the creation of the sparse graph)
    */
{
    int        v;

    for (v = 0; v < n; v++)
    {
        int      start, next;

        if (user_level)
            fprintf(stdout, "%d:\t", v + 1);
        else
            fprintf(stdout, "%d:\t", v);

        if (V_e[v].first_edge == NIL)
        {
            fprintf(stdout, "\n");
            continue;
        }
        start = next = V_e[v].first_edge;

        if (user_level)
            fprintf(stdout, "%d ", A[ E[next].in_adjl ].end_vertex + 1);
        else
            fprintf(stdout, "%d ", A[ E[next].in_adjl ].end_vertex);

        next = E[next].next;
        
        while (next != start)
            /*
              recall that in E edges are linked into a circular list
            */
        {
            if (user_level)
                fprintf(stdout, "%d ", A[ E[next].in_adjl ].end_vertex + 1);
            else
                fprintf(stdout, "%d ", A[ E[next].in_adjl ].end_vertex);

            next = E[next].next;
        }
        fprintf(stdout, "\n");
    }
}

graph *
sparseg_adjl_to_nauty_graph (t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A)
    /*
      write the sparse graph as a nauty graph
    */
{
    int          m, v, e, i;
    graph        *g;
 
    m = (n + WORDSIZE - 1) / WORDSIZE;
    g = (graph *) mem_malloc(n * m * sizeof(graph));
    for (i = (long) m * n; --i >= 0;)
        g[i] = 0;
 
    /*
      we first copy V and A's information into g
    */
    for (v = 0; v < n; v++)
    {
        e = V[v].first_edge;
        while (e != NIL)
            /*
              A[e].end_vertex is the next neighbour in the list,
              A[e].next points to the next edge in the list
            */
        {
            if (A[e].end_vertex != v)  /* no loops */
            {
                ADDELEMENT(GRAPHROW(g, v, m), A[e].end_vertex);
            }
            e = A[e].next;
        }
    }

    return g;
}



#if 0
t_edge_sparse_rep *
sparseg_adjl_edges (t_ver_sparse_rep *V, int n,
	t_adjl_sparse_rep *A, int e, boolean digraph)
    /*
      e is the number of edges
    */
{
    t_edge_sparse_rep *edges;
    int               m, u, v, pos_e;
    graph             *g;

    edges = (t_edge_sparse_rep *) mem_malloc(sizeof(t_edge_sparse_rep) * e);

    m = (n + WORDSIZE - 1) / WORDSIZE;
    g = sparseg_adjl_to_nauty_graph(V, n, A);

    pos_e = 0;
    for (u = 0; u < n; u++)
    {
        v = digraph == TRUE ? 0 : u + 1;
        for (; v < n; v++)
        {
            if (ISELEMENT(GRAPHROW(g, u, m), v))
            {
                t_edge_sparse_rep edge;

                edge.ends[0] = u;
                edge.ends[1] = v;
                edges[pos_e++] = edge;
            }
        }
    }
    ASSERT(pos_e == e);
    mem_free(g);
    
    return edges;
}
#endif



t_edge_sparse_rep *
sparseg_adjl_edges (t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A,
	int e, boolean digraph)
    /*
      e is the number of edges
    */
{
#if 0
    t_edge_sparse_rep   *edges;
    int                 u, v, pos_e, *loops, *foot_print;
    graph               *g;
 
    loops = (int *) mem_malloc(sizeof(int) * n);
    for (v = 0; v < n; v++)
    {
        loops[v] = 0;
    }
 
    edges = (t_edge_sparse_rep *) mem_malloc(sizeof(t_edge_sparse_rep) * e);
    pos_e = 0;

    foot_print = (int *) mem_malloc(sizeof(int) * n);
    for (u = 0; u < n; u++)
        foot_print[u] = NIL;
    
    for (v = 0; v < n; v++)
    {
        int               ne;
        t_edge_sparse_rep edge;

        ne = V[v].first_edge;
        while (ne != NIL)
        {
            u = A[ne].end_vertex;
            if (digraph
                || (!digraph && u > v))
            {
                foot_print[u] = v;
            }
            else if (!digraph && u == v)
            {
                if (loops[v] == 0)
                {
                    foot_print[u] = v;
                }
 
                loops[v] ^= 1;
            }

            ne = A[ne].next;
        }

        for (u = 0; u < n; u++)
            if (foot_print[u] == v)
            {
                edge.ends[0] = v;
                edge.ends[1] = u;
                edges[pos_e++] = edge;
            }
    }
    ASSERT(pos_e == e);
    mem_free(loops);
    mem_free(foot_print);
    
    return edges;
    
#endif   
    /*
      there must be a simpler way
    */
#if 0
    typedef struct edge_list {
        int               size;
        t_edge_sparse_rep *edges;
    } t_edge_list;
 
    t_edge_list         *edge_table;        
    t_edge_sparse_rep   *edges;
    int                 u, v, nbr_e, pos_e, *loops;
    graph               *g;

    loops = (int *) mem_malloc(sizeof(int) * n);
    for (v = 0; v < n; v++)
    {
        loops[v] = 0;
    }

    /*
      now create an edge table as follows:
      - there are n lists in total
      - their respective size is given by size
      - their contents by *edges:

      edge_table[i] will contain all the edges whose end-point is i:
      these edges, by construction, will be sorted according to their
      starting point

      what for? to finish off each start-vertex processing
      with a bucket sort so that
      the edges are sorted wrt start- & end-point

      bucket sort is linear, hence why...
    */
    edge_table = (t_edge_list *) mem_malloc(sizeof(t_edge_list) * n);
    for (v = 0; v < n; v++)
    {
        edge_table[v].size = 0;
        edge_table[v].edges = NP;
    }

    edges = (t_edge_sparse_rep *) mem_malloc(sizeof(t_edge_sparse_rep) * e);
    
    nbr_e = 0;
    pos_e = 0;
    for (v = 0; v < n; v++)
    {
        int    ne, w, u;

        ne = V[v].first_edge;
        while (ne != NIL)
        {
            u = A[ne].end_vertex;
            if (digraph
                || (!digraph && u > v))
            {
                t_edge_sparse_rep edge;
 
                edge.ends[0] = v;
                edge.ends[1] = u;

                /*
                  now stick this edge into the table: one may ponder
                  as to the cost of constantly reallocating memory...
                  some cursory tests in another context tell me that
                  this is pretty much ok
                  (and certainly better than allocating n^2 storage space)
                */
                if (edge_table[u].size == 0)
                {
                    edge_table[u].edges = (t_edge_sparse_rep *)
                        mem_malloc(sizeof(t_edge_sparse_rep));
                }
                else
                {
                    edge_table[u].edges = (t_edge_sparse_rep *)
                         mem_realloc(edge_table[u].edges,
                                     sizeof(t_edge_sparse_rep)
                                     * (edge_table[u].size + 1));
                }

                (edge_table[u].edges)[edge_table[u].size] = edge;
                edge_table[u].size += 1;
                nbr_e++;
            }
            else if (!digraph && u == v)
            {
                if (loops[v] == 0)
                {
                    t_edge_sparse_rep edge;
 
                    edge.ends[0] = v;
                    edge.ends[1] = u;

                    if (edge_table[u].size == 0)
                    {
                        edge_table[u].edges = (t_edge_sparse_rep *)
                            mem_malloc(sizeof(t_edge_sparse_rep));
                    }
                    else
                    {
                        edge_table[u].edges = (t_edge_sparse_rep *)
                            mem_realloc(edge_table[u].edges,
                                        sizeof(t_edge_sparse_rep)
                                        * (edge_table[u].size + 1));
                    }
                    
                    (edge_table[u].edges)[edge_table[u].size] = edge;
                    edge_table[u].size += 1;
                    nbr_e++;
                }
                
                loops[v] ^= 1;
            }
            
            ne = A[ne].next;
        }

        /*
          bucket sort must take place here:
          of course the whole lot is not exactly linear!
          since we perform the sort n times; but we can hope for
          a "good"  ?? average behaviour:

          in any case this must be better that checking adjacencies
          n^2 times in a sparse rep. (see edge_set_iset_assure)
        */
        for (w = 0; w < n; w++)
        {
            if (edge_table[w].size > 0)
            {
                for (u = 0; u < edge_table[w].size; u++)
                {
                    ASSERT((edge_table[w].edges)[u].ends[0] == v);
                    edges[pos_e++] = (edge_table[w].edges)[u];
                }
                mem_free(edge_table[w].edges);
                edge_table[w].size = 0;
                edge_table[w].edges = NP;
            }
        }
    }
    ASSERT(nbr_e == e);
    ASSERT(pos_e == e);
    mem_free(loops);
    mem_free(edge_table);
    
    return edges;
#endif

    t_edge_sparse_rep *edges;
    int               v, pos_e, *loops;

    edges = (t_edge_sparse_rep *) mem_malloc(sizeof(t_edge_sparse_rep) * e);
    loops = (int *) mem_malloc(sizeof(int) * n);
    for (v = 0; v < n; v++)
    {
        loops[v] = 0;
    }
    
    pos_e = 0;
    for (v = 0; v < n; v++)
    {
        int    ne;

        ne = V[v].first_edge;
        while (ne != NIL)
        {
            int      u;

            u = A[ne].end_vertex;
            if (digraph
                || (!digraph && u > v))
            {
                t_edge_sparse_rep edge;
 
                edge.ends[0] = v;
                edge.ends[1] = u;
                edges[pos_e++] = edge;
            }
            else if (!digraph && u == v)
            {
                if (loops[v] == 0)
                {
                    t_edge_sparse_rep edge;
 
                    edge.ends[0] = v;
                    edge.ends[1] = u;
                    edges[pos_e++] = edge;
                }

                loops[v] ^= 1;
            }
            ne = A[ne].next;
        }
    }
    ASSERT(pos_e == e);
    mem_free(loops);
    
    return edges;

}

/*
 *  sparseg_adjl_modify.c
 */
 
/*
  What:
  *****
  
  Implementing:
 
  Some high-level functions on the sparse graph as
  an adjacency list.
  In particular, adding/removing vertices/edges.


  NOTE: Most of the functions implicitely assume that the
  graph is undirected;
  this must be slightly rewritten for the general case
  -- just haven't got the time right now...


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

#include "planarity.h"

#define IF_DEB(x)    {}
#define IF_VERB(x)   {}


/* aproto: header embed_graph_protos.h */



#ifndef PLANAR_IN_MAGMA
#endif



boolean 
sparseg_adjl_add_edge (t_ver_sparse_rep *V, int n, t_adjl_sparse_rep **A,
	int *size_A, int *pos, int u, int v, boolean CHECK)
    /*
      add the UNDIRECTED edge to the sparse graph (V, n, A)
      - pos records where to add the next edge in A
      - if pos + 1 == size_A, we must extend A

      we check if the edge is already in the graph iff CHECK true

      also we assume that the graph (V, n, A) is undirected
    */
{
    boolean             edge_exists;

    edge_exists = FALSE;
    if (CHECK)
    {
        edge_exists = sparseg_adjl_dir_edge_exists(V, n, *A, u, v);

        if (edge_exists)
            return FALSE;
    }
    
    if (*pos == *size_A)
    {
        IF_DEB(
               fprintf(stdout, "realloc \n");
               )

        *size_A += 2;    /* add two directed edges */
        *A = (t_adjl_sparse_rep *)
            mem_realloc(*A, sizeof(t_adjl_sparse_rep) * *size_A);
    }
    else if (*pos + 1 == *size_A)
    {
        IF_DEB(
               fprintf(stdout, "realloc \n");
               )

        *size_A += 1;    /* add two directed edges */
        *A = (t_adjl_sparse_rep *)
            mem_realloc(*A, sizeof(t_adjl_sparse_rep) * *size_A);
    }
    ASSERT(*pos + 1 < *size_A);
    
    sparseg_adjl_add_dir_edge(V, n, A, size_A, pos, u, v, FALSE);
    sparseg_adjl_add_dir_edge(V, n, A, size_A, pos, v, u, FALSE);

    return TRUE;
}
 
boolean 
sparseg_adjl_add_edge_no_extend (t_ver_sparse_rep *V, int n,
    t_adjl_sparse_rep *A, int size_A, int *pos, int u, int v, boolean CHECK)
    /*
      like sparseg_adjl_add_edge but here we are guaranteed
      that pos + 1 < size_A
      (unless that for some reason we attempt to add
      an edge which is already there)
      
      this feature is required when A is part of a Magma block:
      we do not want to reallocate A here
      (would be done at a higher level)
 
      we check if the edge is already in the graph iff CHECK true

      also, we assume that we use this procedur only when dealing
      with an undirected graph
    */
{
    boolean   edge_added;
   
    edge_added =
        sparseg_adjl_add_dir_edge_no_extend(V, n, A, size_A, pos, u, v,
                                            CHECK);

    if (edge_added)
        sparseg_adjl_add_dir_edge_no_extend(V, n, A, size_A, pos, v, u,
                                            FALSE);

    return edge_added;
}


boolean 
sparseg_adjl_add_dir_edge (t_ver_sparse_rep *V, int n,
	t_adjl_sparse_rep **A, int *size_A, int *pos, int u, int v,
	boolean CHECK)
    /*
      add the DIRECTED edge to the sparse graph (V, n, A)
      - pos records where to add the next edge in A
      - if pos >= size_A, we must extend A
 
      we check if the edge is already in the graph iff CHECK true
    */
{
    boolean             edge_exists;

    edge_exists = FALSE;
    if (CHECK)
    {
        edge_exists = sparseg_adjl_dir_edge_exists(V, n, *A, u, v);
 
        if (edge_exists)
            return FALSE;
    }

    if (*pos == *size_A)
    {
        *size_A += 1;    /* add one directed edge */
        *A = (t_adjl_sparse_rep *)
            mem_realloc(*A, sizeof(t_adjl_sparse_rep) * *size_A);
    }
    ASSERT(*pos < *size_A);

    sparseg_adjl_add_dir_edge_no_extend(V, n, *A, *size_A, pos, u, v,
                                        FALSE);

    return TRUE;
}
 
boolean 
sparseg_adjl_add_dir_edge_no_extend (t_ver_sparse_rep *V, int n,
     t_adjl_sparse_rep *A, int size_A, int *pos, int u, int v, boolean CHECK)
    /*
      add an edge where A is guaranteed to be be big enough
      (unless that for some reason we attempt to add
      an edge which is already there)

      this feature is required when A is part of a Magma block:
      we do not want to reallocate A here
      (would be done at a higher level)
 
      we check if the edge is already in the graph iff CHECK true
    */
{
    /*
      given the way V and A represent the graph, it is simplest
      to add the new edge at the beginning of i's adj. list
    */
    int                  i_v;
    t_adjl_sparse_rep    a;

    if (CHECK && sparseg_adjl_dir_edge_exists(V, n, A, u, v))
        return FALSE;
    
    if (*pos >= size_A)
        DIE();

    /*
      otherwise always add the edge
    */
    i_v = *pos;
    a.end_vertex = v;
    a.next = V[u].first_edge;
    A[(*pos)++] = a;
    V[u].first_edge = i_v;

    return TRUE;
}
 


boolean 
sparseg_adjl_remove_edge_no_red (t_ver_sparse_rep *V, t_adjl_sparse_rep *A,
	int u, int v)
    /*
      remove the UNDIRECTED edge from sparse graph (V, A)
      if (u, v) is not an edge then nothing changes (and return FALSE)

      A will be left with "holes"
    */
{
    sparseg_adjl_remove_dir_edge_no_red(V, A, u, v);
    return sparseg_adjl_remove_dir_edge_no_red(V, A, v, u);
}
 

boolean 
sparseg_adjl_remove_dir_edge_no_red (t_ver_sparse_rep *V,
	t_adjl_sparse_rep *A, int u, int v)
    /*
      remove the DIRECTED edge from the sparse graph (V, n, A)
      if (u, v) is not an edge then nothing changes  (and return FALSE)

      A will be left with "holes"
    */
{
    int         cur_e, prev_e;

    cur_e = V[u].first_edge;
    if (cur_e == NIL)
        /*
          (u, v) is not an edge
        */
        return FALSE;

    if (A[cur_e].end_vertex == v)
    {
        V[u].first_edge = A[cur_e].next;
        return TRUE;   /* done */
    }

    while (A[cur_e].end_vertex != v)
        /*
          if (u, v) is an edge then this loop will terminate
        */
    {
        prev_e = cur_e;
        cur_e = A[cur_e].next;
        if (cur_e == NIL)
            /*
              (u, v) is not an edge
            */
            return FALSE;
    }
    ASSERT(A[cur_e].end_vertex == v);

    A[prev_e].next = A[cur_e].next;
    return TRUE;
}
 
int 
sparseg_adjl_remove_all_dir_edge_no_red (t_ver_sparse_rep *V,
	t_adjl_sparse_rep *A, int u, int v)
    /*
      remove all DIRECTED edges [u, v] from the non-simple
      sparse graph (V, n, A)
      if (u, v) is not an edge then nothing changes;
      we return the number of edges removed

      A will be left with "holes"
    */
{
    int         cur_e, prev_e, e_removed;

    if (V[u].first_edge == NIL)
        /*
          (u, v) is not an edge
        */
        return 0;

    e_removed = 0;
    while (A[V[u].first_edge].end_vertex == v)
    {
        V[u].first_edge = A[V[u].first_edge].next;
        e_removed++;

        if (V[u].first_edge == NIL)
            return e_removed;
    }
    ASSERT(A[V[u].first_edge].end_vertex != v);
    
    prev_e = V[u].first_edge;
    cur_e = A[prev_e].next;
    while (cur_e != NIL)
    {
        if (A[cur_e].end_vertex == v)
        {
            A[prev_e].next = A[cur_e].next;
            e_removed++;
            cur_e = A[cur_e].next;
        }
        else
        {
            prev_e = cur_e;
            cur_e = A[cur_e].next;
        }
    }

    return e_removed;
}
 


void 
sparseg_adjl_add_vertices (t_ver_sparse_rep **V, int n, int nmore)
    /*
      add nmore vertices
      V is assumed to have length n 
    */
{
    *V = (t_ver_sparse_rep *)
        mem_realloc(*V, sizeof(t_ver_sparse_rep) * (n + nmore));

    sparseg_adjl_add_vertices_no_extend(*V, n, nmore);
}
 
void 
sparseg_adjl_add_vertices_no_extend (t_ver_sparse_rep *V, int n, int nmore)
    /*
      add nmore vertices,
      here V is assumed to have length n + nmore (ie V has already
      been made bigger)
    */
{
    int                  v;

    for (v = n; v < n + nmore; v++)
    {
        V[v].first_edge = NIL;
    }
}
 
void 
sparseg_adjl_remove_vertex (t_ver_sparse_rep **V, int n,
	t_adjl_sparse_rep *A, int pos_A, int w, int *e)
    /*
      V is assumed to have length n: we will reallocate
      V so that V will have length n-1

      A is occupied from [0..pos-1], A will be left with holes

      we also assume that the graph can have loops and multiple edges;
      further, we the edge counting implicitely assumes that graph
      is undirected!!!

      this must be eventually fixed
    */
{
    int                  v, nv, edge, loops;
    t_ver_sparse_rep     *new_V;

    /*
      we first count the loops if any 
    */
    loops = 0;
    edge = (*V)[w].first_edge;
    while (edge != NIL)
    {
        loops = A[edge].end_vertex == w ? loops + 1 : loops;
        edge = A[edge].next;
    }
    ASSERT(loops % 2 == 0);
    loops /= 2;

    /*
      we recreate the vertices array
    */
    new_V = (t_ver_sparse_rep *)
        mem_malloc(sizeof(t_ver_sparse_rep) * (n - 1));

    for (v = 0, nv = 0; v < n; v++, nv++)
    {
        if (v == w)
        {
            nv--;
        }
        else
        {
            new_V[nv].first_edge = (*V)[v].first_edge;
        }
    }
    mem_free(*V);
    *V = new_V;

    *e -= loops;
    sparseg_adjl_remove_vertex_no_red(*V, n, A, w, e);
    
    /*
      oops! not relabelling vertices can wreck havock!
    */
    sparseg_adjl_relabel_vertex(A, pos_A, w);
}
 
void 
sparseg_adjl_remove_vertex_no_red (t_ver_sparse_rep *V, int n,
	t_adjl_sparse_rep *A, int w, int *e)
    /*
      here V has already size n - 1 and has been initialised,
      all what remains to do is to remove the edges incident
      from w in A

      A will be left with holes
    */
{
    int                  v, nbr_e_removed;

    nbr_e_removed = 0;
    for (v = 0; v < n - 1; v++)
    {
        nbr_e_removed += sparseg_adjl_remove_all_dir_edge_no_red(V, A, v, w);
    }

    *e= *e - nbr_e_removed;
}
 
void 
sparseg_adjl_relabel_vertex (t_adjl_sparse_rep *A, int pos, int u)
    /*
      relabel all vertices v > u as v-1
      (required when removing a vertex)
    */
{
    int                  i;

    for (i = 0; i < pos; i++)
    {
        A[i].end_vertex = A[i].end_vertex > u ?
            A[i].end_vertex - 1 : A[i].end_vertex;
    }
}
 
/*
 *  sparseg_adjl_pred.c
 */
 
/*
  What:
  *****
  
  Implementing:
 
  Some high-level functions on the sparse graph as
  an adjacency list: predicates.


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

#include "planarity.h"

#define IF_DEB(x)    {}
#define IF_VERB(x)   {}



/* aproto: header embed_graph_protos.h */


#ifndef PLANAR_IN_MAGMA
#endif

boolean 
sparseg_adjl_dir_edge_exists (t_ver_sparse_rep *V, int n,
	t_adjl_sparse_rep *A, int u, int v)
    /*
      does the directed edge [u, v] already exist in the graph
    */
{
    int         cur_e, prev_e;

    cur_e = V[u].first_edge;
    if (cur_e == NIL)
        return FALSE;

    if (A[cur_e].end_vertex == v)
    {
        return TRUE; 
    }

    while (A[cur_e].end_vertex != v)
    {
        prev_e = cur_e;
        cur_e = A[cur_e].next;
        if (cur_e == NIL)
            /*
              (u, v) is not an edge
            */
            return FALSE;
    }
    ASSERT(A[cur_e].end_vertex == v);
    return TRUE;
}



boolean 
sparseg_adjl_u_adj_v (t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A,
	int u, int v)
    /*
      is u adj. to v 
    */
{
    return sparseg_adjl_dir_edge_exists(V, n, A, u, v);
}


boolean 
sparseg_adjl_sub (t_ver_sparse_rep *V1, int n1, t_adjl_sparse_rep *A1,
	t_ver_sparse_rep *V2, int n2, t_adjl_sparse_rep *A2)
    /*
      test if the (V1, n1, A1) sparse graph is a subgraph of
      the (V2, n2, A2) graph
    */
{
    int             v, *fp, n, bign, i;

    n = n1 > n2 ? n2 : n1;
    bign = n1 > n2 ? n1 : 0;
    fp = (int *) mem_malloc(sizeof(int) * n);
    for (i = 0; i < n; i++)
        fp[i] = NIL;

    for (v = 0; v < n; v++)
    {
        int      ne1, ne2;

        ne1 = V1[v].first_edge;
        ne2 = V2[v].first_edge;
        if (ne1 == NIL)
        {
            continue;
        }
        else if (ne2 == NIL)
        {
            mem_free(fp);
            return FALSE;
        }
        
        while (ne2 != NIL)
        {
            int u2;

            u2 = A2[ne2].end_vertex;
            fp[u2] = v;
            ne2 = A2[ne2].next;
        }

        while (ne1 != NIL)
        {
            int u1;

            u1 = A1[ne1].end_vertex;
            if (fp[u1] != v)
            {
                mem_free(fp);
                return FALSE;
            }
            ne1 = A1[ne1].next;
        }
    }
    mem_free(fp);

    for (v = n; v < bign; v++)
        /*
          those vertices must not be end points of edges:
          this chcek is only necessary in the digraph case
        */
    {
        if (V1[v].first_edge != NIL)
            return FALSE;
    }
    
    return TRUE;
}



boolean 
sparseg_adjl_eq (t_ver_sparse_rep *V1, int n1, t_adjl_sparse_rep *A1,
	t_ver_sparse_rep *V2, int n2, t_adjl_sparse_rep *A2)
    /*
      compare the two sparse graphs (V1, n1, A1) & (V2, n2, A2)
      we don't know their number of edges
    */
{
    if (n1 != n2)
        return FALSE;

    return sparseg_adjl_sub(V1, n1, A1, V2, n2, A2)
        && sparseg_adjl_sub(V2, n2, A2, V1, n1, A1);
}



/*
 *  sparseg_dlcl_misc.c
 */
 
/*
  What:
  *****
  
  Implementing:
 
  Housekeeping for an internal sparse graph representation
  internal to the planarity tester and obstruction isolator.

  This sparse graph consists of an array of doubly linked circular lists
  (the neighbour lists for each vertex).
 
 
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********
 
  Paulette Lieby (Magma), Brendan McKay (ANU)
 
  Started October 2001
*/

#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
 
 
/* aproto: header embed_graph_protos.h */

/* aproto: beginstatic -- don't touch this!! */
static boolean sparseg_dlcl_is_present (t_dlcl *, int, t_dlcl **);
/* aproto: endstatic -- don't touch this!! */
 
 
#ifndef PLANAR_IN_MAGMA
#endif
 

void 
sparseg_dlcl_delete (t_dlcl **g, int n)
{
    int      i;

    for (i = 0; i < n; i++)
    {
       embedg_dlcl_delete(g[i]);
    }
    mem_free(g);
}

void 
sparseg_dlcl_print (t_dlcl **g, int n)
{
    int      i;

    for (i = 0; i < n; i++)
    {
        fprintf(stdout,"%d:\t", i);
        embedg_dlcl_print(g[i]);
    }
}


static boolean 
sparseg_dlcl_is_present (t_dlcl *l, int label, t_dlcl **p)
{
    *p = embedg_dlcl_find(l, label);
    return *p == NP ? FALSE : TRUE;
}


boolean 
sparseg_dlcl_is_adjacent (t_dlcl **g, int n, int v, int u, t_dlcl **p)
    /*
      is u adjacent to v
    */
{
    ASSERT(v >= 0 && v < n && u >= 0 && u < n);
    return sparseg_dlcl_is_present(g[v], u, p);
}

void 
sparseg_dlcl_append_to_neigh_list (t_dlcl **g, int n, int v, int u, int in_adjl)
    /*
      append u to the neighbour list of v
    */
{
    t_dlcl   *u_rec;
 
    u_rec = embedg_dlcl_rec_new(u);
    u_rec->in_adjl = in_adjl;
    g[v] = embedg_dlcl_rec_append(g[v], u_rec);
}




void 
sparseg_dlcl_to_sparseg (t_dlcl **g, int n, int e,
	t_ver_sparse_rep **V, t_adjl_sparse_rep **A)
    /*
      e is the number of undirected edges of g

      convert a dlcl into the standard sparseg rep. as an
      adjacency list
    */
{
    int                 i_e, v;

    *V = (t_ver_sparse_rep *) mem_malloc(sizeof(t_ver_sparse_rep) * n);
    *A = (t_adjl_sparse_rep *) mem_malloc(sizeof(t_adjl_sparse_rep) * 2 * e);
 
    for (v = 0; v < n; v++)
        (*V)[v].first_edge = NIL;
    
    i_e = 0;
    for (v = 0; v < n; v++)
    {
        t_dlcl     *l, *p;

        l = p = g[v];
        if (!embedg_dlcl_is_empty(p))
        {
            t_adjl_sparse_rep   a;
            
            ASSERT((*V)[v].first_edge == NIL);
            (*V)[v].first_edge = i_e;
            a.end_vertex = p->info;
            a.next = i_e + 1;
            (*A)[i_e++] = a;
            
            p = embedg_dlcl_list_next(p);
            while (p != l)
            {
                a.end_vertex = p->info;
                a.next = i_e + 1;
                (*A)[i_e++] = a;

                p = embedg_dlcl_list_next(p);
            }

            /*
              end of list for v
            */
            (*A)[i_e - 1].next = NIL;
        }
    }
    ASSERT(i_e == 2 * e);
}

boolean 
sparseg_dlcl_sub (t_dlcl **g1, int n1, t_dlcl **g2, int n2)
    /*
      is g2 a subgraph of g1

      I request that both graphs have same order

      This is not used anywhere... do we need it???
    */
{
    int           n, v, *fp;

    if (n1 != n2)
        return FALSE;

    n = n1;
    fp = (int *) mem_malloc(sizeof(int) * n);
    for (v = 0; v < n; v++)
        fp[v] = NIL;

    for (v = 0; v < n; v++)
    {
         t_dlcl     *l1, *p1, *l2, *p2;
 
        l1 = p1 = g1[v];
        l2 = p2 = g2[v];
        if (embedg_dlcl_is_empty(p1) && !embedg_dlcl_is_empty(p2))
        {
            mem_free(fp);
            return FALSE;
        }
        if (embedg_dlcl_is_empty(p2))
        {
            continue;
        }
        
        fp[p1->info] = v;
        p1 = embedg_dlcl_list_next(p1);
        while (p1 != l1)
        {
            fp[p1->info] = v;
            p1 = embedg_dlcl_list_next(p1);
        }

        if (fp[p2->info] != v)
        {
            mem_free(fp);
            return FALSE;
        }
        p2 = embedg_dlcl_list_next(p2);
        while (p2 != l2)
        {
            if (fp[p2->info] != v)
            {
                mem_free(fp);
                return FALSE;
            }
        }
    }
    mem_free(fp);
    
    return TRUE;
}
/*
 *  VES_misc.c
 */
 
/*
  What:
  *****
  
  Implementing:

  All low-level routines for the VES structure:

  - the VES structure is solely used within the planarity tester
    and obstruction isolator

  - it stores vertices, virtual vertices and edges
       --more on this later--

  - it allows for circular doubly linked lists, hence
    enabling us -among other things- to store the
    graph embedding if the tester is successful

  - basic features:
    + the VES has exactly size 2n + 2(3n-5) :
      we add at most one more edge than the max for a planar graph
      (need to x by 2: we store directed edges)
    + a vertex and the edges incident FROM it are linked in a doubly
      linked circular list
    + where a vertex is inserted between two of its outcoming edges
      determines an external face walk for a bicomponent
    + the twin edge is more commonly known as the inverse edge
    + we have tree and back edges (from the DFS), and short-cut edges
      which are added by the tester
      -but short-cut edges are added in such a way as to maintain
       planarity (in a local sense)
    + vertices and edges can be marked (visited for example)
    + they have an orientation which must be eventuall recovered
      and which is set in the merge_bicomp routine
    + vertices are essentially known via their DFI or DFS index
      (though their label is stored too)

     blah, blah.... later then.
     Have a look at embedg_planar_alg_init which initialises the VES
     structure

     
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_DEB_SCE(x)    {}
#define IF_DEB_PROPER_FACE(x) {}
#define IF_VERB(x)   {}
 
 
/* aproto: header embed_graph_protos.h */

boolean 
embedg_VES_is_vertex (int n, int i)
    /*
      is this a vertex
      (relative to the "big" array of size 2n + 2(3n-5))
    */
{
    return i < n ? TRUE : FALSE;
}

boolean 
embedg_VES_is_virtual_vertex (int n, int i)
    /*
      is this a virtual vertex
      (relative to the "big" array of size 2n + 2(3n-5))

      a virtual vertex is a vertex v^c which denotes the
      DFS parent of the child c

      see embedg_planar_alg_init for more
    */
{
    return i >= n && i < 2*n ? TRUE : FALSE;
}

boolean 
embedg_VES_is_edge (int n, int i)
    /*
      is this an edge
      (relative to the "big" array of size 2n + 2(3n-5))
    */
{
    return i >= 2*n ? TRUE : FALSE;
}

boolean 
embedg_VES_is_tree_edge (t_ver_edge *embed_graph, int n, int i)
    /*
      is this s tree edge
    */
{
    return embedg_VES_is_edge(n, i)
        && embed_graph[i].type == TE;
}

boolean 
embedg_VES_is_back_edge (t_ver_edge *embed_graph, int n, int i)
    /*
      is this a back edge
    */
{
    return embedg_VES_is_edge(n, i)
        && embed_graph[i].type == BE;
}

boolean 
embedg_VES_is_short_cut_edge (t_ver_edge *embed_graph, int n, int i)
    /*
      as the name indicates...
    */
{
    return embedg_VES_is_edge(n, i)
        && embed_graph[i].type == SCE;
}

void 
embedg_VES_print_vertex (int n, int v)
{
    ASSERT(embedg_VES_is_vertex(n, v));
    fprintf(stdout, "%d  ", v);
}

void 
embedg_VES_print_virtual_vertex (t_ver_edge *embed_graph, int n, int v)
{
    int          c;
    
    ASSERT(embedg_VES_is_virtual_vertex(n, v));
    c = v - n;
    fprintf(stdout, "%d^%d  ", embed_graph[c].DFS_parent, c);
}

void 
embedg_VES_print_any_vertex (t_ver_edge *embed_graph, int n, int v)
{
    if (embedg_VES_is_vertex(n, v))
    {
        embedg_VES_print_vertex(n, v);
    }
    else
    {
        embedg_VES_print_virtual_vertex(embed_graph, n, v);
    }
}

void 
embedg_VES_print_any_rec (t_ver_edge *embed_graph, int n, int r)
{
    if (embedg_VES_is_edge(n, r))
    {
        embedg_VES_print_edge(embed_graph, n, r);
    }
    else
    {
        embedg_VES_print_any_vertex(embed_graph, n, r);
    }
}

void 
embedg_VES_print_edge (t_ver_edge *embed_graph, int n, int e)
{
    int          v, prev, cur;
    
    ASSERT(embedg_VES_is_edge(n, e));

    /*
      must find the vertex in the doubly linked circular list
      of vertices/edges
    */

    prev = e;
    cur = v = embed_graph[e].link[0];
    if (embedg_VES_is_vertex(n, v)
        || embedg_VES_is_virtual_vertex(n, v))
    {
        embedg_VES_print_any_vertex(embed_graph, n, v);
        fprintf(stdout, ", ");
        embedg_VES_print_any_vertex(embed_graph, n,
                                         embed_graph[e].neighbour);
        fprintf(stdout, "):0\n");
    }
    else while (!embedg_VES_is_vertex(n, v)
                && !embedg_VES_is_virtual_vertex(n, v))
    {
        v = embedg_VES_get_next_in_dlcl(embed_graph, n,
                                                   cur, prev);

        if (embedg_VES_is_vertex(n, v)
            || embedg_VES_is_virtual_vertex(n, v))
        {
            embedg_VES_print_any_vertex(embed_graph, n, v);
            fprintf(stdout, ", ");
            embedg_VES_print_any_vertex(embed_graph, n,
                                             embed_graph[e].neighbour);
            fprintf(stdout, "):0\n");
        }
        else
        {
            prev = cur;
            cur = v;
        }
    }
}

void 
embedg_VES_print_flipped_edges (t_ver_edge *embed_graph, int n, int edge_pos)
    /*
      print those edges in the structure whose sign is CLOCKW,
      ie which have been flipped at some stage
    */
{
    int          e;
    
    for (e = 2*n; e <= edge_pos; e++)
    {
        if (!embedg_VES_is_short_cut_edge(embed_graph, n, e))
            /*
              we don't care about the short-cut edges
            */
        {
            if (embed_graph[e].sign != CCLOCKW)
            {
                embedg_VES_print_edge(embed_graph, n, e);
            }
        }
    }
}

#if 0
int 
embedg_VES_get_edge_from_ver (t_ver_edge *embed_graph, int n, int v)
    /*
      not used anywhere; why is this here???
    */
{
    int          in, e;

    ASSERT(embedg_VES_is_vertex(n, v)
           || embedg_VES_is_virtual_vertex(n, v));
    
    in = embedg_VES_is_edge(n, embed_graph[v].link[0]) ? 0 : 1;
    e = embed_graph[v].link[in];
    ASSERT(embedg_VES_is_edge(n, e));

    return e;
}

int 
embedg_VES_get_ver_from_edge (t_ver_edge *embed_graph, int n, int e)
{
    int          in, v;

    ASSERT(embedg_VES_is_edge(n, e));

    in = embedg_VES_is_vertex(n, embed_graph[e].link[0])    
        || embedg_VES_is_virtual_vertex(n, embed_graph[e].link[0])
        ?
        0 : 1;

    v = embed_graph[e].link[in];
    ASSERT(embedg_VES_is_vertex(n, v)
           || embedg_VES_is_virtual_vertex(n, v));

    return v;
}
#endif

int 
embedg_VES_get_twin_edge (t_ver_edge *embed_graph, int n, int e)
    /*
      the twin edge is understood as being the inverse edge
    */
{
    int          twin;

    ASSERT(embedg_VES_is_edge(n, e));

    twin = e % 2 == 0 ? e + 1 : e - 1;
    ASSERT(embedg_VES_is_edge(n, twin));

    return twin;
}

int 
embedg_VES_get_ver_from_virtual (t_ver_edge *embed_graph, int n, int vv)
    /*
      get v from the virtual vertex v^c
    */
{
    int          v;

    ASSERT(embedg_VES_is_virtual_vertex(n, vv));
    v = embed_graph[vv - n].DFS_parent;

    return v;
}

int 
embedg_VES_get_ver (t_ver_edge *embed_graph, int n, int v)
{
    if (embedg_VES_is_virtual_vertex(n, v))
        return embedg_VES_get_ver_from_virtual(embed_graph, n, v);

    return v;
}
    

int 
embedg_VES_get_next_in_dlcl (t_ver_edge *embed_graph, int n, int r, int prev)
    /*
      r is a (virtual) vertex or edge record in embed_graph:
      get the next in the list (formed by the .link[] fields)
      in the doubly linked circular list

      so that prev != next
      -- NOTE: a priori these lists always contain 2 elts at least
      so that there shouldn't be any problem...
      --> huh? is that true?
    */
{
    return embed_graph[r].link[0] == prev ?
        embed_graph[r].link[1] : embed_graph[r].link[0];
}


void 
embedg_VES_walk_bicomp (t_ver_edge *embed_graph, int n, int v, int vin)
    /*
      walk the external face of the bicomp starting
      at VIRTUAL vertex v entered via vin

      this of course assumes that the "thing" rooted at
      v is a bicomponent -- depending where we are at in the
      tester this is not necessarily the case
      -- I comment upon this in merge_bicomps.c:
               embedg_VES_merge_pertinent_bicomps
    */
{
    int        start, startin, s, sin;

    ASSERT(embedg_VES_is_virtual_vertex(n, v));

    embedg_VES_print_virtual_vertex(embed_graph, n, v);

    s = NIL;
    start = v;
    startin = vin;
    while (s != v)
    {
        embedg_VES_get_succ_on_ext_face(embed_graph, n, start, startin,
                                             FALSE, 0, &s, &sin);
        if (embedg_VES_is_virtual_vertex(n, s))
        {
            embedg_VES_print_virtual_vertex(embed_graph, n, s);
        }
        else
        {
            embedg_VES_print_vertex(n, s);
        }
        start = s;
        startin = sin;
    }
    fprintf(stdout, "\n");
}

void 
embedg_VES_print_adj_list (t_ver_edge *embed_graph, int n, int r,
	boolean consistent)
    /*
      print r's adjacency list - r can be a vertex or edge

      the boolean <consistent> if true assumes that
      the list is consistent (will determine the way we traverse the list)

      a priori we should get the same result either way
    */
{
    if (consistent)
    {
        int        next;

        embedg_VES_print_any_rec(embed_graph, n, r);
        
        next = embed_graph[r].link[0];
        while (next != r)
        {
            embedg_VES_print_any_rec(embed_graph, n, next);
            next = embed_graph[next].link[0];
        }
    }
    else
    {
        int          prev, cur, next;
        
        embedg_VES_print_any_rec(embed_graph, n, r);
        
        prev = r;
        cur = embed_graph[r].link[0];
        
        while (cur != r)
        {
            embedg_VES_print_any_rec(embed_graph, n, cur);
            next = embedg_VES_get_next_in_dlcl(embed_graph, n,
                                                          cur, prev);
            prev = cur;
            cur = next;
        }
    }
}

boolean 
embedg_VES_is_adj_list_consistent (t_ver_edge *embed_graph, int n, int r)
    /*
      checks that r's adjacency list is consistent:
      ie, that either traversing it using link[0] always
      or traversing it using embedg_VES_get_next_in_dlcl
      gives the SAME result
    */
{
    int          *list_link, *list_n_dldl, il, id, i;

    list_link = (int *) mem_malloc(sizeof(int) * 2 * n);
    list_n_dldl = (int *) mem_malloc(sizeof(int) * 2 * n);
    /*
      must allocate 2*n space: I could have TE and SCE with same neighbour
      (or BE and SCE as well)
    */
    il = id = -1;

    /*
      traversing the list via link[0]
    */
    {
        int        next;

        list_link[++il] = r;
        
        next = embed_graph[r].link[0];
        while (next != r)
        {
            list_link[++il] = next;
            next = embed_graph[next].link[0];
        }
    }

    /*
       traversing the list using embedg_VES_get_next_in_dlcl
    */
    {
        int          prev, cur, next;
        
        list_n_dldl[++id] = r;
        prev = r;
        cur = embed_graph[r].link[0];
        
        while (cur != r)
        {
            list_n_dldl[++id] = cur;
            next = embedg_VES_get_next_in_dlcl(embed_graph, n,
                                                          cur, prev);
            prev = cur;
            cur = next;
        }
    }

    if (il != id)
    {
        mem_free(list_link);
        mem_free(list_n_dldl);
        return FALSE;
    }

    for (i = 0; i <= il; i++)
    {
        if (list_link[i] != list_n_dldl[i])
        {
            mem_free(list_link);
            mem_free(list_n_dldl);
            return FALSE;
        }
    }

    mem_free(list_link);
    mem_free(list_n_dldl);
    return TRUE;
}


boolean 
embedg_VES_are_adj_lists_consistent (t_ver_edge *embed_graph, int n)
    /*
      checks that the adjacency list of each vertex is consistent
      in the manner of embedg_VES_is_adj_list_consistent
    */
{
    int          i;

    /*
      it is enough to visit the vertices and virtual vertices only
      (I don't think it is enough to do the vertices only --??)
    */
    for (i = 0; i < 2*n; i++)
        if (!embedg_VES_is_adj_list_consistent(embed_graph, n, i))
            return FALSE;

    return TRUE;
}



void 
embedg_VES_remove_edge (t_ver_edge *embed_graph, int n, int e)
    /*
      remove edge e from the embedding
    */
{
    int          r1, r2, r1out, r2in, twin;

    ASSERT(embedg_VES_is_edge(n, e));

    IF_DEB_SCE(
               fprintf(stdout, "removing an SCE, enter\n");
               embedg_VES_print_edge(embed_graph, n, e);
               )
    
    r1 = embed_graph[e].link[0];
    r2 = embed_graph[e].link[1];

    /*
      disable e and link r1 and r2 together:
      we had r1 -> e -> r2
    */
    embed_graph[e].link[0] = embed_graph[e].link[1] = e;
    
    r1out = embed_graph[r1].link[0] == e ? 0 : 1;
    r2in = embed_graph[r2].link[0] == e ? 0 : 1;

    if (r1 == r2)
        /*
          this I think should never happen, but one never knows...
        */
    {
        embed_graph[r1].link[0] = embed_graph[r1].link[1] = r1;
    }
    else
    {
        embed_graph[r1].link[r1out] = r2;
        embed_graph[r2].link[r2in] = r1;
    }

    ASSERT(embedg_VES_is_adj_list_consistent(embed_graph, n, r1));

    /*
      now we must do a similar thing for the twin
      (which must get reomved as well)
    */
    twin = embedg_VES_get_twin_edge(embed_graph, n, e);

    IF_DEB_SCE(
               fprintf(stdout, "removing an SCE, the twin\n");
               embedg_VES_print_edge(embed_graph, n, twin);
               )
    
    r1 = embed_graph[twin].link[0];
    r2 = embed_graph[twin].link[1];

    embed_graph[twin].link[0] = embed_graph[twin].link[1] = twin;

    r1out = embed_graph[r1].link[0] == twin ? 0 : 1;
    r2in = embed_graph[r2].link[0] == twin ? 0 : 1;

    if (r1 == r2)
    {
        embed_graph[r1].link[0] = embed_graph[r1].link[1] = r1;
    }
    else
    {
        embed_graph[r1].link[r1out] = r2;
        embed_graph[r2].link[r2in] = r1;
    }

    ASSERT(embedg_VES_is_adj_list_consistent(embed_graph, n, r1));
}


void 
embedg_VES_set_orientation (t_ver_edge *embed_graph, int n, int *ver_orient)
    /*
      using the vertices' orientation as given in ver_orient
      we set the orientation for each edge in the adjacency list
      for each vertex

      to do this we use the field sign which is NOT needed
      anymore by the tester since by the time we call this
      function we would have finished with that bit (the tester)

      sign is only set when merging bicomps
      - even though we'll perform another walkdown when
      recovering an obstruction (if any) no bicomp merging will occur,
      so we are safe
    */
{
    int          v;

    for (v = 0; v < n; v++)
    {
        int      o, e;

        o = ver_orient[v];
        embed_graph[v].sign = o;

        e = embed_graph[v].link[0];

        while (e != v)
            /*
              just as a note: note the way I get the next in the list
              here (as opposed to using
              embedg_VES_get_next_in_dlcl):
              this is because I implicitely assume that
              the adjacency lists are consistent
              
              Also note that edges can be SCE, it doesn't really matter
              anyway (they may not have been removed yet
              -- see the way we recover the obstruction:
              embedg_mark_obstruction)
            */
        {
            embed_graph[e].sign = o;
            e = embed_graph[e].link[0];
        }       
    }
}


/*
 *  dlcl_misc.c
 */
 
/*
  What:
  *****
  
  Implementing:

  Housekeeping for a simple doubly linked circular list:
  this is a data structure ONLY used WITHIN
  the planarity tester and obstruction isolator and is not to be
  confused with the VES structure mentionned elsewhere.

  The VES structure is an array, while the dlcl one is a list of
  pointers.

  The  dlcl is especially useful as it allows for the storage
  of an ordered list.


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/


#include "planarity.h"

#define IF_DEB(x)    {}
#define IF_VERB(x)   {}


/* aproto: header embed_graph_protos.h */

/* aproto: beginstatic -- don't touch this!! */
static void embedg_dlcl_rec_free (t_dlcl *);
static void embedg_dlcl_rec_insert_right (t_dlcl *, t_dlcl *);
static void embedg_dlcl_rec_insert_left (t_dlcl *, t_dlcl *);
static void embedg_dlcl_rec_retrieve (t_dlcl *);
static void embedg_dlcl_rec_delete (t_dlcl *);
static boolean embedg_dlcl_is_singleton (t_dlcl *);
/* aproto: endstatic -- don't touch this!! */


#ifndef PLANAR_IN_MAGMA
#endif


t_dlcl *
embedg_dlcl_rec_new (int info)
    /*
      create a new record with info <info> in the global array
      to insert in the list
    */
{
    t_dlcl    *r;
 
    r = (t_dlcl *) mem_malloc(sizeof(t_dlcl));
    r->info = info;
    r->in_adjl = r->twin_in_adjl = NIL;
    r->mult = 1;
    r->right = r;
    r->left = r;
    return r;
}

static void 
embedg_dlcl_rec_free (t_dlcl *r)
    /*
      free
    */
{
    mem_free(r);
}

void 
embedg_dlcl_rec_print (t_dlcl *r)
{
    fprintf(stdout,"%d ", r->info);
}

void 
embedg_dlcl_print (t_dlcl *l)
{
    t_dlcl    *p = l;
    
    if (!embedg_dlcl_is_empty(p))
    {
        embedg_dlcl_rec_print(p);
        p = embedg_dlcl_list_next(p);
        while (p != l)
        {
            embedg_dlcl_rec_print(p);
            p = embedg_dlcl_list_next(p);
        }
    }
    fprintf(stdout,"\n");
}


static void 
embedg_dlcl_rec_insert_right (t_dlcl *l, t_dlcl *r)
{
    t_dlcl    *tmp_r, *tmp_l;

    tmp_r = l->right;
    tmp_l = tmp_r->left;

    l->right = r;
    r->right = tmp_r;

    r->left = tmp_l;
    tmp_r->left = r;
}


static void 
embedg_dlcl_rec_insert_left (t_dlcl *l, t_dlcl *r)
{
    t_dlcl    *tmp_r, *tmp_l;

    tmp_l = l->left;
    tmp_r = tmp_l->right;

    l->left = r;
    r->left = tmp_l;

    r->right = tmp_r;
    tmp_l->right = r;
}

t_dlcl *
embedg_dlcl_rec_append (t_dlcl *l, t_dlcl *r)
{
    if (embedg_dlcl_is_empty(l))
        return r;
    
    embedg_dlcl_rec_insert_left(l, r);
    return l;
}

t_dlcl *
embedg_dlcl_rec_prepend (t_dlcl *l, t_dlcl *r)
{
    if (embedg_dlcl_is_empty(l))
        return r;

    embedg_dlcl_rec_insert_left(l, r);
    return r;
}

t_dlcl *
embedg_dlcl_cat (t_dlcl *l, t_dlcl *m)
    /*
      concatenate m to the RIGHT of the end of l
      WITHOUT copying m 
    */
{
    t_dlcl    *h1, *h2, *e1, *e2;

    if (embedg_dlcl_is_empty(l))
        return m;
    if (embedg_dlcl_is_empty(m))
        return l;
    
    h1 = l;
    e1 = l->left;
    h2 = m;
    e2 = m->left;

    e1->right = h2;
    h2->left = e1;
    e2->right = h1;
    h1->left = e2;

    return l;
}

t_dlcl *
embedg_dlcl_find (t_dlcl *l, int info)
{
    t_dlcl    *p = l;
    
    if (!embedg_dlcl_is_empty(p))
    {
        if (p->info == info)
        {
            return p;
        }
        p = embedg_dlcl_list_next(p);
        while (p != l)
        {
            if (p->info == info)
            {
                return p;
            }
            p = embedg_dlcl_list_next(p);
        }
    }
    return NP;
}

t_dlcl *
embedg_dlcl_find_with_NIL_twin_in_adjl (t_dlcl *l, int info)
{
    t_dlcl    *p = l;
    
    if (!embedg_dlcl_is_empty(p))
    {
        if (p->info == info && p->twin_in_adjl == NIL)
        {
            return p;
        }
        p = embedg_dlcl_list_next(p);
        while (p != l)
        {
            if (p->info == info && p->twin_in_adjl == NIL)
            {
                return p;
            }
            p = embedg_dlcl_list_next(p);
        }
    }
    return NP;
}



static void 
embedg_dlcl_rec_retrieve (t_dlcl *r)
{
    t_dlcl    *right, *left;
 
    right = r->right;
    left = r->left;
 
    left->right = right;
    right->left = left;
 
    r->right = r;
    r->left = r;
}

static void 
embedg_dlcl_rec_delete (t_dlcl *r)
{
    embedg_dlcl_rec_retrieve(r);
    embedg_dlcl_rec_free(r);
}


t_dlcl *
embedg_dlcl_delete_first (t_dlcl *l)
    /*
      prune the list from the head:
      - set new head to right of old head
      - delete old head
    */
{
    t_dlcl    *new_head;

    ASSERT(!embedg_dlcl_is_empty(l));
    if (embedg_dlcl_is_singleton(l))
    {
        new_head = NP;
    }
    else
    {
        new_head = l->right;
    }
    embedg_dlcl_rec_delete(l);
    return new_head;
}


t_dlcl *
embedg_dlcl_delete_rec (t_dlcl *l, t_dlcl *r)
    /*
      delete r from l;
      if r == l, set new head to right of old head
    */
{
    if (r == l)
    {
        return embedg_dlcl_delete_first(l);
    }
    embedg_dlcl_rec_delete(r);
    return l;
}


boolean 
embedg_dlcl_is_empty (t_dlcl *l)
{
    return (l == NP) ? TRUE : FALSE;
}


static boolean 
embedg_dlcl_is_singleton (t_dlcl *l)
{
    return (l->right == l) ? TRUE : FALSE;
    /*
      same as l->left == l
    */
}

t_dlcl *
embedg_dlcl_list_next (t_dlcl *l)
    /*
      this assumes no choice in the direction of the walking
      (always to the right)
      -- good enough when deleting for example or when
      the direction of the walking does not matter
    */
{
    return l->right;
}
 

t_dlcl *
embedg_dlcl_list_prev (t_dlcl *l)
    /*
      this assumes no choice in the direction of the walking
      (always to the right)
    */
{
    return l->left;
}
 
t_dlcl *
embedg_dlcl_list_last (t_dlcl *l)
{
    return embedg_dlcl_list_prev(l);
}
 


void 
embedg_dlcl_delete (t_dlcl *l)
{
    if (!embedg_dlcl_is_empty(l))
    {
        while (!embedg_dlcl_is_singleton(l))
        {
            t_dlcl    *next;
 
            next = embedg_dlcl_list_next(l);
            embedg_dlcl_rec_delete(next);
        }
        embedg_dlcl_rec_delete(l);
    }
}

t_dlcl *
embedg_dlcl_copy (t_dlcl *l)
{
    t_dlcl    *p, *c;
    
    if (embedg_dlcl_is_empty(l))
        return NP;

    c = embedg_dlcl_rec_new(l->info);
    
    p = embedg_dlcl_list_next(l);
    while (p != l)
    {
        t_dlcl     *temp;

        temp = embedg_dlcl_rec_new(p->info);
        temp->in_adjl = p->in_adjl;
        temp->twin_in_adjl = p->twin_in_adjl;
        temp->mult = p->mult;
        c = embedg_dlcl_rec_append(c, temp);
        p = embedg_dlcl_list_next(p);
    }
    return c;
}


int 
embedg_dlcl_length (t_dlcl *l)
{
    t_dlcl    *p;
    int       n;
    
    if (embedg_dlcl_is_empty(l))
        return 0;

    p = embedg_dlcl_list_next(l);
    n = 1;
    while (p != l)
    {
        n++;
        p = embedg_dlcl_list_next(p);
    }
    return n;
}
/*
 *  planar_by_edge_addition.c
 */
 
/*
  What:
  *****
  
  Implementing:

  The top level for the planarity tester.


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
#define IF_DEB_TREE(x)    {}
#define IF_DEB_EDGES(x) {}
#define IF_CPU(x) {}
 
 
/* aproto: header embed_graph_protos.h */


#ifndef PLANAR_IN_MAGMA
#endif
 

boolean 
sparseg_adjl_is_planar (
    t_ver_sparse_rep *V,
    int n,
    t_adjl_sparse_rep *A,        /* input sparse graph */
    int *nbr_c,        /* size of the graph, #components
                                    */
    t_dlcl ***dfs_tree,      /* a sparse graph rep. for the dfs tree
                                      -- vertices are as DFIs
                                      -- and children are ordered wrt
                                         lowpoint value
                                   */
    t_dlcl ***back_edges,    /* for each vertex v, a dlcl
                                      of the back edges [v, x] incident to v
                                      where x is a DESCENDANT of v
                                      (vertices are given as DFIs)
                                   */
    t_dlcl ***mult_edges,    /* for each vertex v, a dlcl
                                      of the back edges [v, x] incident to v
                                      where x is a DESCENDANT of v
                                      (vertices are given as DFIs)
                                   */
    t_ver_edge **embed_graph,    /* output graph embedding -- more on that
                                      later
                                   */
    int *edge_pos,        /* pos. in embed_graph for addition
                                      of the next edge */
    int *vr,
    int *wr         /* if graph is non planar, return
                                      the unembedded edge
                                      (where wr descendant of vr)
                                   */
)
    /*
      as the name indicates: is the graph planar?
    */
{
    int          v;

    IF_CPU(
    float      sttime; float time_to_now;
    )


    *embed_graph =
        embedg_planar_alg_init(V, n, A, nbr_c,
                               edge_pos, dfs_tree, back_edges, mult_edges);
    IF_CPU(
    sttime = time_current_user();
          )

    for (v = n - 1; v >= 0; v--)
        /*
          visit all vertices in descending DFI order
        */
    {
        t_dlcl     *be_l, *te_l, *p;

        IF_DEB(
               fprintf(stdout, "top level, vertex   %d\n", v);
               )
            
        /*
          find all the back edges [w, v] where w is a descendant of v
          and perform a walkup from w to v
          (ie determine which bicomps are pertinent)
        */
        be_l = (*back_edges)[v];
        p = be_l;

        if (!embedg_dlcl_is_empty(p))
        {
            int       w;
            
            w = p->info;
            IF_DEB(
                   fprintf(stdout, "top level, before walkup for w %d\n", w);
                   )
            embedg_walkup(*embed_graph, n, v, p);

            p = embedg_dlcl_list_next(p);
            while (p != be_l)
            {
                w = p->info;
                IF_DEB(
                       fprintf(stdout, "top level, before walkup for w %d\n", w);
                       )
                embedg_walkup(*embed_graph, n, v, p);
                
                p = embedg_dlcl_list_next(p);
            }
        }

        /*
          perform a walkdown for each tree edge [v, c], c a descendant of v
          (ie attempt to embed all back edges on the pertinent bicomps)
        */
        te_l = (*dfs_tree)[v];
        p = te_l;

        if (!embedg_dlcl_is_empty(p))
        {
            int             c, vv;
            t_merge_queue   q;

            c = p->info;
            vv = c + n;
            IF_DEB(
                   fprintf(stdout, "top level, before walkdown for c %d\n", c);
                   )
            q = embedg_walkdown(*embed_graph, n, edge_pos, vv);

            IF_DEB(
                   fprintf(stdout, "top level, after walkdown for c %d, state of edges'sign\n", c);
                   embedg_VES_print_flipped_edges(*embed_graph,
                                                       n, *edge_pos);
                   )
 
            /*
              temp only
            */
            embedg_merge_queue_delete(q);
            p = embedg_dlcl_list_next(p);
            while (p != te_l)
            {
                c = p->info;
                vv = c + n;
                IF_DEB(
                       fprintf(stdout, "top level, before walkdown for c %d\n", c);
                       )
                q = embedg_walkdown(*embed_graph, n, edge_pos, vv);

                IF_DEB(
                       fprintf(stdout, "top level, after walkdown for c %d, state of edges'sign\n", c);
                       embedg_VES_print_flipped_edges(*embed_graph,
                                                           n, *edge_pos);
                       )
 
                /*
                  temp only
                */
                embedg_merge_queue_delete(q);
                
                p = embedg_dlcl_list_next(p);
            }
        }


        /*
          check that each back edge [w, v], w a descendant of v,
          has been embedded
        */
        be_l = (*back_edges)[v];
        p = be_l;

        if (!embedg_dlcl_is_empty(p))
        {
            int       w;
            
            w = p->info;
            IF_DEB(
                   fprintf(stdout, "top level, before checking embedding for w %d\n",
                           w);
                   )
            if ((*embed_graph)[w].adjacent_to == v)
                /*
                  this edge hasn't been embedded:
                  the graph is non-planar
                */
            {
                /*
                  before returning we really want to ensure that
                  the vertices' adjacency lists are consistent
                */
                ASSERT(embedg_VES_are_adj_lists_consistent(
                                                       *embed_graph, n));

                IF_CPU(
                       fprintf(stdout, "CPU for tester only %f\n",
                               (time_current_user() - sttime));
                       )

                *vr = v;
                *wr = w;
                return FALSE;
            }

            p = embedg_dlcl_list_next(p);
            while (p != be_l)
            {
                w = p->info;
                IF_DEB(
                       fprintf(stdout, "top level, before checking embedding for w %d\n",
                               w);
                       )
                if ((*embed_graph)[w].adjacent_to == v)
                {
                    /*
                      before returning we really want to ensure that
                      the vertices' adjacency lists are consistent
                    */
                    ASSERT(embedg_VES_are_adj_lists_consistent(
                                                       *embed_graph, n));

                    IF_CPU(
                           fprintf(stdout, "CPU for tester only %f\n",
                                   (time_current_user() - sttime));
                           )

                    *vr = v;
                    *wr = w;
                    return FALSE;
                }
                
                p = embedg_dlcl_list_next(p);
            }
        }
    }
    IF_DEB_EDGES(
                 fprintf(stdout, "top level, total number of edges in embedding %d\n",
                         *edge_pos - 2 * n + 1);
                 )


    /*
      before returning we really want to ensure that
      the vertices' adjacency lists are consistent
    */
    ASSERT(embedg_VES_are_adj_lists_consistent(*embed_graph, n));

    IF_CPU(
           fprintf(stdout, "CPU for tester only %f\n",
                   (time_current_user() - sttime));
           )

    return TRUE;    
}



/*
 *  walkup.c
 */
 
/*
  What:
  *****
  
  Implementing:
 
  The walkup routine within the VES structure:

  Walking up from w where [w, v^c] is a (directed)
  back edge to be embeeding later.
  Along the way collect all the pertinent bicomps that
  will need to be merged before embedding the back edges
  to v^c.


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/


#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
 
 
 
/* aproto: header embed_graph_protos.h */
 
 
#ifndef PLANAR_IN_MAGMA
#endif


void 
embedg_walkup (t_ver_edge *embed_graph, int n, int v, t_dlcl *p)
    /*
      walkup from w = p->info to v: [w, v] is a back edge where w is a DFS
      descendant of v
    */
{
    int          w, x, xin, y, yin;

    w = p->info;
    
    IF_DEB(
           fprintf(stdout, "walkup from %d to %d, enter\n", w, v);
           )
        
    embed_graph[w].adjacent_to = v;
    /*
      dirty trick to record some information about the BE [w, v]
      which will be useful at the time of creation and insertion of
      this BE: this happens in the walkdown procedure

      note that what I am doing here is safe: [w].in_adjl,
      [w].twin_in_adjl, [w].mult had no use so far since w is a vertex
      (and not an edge...)
    */
    embed_graph[w].in_adjl = p->in_adjl;
    embed_graph[w].twin_in_adjl = p->twin_in_adjl;
    embed_graph[w].mult = p->mult;
    
    /*
      set up the traversal contexts for w: one in each direction
    */
    x = w;
    xin = 1;
    y = w;
    yin = 0;

    while (x != v)
    {
        int          vz, z, c;

        IF_DEB(
               fprintf(stdout, "walkup, x %d and y %d\n", x, y);
               )
            
        if (embed_graph[x].visited == v
            || embed_graph[y].visited == v)
        {
            IF_DEB(
                   if (embed_graph[x].visited == v)
                       fprintf(stdout, "walkup, x visited\n");
                   else
                       fprintf(stdout, "walkup, y visited\n");
                   )
            break;
        }

        /*
          set x and y as visited!
        */
        embed_graph[x].visited = embed_graph[y].visited = v;
        
        vz = embedg_VES_is_virtual_vertex(n, x) ? x : NIL;
        vz = embedg_VES_is_virtual_vertex(n, y) ? y : vz;

        if (vz != NIL)
            /*
              that is, x (or y) is a virtual vertex
              -- in other words, we are set to find the root of the bicomp
              containing w, or of the bicomp r^c such that w is in the tree
              rooted by c

              consequently, by definition, vz is PERTINENT
            */
        {
            c = vz - n;
            z = embed_graph[c].DFS_parent;

            IF_DEB(
                   fprintf(stdout, "walkup, vz is virtual, %d^%d\n",
                           z, c);
                   )

            if (z != v)
                /*
                  determine if vz externally or internally active
                */
            {
                if (embed_graph[c].lowpoint < v)
                    /*
                      vz is externally active: APPEND to the list
                      of pertinent bicomps
                    */
                {
                    IF_DEB(
                           fprintf(stdout, "walkup, vz is ext. active\n");
                           )

                    embed_graph[z].pertinent_bicomp_list =
                        embedg_dlcl_rec_append(
                                embed_graph[z].pertinent_bicomp_list,
                                embedg_dlcl_rec_new(vz));
                }
                else
                    /*
                      vz is internally active: PREPEND to the list
                      of pertinent bicomps
                    */
                {
                    IF_DEB(
                           fprintf(stdout, "walkup, vz is pertinent\n");
                           )

                    embed_graph[z].pertinent_bicomp_list =
                        embedg_dlcl_rec_prepend(
                                embed_graph[z].pertinent_bicomp_list,
                                embedg_dlcl_rec_new(vz));
                }
            }
            
            /*
              continue the walkup, look if there are any other
              pertinent bicomps
              -- here "jump" to the next bicomp "up"
            */
            x = z;
            xin = 1;
            y = z;
            yin = 0;
        }
        else
            /*
              continue the traversal of the bicomp until one finds
              its (virtual) root
            */
        {
            embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                                 x, xin, FALSE, 0, &x, &xin);
            embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                                 y, yin, FALSE, 0, &y, &yin);
        }
    }
}

/*
 *  walkdown.c
 */
 
/*
  What:
  *****
  
  Implementing:

  The walkdown routine within the VES structure:

  walking down a bicomp rooted by a virtual vertex v^c
  and attempting to embed the back edges.
  This cannot be done if the walk has to stop due to the
  presence of externally active vertices on both
  the clockwise and the anticlockwise side of the bicomp.


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_DEB_EMBED(x)    {}
#define IF_DEB_BE(x) {}
#define IF_DEB_SCE(x) {}
#define IF_VERB(x)   {}
 
 
 
/* aproto: header embed_graph_protos.h */
 
 
#ifndef PLANAR_IN_MAGMA
#endif
 


t_merge_queue 
embedg_walkdown (t_ver_edge *embed_graph, int n, int *edge_pos, int vv)
    /*
      walkdown from the virtual vertex:
      embed any back edges incident to vv if any
      and merge the encountered bicomps while walking down
      (very informative isn't it? :))

      ... and return the merge queue: will be useful when
      isolating the Kuratowski subgraphs
    */
{
    t_merge_queue    q;
    int              v, c, vvout;
    
    ASSERT(embedg_VES_is_virtual_vertex(n, vv));

    /*
      find v and c such that v^c = vv
    */
    c = vv - n;
    v = embed_graph[c].DFS_parent;

    IF_DEB(
           fprintf(stdout, "walkdown from %d^%d, enter\n", v, c);
           )

    IF_DEB_EMBED(
                 fprintf(stdout, "walkdown, embedding at start\n");
                 embedg_VES_print_bigcomps(embed_graph, n);
                 )
        
    /*
      create an empty merge queue
    */
    q = embedg_merge_queue_new(n);

    for (vvout = 0; vvout <= 1; vvout++)
        /*
          chose a direction for the walk, but walk in both
          directions unless a stopping vertex is encountered
          and other conditions are satisfied (see below)
        */
    {
        int    w, win;

        embedg_VES_get_succ_on_ext_face(embed_graph, n, vv, vvout ^ 1,
                                             FALSE, 0, &w, &win);

        IF_DEB(
               fprintf(stdout, "walkdown, successor (outside while loop) from %d^%d:%d is %d:%d\n",
                       embed_graph[vv-n].DFS_parent, vv-n, vvout ^ 1,
                       w, win);
               )

        while (w != vv)
            /*
              is there no danger we walk the whole way back to vv
              and that all the vertices along the walk are inactive?

              answer: no, because of the short-cut edges.

              Short-cut edges are precisely inserted to remove the inactive
              vertices from the external face (ie they are "pushed"
              to the internal face of the bicomp)
            */
        {
            if (embed_graph[w].adjacent_to == v)
                /*
                  ie there is a (directed) back edge [w, v]
                  (would have been set in the previous walkup routine):
                  embed this edge, but before that, merge all the bicomps
                  previouslsy collected 
                */
            {
                IF_DEB(
                       fprintf(stdout, "walkdown, embed BE (%d^%d:%d, %d:%d)\n",
                               embed_graph[vv-n].DFS_parent, vv - n, vvout,
                               w, win);
                       fprintf(stdout, "walkdown, queue before pulling elts\n");
                       embedg_merge_queue_print(q);
                       )
            
                while (!embedg_merge_queue_empty(q))
                {
                    int     u, uin, vu, vuout;

                    embedg_merge_queue_get(&q, &u, &uin, &vu, &vuout);

                    IF_DEB(
                           fprintf(stdout, "walkdown, pull from queue (%d:%d, %d^%d:%d)\n",
                                   u, uin,
                                   embed_graph[vu-n].DFS_parent, vu-n,
                                   vuout);
                       )

                    embedg_VES_merge_pertinent_bicomps(
                                                        embed_graph, n,
                                                        vu, vuout, u, uin);
                }
                IF_DEB_BE(
                          fprintf(stdout, "walkdown, before embed BE [%d^%d:%d, %d:%d]\n",
                                  embed_graph[vv-n].DFS_parent, vv - n,
                                  vvout, w, win);
                          embedg_VES_print_adj_list(
                                                    embed_graph, n, vv,
                                                    TRUE);
                          fprintf(stdout, "\n");
                          embedg_VES_print_adj_list(
                                                    embed_graph, n, vv,
                                                    FALSE);
                          )
                    
                embedg_VES_embed_edge(embed_graph, n, edge_pos,
                                                 BE, vv, vvout, w, win);

                IF_DEB_BE(
                          fprintf(stdout, "walkdown, after embed BE [%d^%d:%d, %d:%d]\n",
                                  embed_graph[vv-n].DFS_parent, vv - n,
                                  vvout, w, win);
                          embedg_VES_print_adj_list(
                                                    embed_graph, n, vv,
                                                    TRUE);
                          fprintf(stdout, "\n");
                          embedg_VES_print_adj_list(
                                                    embed_graph, n, vv,
                                                    FALSE);
                          )
                IF_DEB_EMBED(
                       fprintf(stdout, "walkdown, embedding after bicomp merge & back edge embedding\n");
                       embedg_VES_print_bigcomps(embed_graph, n);
                       )
                        
                /*
                  clear the adjacent_to flag
                */
                embed_graph[w].adjacent_to = n;  /* "invalid" value */
            }

            if (!embedg_dlcl_is_empty(embed_graph[w].pertinent_bicomp_list))
                /*
                  each pertinent child bicomp of w
                  (pertinent: contains active (ie more back edges to embed)
                  elts)
                  must be traversed
                  and pushed onto the queue for later bicomp merging
                */
            {
                int           vw, vwout, x, xin, y, yin, s, sin;

                IF_DEB(
                       fprintf(stdout, "walkdown, pertinent list for %d\n",
                               w);
                       embedg_dlcl_print(embed_graph[w].pertinent_bicomp_list);
                       )
                
                /*
                  get the first child in the pertinent list
                  (see how the list is built in embedg_walkup)

                  the child will eventually be removed from that list
                  when merging the bicomps, and surely
                  this bicomp (rooted at vw) will be merged (later)
                  because it is active and hence pushed on
                  the merge queue
                */

                /*
                  we can start by pushing the vertex (w, win) on
                  the merge queue
                */
                embedg_merge_queue_append_vertex(&q, embed_graph, n, w, win);

                IF_DEB(
                       fprintf(stdout, "walkdown, push 1rst 2-tuple on queue\n");
                       embedg_merge_queue_print(q);
                       )
                    
                /*
                  get the first child in the pertinent list
                */
                vw = (embed_graph[w].pertinent_bicomp_list)->info;

                IF_DEB(
                       fprintf(stdout, "walkdown, get pertinent %d^%d\n",
                               embed_graph[vw - n].DFS_parent, vw - n);
                       )
                    
                /*
                  start two walks starting at vw
                */
                embedg_VES_get_succ_active_on_ext_face(embed_graph, n,
                                                     v , vw, 1,
                                                     FALSE, 0, &x, &xin);
                embedg_VES_get_succ_active_on_ext_face(embed_graph, n,
                                                     v, vw, 0,
                                                     FALSE, 0, &y, &yin);

                /*
                  because of the trick of inserting short-cut edges
                  at previous stages, neighbours of vw are guaranteed
                  to be active
                  
                  (however I'll use the more general 
                  embedg_VES_get_succ_active_on_ext_face
                  instead of the restrictive 
                  embedg_VES_get_succ_on_ext_face
                  because  the walkdown may be used later to isolate
                  Kuratowski minors, in a situation where SCEs could have
                  been removed and thus where the successor on the
                  external face will no longer be guaranteed to be active)
                  (* actually I have decided to remove the SCE at the
                  very last moment hence the above pb
                  does not occur in the present implementation)

                  
                  it only remains to chose the next vertex where from
                  to continue the walk; the choice is made in that order:
                  - an internally active vertex
                    (incident to v via a backedge but whose lowpoint 
                    is NO less than v)
                  - a (externally active) pertinent vertex
                    (incident to v via a backedge but whose lowpoint
                    is less than v: ie which is also externally active)
                  - as a last resort, a non-pertinent externally vertex,
                    which is then a stopping vertex
                */
                IF_DEB(
                       fprintf(stdout, "walkdown, x and y: %d, %d\n", x, y);
                       )
                
                if (embedg_VES_is_ver_int_active(embed_graph, n,
                                                            v, x))
                    /*
                      x is internally active
                    */
                {
                    IF_DEB(
                           fprintf(stdout, "walkdown, x is int. active\n");
                           )
                        
                    s = x;
                    sin = xin;
                }
                else if (embedg_VES_is_ver_int_active(
                                                            embed_graph, n,
                                                            v, y))
                    /*
                      y is internally active
                    */
                {
                    IF_DEB(
                           fprintf(stdout, "walkdown, y is int. active\n");
                           )
                        
                    s = y;
                    sin = yin;
                }
                else if (embedg_VES_is_ver_pertinent(
                                                            embed_graph, n,
                                                            v, x))
                    /*
                      x is pertinent
                    */
                {
                    IF_DEB(
                           fprintf(stdout, "walkdown, x is pertinent\n");
                           )
                        
                    s = x;
                    sin = xin;
                }
                else
                    /*
                      tough luck: y may be externally active
                    */
                {
                    IF_DEB(
                           fprintf(stdout, "walkdown, tough luck\n");
                           )
                        
                    s = y;
                    sin = yin;
                }

                IF_DEB(
                       fprintf(stdout, "walkdown, succ. on pertinent bicomp is %d:%d\n", s, sin);
                       )
                        
                /*
                  set vwout to respect consistency of traversal
                */
                vwout = s == x ? 0 : 1;

                /*
                  now that we know vwout we can push (vw, vwout)
                  on the merge queue, thus completing the 4-tuple
                  (w, win, vw, vwout) describing a bicomp merge
                  to occur at a later stage
                */
                embedg_merge_queue_append_virtual_vertex(&q, embed_graph, n,
                                                         vw, vwout);

                IF_DEB(
                       fprintf(stdout, "walkdown, push on queue (%d:%d, %d^%d:%d)\n",
                               w, win, embed_graph[vw-n].DFS_parent, vw - n,
                               vwout);
                       embedg_merge_queue_print(q);
                       )
                    
                /*
                  we continue the walk
                */
                w = s;
                win = sin;
            }
            /*
              at this point, w is either inactive or externally active
              (w can't be pertinent: its pertinent bicomp list is empty,
              and the back edge [w, v], if any, has already been embedded)
            */
            else if (embedg_VES_is_ver_inactive(embed_graph, n,
                                                           v, w))
                /*
                  w is inactive: continue with the walk on the external face
                  and, insert a short cut edge so that w is removed
                  from the external face
                */
            {
                int   s, sin;

                IF_DEB(
                       fprintf(stdout, "walkdown, %d has no pertinent bicomps and is inactive\n", w);
                       )
                    
                embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                                     w, win,
                                                     FALSE, 0, &s, &sin);

                IF_DEB(
                       fprintf(stdout, "walkdown, successor from %d:%d is %d:%d\n",
                               w, win, s, sin);
                       )

                /*
                  s is the successor of w: we embed a short circuit edge
                  [vv, s] if
                  - the bicomp is externally active (to ensure that
                    at a later stage this new face gets bisected:
                    so that we don't end up with a face of degree 2
                    (parallel edges))
                  - if [s, vv] is not a back edge

                  CONSEQUENTLY, adding SCE edges
                  + does not destroy the planarity of the graph
                  + ensures that each face has degree > 2 so that
                    |E| <= 3 * |V| - 6 remains valid at all times
                  + that the space allocated to the edges in embed_graph
                    (via MAXDE(n)) is sufficient

                  NOTE:
                  the above still allows to embed a short-cut edge
                  as an edge parallel to a tree edge OR a back edge
                  (which then has been embedded previously
                  so that [w].adjacent has been cleared)
                  
                  but again, since the degree of the face will be
                  > 2, that's ok 

                  recall that c = vv - n
                */
                if (embed_graph[c].lowpoint < v
                    /*
                      bicomp rooted at vv is externally active
                    */
                    && embed_graph[s].adjacent_to != v)
                    /*
                      [s, vv] is not a back edge
                    */
                {
                    IF_DEB_SCE(
                           fprintf(stdout, "walkdown, before embed SCE [%d^%d:%d, %d:%d]\n",
                                   embed_graph[vv-n].DFS_parent, vv - n,
                                   vvout, s, sin);
                           embedg_VES_print_adj_list(
                                                    embed_graph, n, vv,
                                                    TRUE);
                           fprintf(stdout, "\n");
                           embedg_VES_print_adj_list(
                                                    embed_graph, n, vv,
                                                    FALSE);
                           
                           )

                    embedg_VES_embed_edge(embed_graph,
                                                     n, edge_pos,
                                                     SCE, vv, vvout, s, sin);
                    /*
                      note also that the addition of short cut edges
                      does not change the fact that the graph is planar
                      (when it is, so we never run into the problem
                      of creating/adding too many edges to embed-graph)
                    */
                    IF_DEB_SCE(
                           fprintf(stdout, "walkdown, after embed SCE [%d^%d:%d, %d:%d]\n",
                                   embed_graph[vv-n].DFS_parent, vv - n,
                                   vvout, s, sin);
                           embedg_VES_print_adj_list(
                                                    embed_graph, n, vv,
                                                    TRUE);
                           fprintf(stdout, "\n");
                           embedg_VES_print_adj_list(
                                                    embed_graph, n, vv,
                                                    FALSE);
                           
                           )
                    IF_DEB(
                           fprintf(stdout, "walkdown, embed SCE [%d^%d:%d, %d:%d]\n",
                                   embed_graph[vv-n].DFS_parent, vv - n,
                                   vvout, s, sin);
                           )

                }
                /*
                  continue the walk
                */
                w = s;
                win = sin;
            }
            else
                /*
                  w is non-pertinent and externally active:
                  it is a stopping vertex:
                  we stop here and see if we can walk in the other direction
                */
            {
                IF_DEB(
                       fprintf(stdout, "walkdown, %d is externally active\n", w);
                       )
                break;
            }
        }
        if (!embedg_merge_queue_empty(q))
            /*
              mumm.... don't understand this one... let's see:
              the queue constains pertinent bicomps collected during one of
              the traversal of the external face, so that once
              a stopping vertex has been encountered and the queue
              is not empty, this means that we will be unable
              to embed any remaining back edges:

              it is important to remember that when w is a stopping vertex
              there is no choice left, since we walk the pertinent
              bicomp in both directions at once, and always choose
              the "best" possible vertex
              (see the choice strategy: (a) internally active, (b) pertinent,
              (c) the rest)
            */
        {
            IF_DEB(
                   fprintf(stdout, "walkdown, merge queue is not empty\n");
                   )
            break;
        }
    }

    /*
      and return the merge queue
    */
    return q;
}




/*
 *  merge_queue_misc.c
 */
 
/*
  What:
  *****
  
  Implementing:

  The merge queue stores the pertinent bicomps waiting to
  be merged before a subsequent back edge embedding.
  See walkdown.c

  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
 
 
 
/* aproto: header embed_graph_protos.h */

 
#ifndef PLANAR_IN_MAGMA
#endif

t_merge_queue 
embedg_merge_queue_new (int n)
    /*
      create a merge queue of 4 * (n-1) elts:
      we can only have at most n-1 virtual vertices,
      and for each of those we need to store 4 bits of info
    */
{
    t_merge_queue   q;

    q.start = q.end = 0;
    q.b = (int *) mem_malloc(sizeof(int) * 4 * (n - 1));

    return q;
}

void 
embedg_merge_queue_delete (t_merge_queue q)
{
    mem_free(q.b);
}


boolean 
embedg_merge_queue_empty (t_merge_queue q)
{
    return q.start == q.end ? TRUE : FALSE;
}

void 
embedg_merge_queue_print (t_merge_queue q)
{
    int        i;

    for (i = q.start; i < q.end; i++)
    {
        fprintf(stdout, "%d:%d ", q.b[i], q.b[i+1]);
	++i;
    }
    fprintf(stdout, "\n");
}

void 
embedg_merge_queue_append (t_merge_queue *q, t_ver_edge *embed_graph,
	int n, int v, int vin, int vv, int vvout)
    /*
      append the 4-tuple (v, vin, vv, vvout)
      where v is a vertex and vv is its virtual counterpart

      we don't do much here, most of the work is done
      when pulling a bicomp/4-tuple from the queue
    */
{
    /*
      is this really necessary? - YES!!!
    */
    ASSERT((*q).end < 4 * (n - 2));
    ASSERT(embedg_VES_is_vertex(n, v));
    ASSERT(embedg_VES_is_virtual_vertex(n, vv));
    ASSERT(embed_graph[vv - n].DFS_parent == v);

    (*q).b[(*q).end++] = v;
    (*q).b[(*q).end++] = vin;
    (*q).b[(*q).end++] = vv;
    (*q).b[(*q).end++] = vvout;
}

void 
embedg_merge_queue_append_vertex (t_merge_queue *q, t_ver_edge *embed_graph,
	int n, int v, int vin)
    /*
      same as above but were we only append the 2-tuple (v, vin),
      appending the 2-tuple (vv, vvout) at a later stage
      (see embedg_merge_queue_append_virtual_vertex)
    */
{
    ASSERT((*q).end < 4 * (n - 2));
    ASSERT(embedg_VES_is_vertex(n, v));

    (*q).b[(*q).end++] = v;
    (*q).b[(*q).end++] = vin;

    IF_DEB(
           fprintf(stdout, "merge_queue_append_vertex, after, end is %d\n",
                   (*q).end);
           )
}

void 
embedg_merge_queue_append_virtual_vertex (t_merge_queue *q,
	t_ver_edge *embed_graph, int n, int vv, int vvout)
    /*
      counterpart to embedg_merge_queue_append_vertex:
      here we append the 2-tuple (vv, vvout), vv = v^c,
      where the 2-tuple (v, vin) is already in the queue
      (see embedg_merge_queue_append_vertex)
    */
{
    ASSERT(!embedg_merge_queue_empty(*q));
    ASSERT(embedg_VES_is_virtual_vertex(n, vv));
    ASSERT(embed_graph[vv - n].DFS_parent == (*q).b[(*q).end - 2]);
    
    (*q).b[(*q).end++] = vv;
    (*q).b[(*q).end++] = vvout;

    IF_DEB(
           fprintf(stdout, "merge_queue_append_virtual_vertex, after, end is %d\n",
                   (*q).end);
           )
}

void 
embedg_merge_queue_get (t_merge_queue *q, int *v, int *vin, int *vv, int *vvout)
    /*
      pulling out a 4-tuple from the beginning of the FIFO queue
    */
{
    ASSERT(!embedg_merge_queue_empty((*q)));

    *v = (*q).b[(*q).start++];
    *vin = (*q).b[(*q).start++];
    *vv = (*q).b[(*q).start++];
    *vvout = (*q).b[(*q).start++];
}

void 
embedg_merge_queue_prune (t_merge_queue *q, int *v,
			  int *vin, int *vv, int *vvout)
    /*
      pulling out a 4-tuple from the end of the FIFO queue
    */
{
    ASSERT(!embedg_merge_queue_empty((*q)));

    *vvout = (*q).b[--((*q).end)];
    *vv = (*q).b[--((*q).end)];
    *vin = (*q).b[--((*q).end)];
    *v = (*q).b[--((*q).end)];
}

/*
 *  vertex_activity.c
 */
 
/*
  What:
  *****
  
  Implementing:

  Determining a vertex's activity. This takes place within
  the VES structure.


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/
 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
 
 
 
/* aproto: header embed_graph_protos.h */
 
 
#ifndef PLANAR_IN_MAGMA
#endif
 
 


boolean 
embedg_VES_is_ver_pertinent (t_ver_edge *embed_graph, int n, int v, int w)
    /*
      is w pertinent (wrt v)
      - the field adjacent_to = v: means there is a back edge [w, v]
      - or w has a non empty pertinent_bicomp_list
    */
{
    boolean       ans;
    
    ans = embed_graph[w].adjacent_to == v ? TRUE : FALSE;

    if (ans)
        return TRUE;
    else
        return embedg_dlcl_is_empty(embed_graph[w].pertinent_bicomp_list) ?
            FALSE : TRUE;
}

boolean 
embedg_VES_is_ver_ext_active (t_ver_edge *embed_graph, int n, int v, int w)
    /*
      is w externally active (wrt v)
      this is the case when either w's least_ancestor < v
      or the first member of w's separated_DFS_child_list has lowpoint < v
      (the vertices in separated_DFS_child_list are ordered by lowpoint)

      why? because w's separated_DFS_child_list may be empty
      (due to prior bicomp merging say) and so its children are in effect
      inactive
    */
{
    boolean       ans;

    ans = embed_graph[w].least_ancestor < v ? TRUE : FALSE;

    if (ans)
        return TRUE;
    else
    {
        if (embedg_dlcl_is_empty(embed_graph[w].separated_DFS_child_list))
        {
            return FALSE;
        }
        else
        {
            int      c;
            
            c = (embed_graph[w].separated_DFS_child_list)->info;
            return embed_graph[c].lowpoint < v ? TRUE : FALSE;
        }
    }
}


boolean 
embedg_VES_is_ver_int_active (t_ver_edge *embed_graph, int n, int v, int w)
    /*
      is w internally active (wrt v):
      this happens when w is pertinent but NOT externally active
    */
{
    return embedg_VES_is_ver_pertinent(embed_graph, n, v, w)
        && !embedg_VES_is_ver_ext_active(embed_graph, n, v, w);
}

boolean 
embedg_VES_is_ver_inactive (t_ver_edge *embed_graph, int n, int v, int w)
    /*
      is w inactive (wrt v), that is w nor pertinent nor externally activ
    */
{
    return !embedg_VES_is_ver_pertinent(embed_graph, n, v, w)
        && !embedg_VES_is_ver_ext_active(embed_graph, n, v, w);
}

/*
 *  merge_bicomps.c
 */
 
/*
  What:
  *****
  
  Implementing:

  In the VES structure, merging two bicomponents.
  That is, merging the virtual vertex v^c with the
  actual vertex v while merging their respective
  adjacency lists.
  This must be done in a very specific manner so as to able to
  determine the subsequent internal/external faces.
  Also, great care must be taken so that the resulting
  adj. list for v is consistent (wrt to the direction
  of traversal).


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_DEB_ADJL(x)    {}
#define IF_VERB(x)   {}

 
 
/* aproto: header embed_graph_protos.h */


void 
embedg_VES_merge_simple_bicomps (t_ver_edge *embed_graph, int n, int vv,
	int vvout, int v, int vin)
    /*
      merge the bicomp rooted at vv (vv a virtual vertex) with
      its counterpart v so that the resulting adjacency list for v
      is consistent and is the union of the adjacency lists for vv and v

      we treat the case that the bicomp may be flipped (vvout == vin)
      here
    */
{
    int       c, edge, twin, root_edge, cur, prev;
    int       vout, vvin, e1, e2, e3, e4, e1out, e3out, e4in;

    /*
      find c such that [v^c, c] is the root edge of the bicomp
      rooted at vv = v^c
    */
    c = vv - n;
    ASSERT(embed_graph[c].DFS_parent == v);

    IF_DEB(
           fprintf(stdout, "merge_simple_bicomp, start: merge\n");
           embedg_VES_print_virtual_vertex(embed_graph, n, vv);
           fprintf(stdout, ":%d & ", vvout);
           embedg_VES_print_vertex(n, v);
           fprintf(stdout, ":%d\n", vin);
           )

    IF_DEB_ADJL(
           fprintf(stdout, "merge_simple_bicomp, adj. list for %d (before)\n", vv);
           embedg_VES_print_adj_list(embed_graph, n, vv,
                                                TRUE);
           fprintf(stdout, "\n");
           embedg_VES_print_adj_list(embed_graph, n, vv,
                                                FALSE);
           fprintf(stdout, "\n");

           fprintf(stdout, "merge_simple_bicomp, adj. list for %d (before)\n", v);
           embedg_VES_print_adj_list(embed_graph, n, v,
                                                TRUE);
           fprintf(stdout, "\n");
           embedg_VES_print_adj_list(embed_graph, n, v,
                                                FALSE);
           )
    /*
      find all edges incident to vv and (re)set all references
      to incidence to vv to incidence to v

      by the same token, find the root_edge [v^c, c]

      MOREVOVER, when vin == vvout, the bicomp (rooted by v^v = vv)
      will be flipped:
      we must invert the links of all the edges incident
      to vv so that their further union with v's adjacency list
      results in a consistent adjacency list for v!

      we do everything in one go
    */

    /*
      very careful here: a root edge must ALSO be a TE
      (because the same edge could have been added as a SCE)
    */
        
    root_edge = NIL;
    edge = embed_graph[vv].link[vvout];
    ASSERT(embedg_VES_is_edge(n, edge));
    if (embed_graph[edge].neighbour == c
        && embedg_VES_is_tree_edge(embed_graph, n, edge))
    {
        root_edge = edge;
    }

    if (vin == vvout)
        /*
          invert the links
        */
    {
        int  in, out;

        in = embed_graph[edge].link[0];
        out = embed_graph[edge].link[1];
        embed_graph[edge].link[0] = out;
        embed_graph[edge].link[1] = in;
    }
    /*
      get the twin and set the neighbour there to v (was vv originally)
    */
    twin =  embedg_VES_get_twin_edge(embed_graph, n, edge);
    ASSERT(embed_graph[twin].neighbour == vv);
    embed_graph[twin].neighbour = v;

    prev = vv;
    cur = edge;
    while (edge != vv)
    {
        edge =
            embedg_VES_get_next_in_dlcl(embed_graph, n,
                                                   cur, prev);

        if (embedg_VES_is_edge(n, edge))
            /*
              get the twin again (and invert the links if need be)
            */
        {
            if (embed_graph[edge].neighbour == c
                && embedg_VES_is_tree_edge(embed_graph, n, edge))
            {
                root_edge = edge;
            }

            if (vin == vvout)
            {
                int  in, out;
                
                in = embed_graph[edge].link[0];
                out = embed_graph[edge].link[1];
                embed_graph[edge].link[0] = out;
                embed_graph[edge].link[1] = in;
            }
            
            twin =
                embedg_VES_get_twin_edge(embed_graph, n, edge);
            ASSERT(embed_graph[twin].neighbour == vv);
            embed_graph[twin].neighbour = v;

            prev = cur;
            cur = edge;
        }
        else
        {
            ASSERT(edge == vv);
            /*
              only one vertex in the whole circular list
            */
        }
    }
    ASSERT(root_edge != NIL);

    /*
      and now union the adjacency lists of v and vv:

      let e1 be the edge record used to enter v
          e2                            exit  v
          e3                            enter vv 
          e4                            exit  vv :

          e1 -> v  -> e2
          e3 -> vv -> e4

      the union of the list is done in such a way that
      - e1 and e4 are consecutive in v's adjacency list:
        they are now in the internal face
      - e3 is now the edge record used to enter v:
        it is on the external face (along with e2) :

          e1 -> e4
          e3 -> v -> e2

      (note that this does not assume that e1 & e2 are distinct
      or that e3 & e4 are distinct)
    */
    /*
      I must not forget the case where v is a lone vertex:
      this is the case where v has no DFS ancestor, ie when
      v is the root of a tree in the DFS forest
    */

    e1 = embed_graph[v].link[vin];
    vout = 1 ^ vin;
    e2 = embed_graph[v].link[vout];

    if (e1 != v)
    {
        ASSERT(e2 != v);
        ASSERT(embedg_VES_is_edge(n, e1));
        ASSERT(embedg_VES_is_edge(n, e2));
    }
    
    e4 = embed_graph[vv].link[vvout];
    ASSERT(embedg_VES_is_edge(n, e4));

    vvin = 1 ^ vvout;
    e3 = embed_graph[vv].link[vvin];
    ASSERT(embedg_VES_is_edge(n, e3));
    
    /*
      must take care of the adjacency list's consistency of traversal
      (will be important only when recovering the embedding)
    */
    if (e1 == e2)
    {
        ASSERT(embed_graph[e1].link[0] == embed_graph[e1].link[1]);
        if (vin == vvout)
            /*
              the bicomp will be flipped:
              must take 1 ^ vvout - difficult to explain -- later...
            */
        {
            e1out = 1 ^ vvout;
        }
        else
        {
            e1out = vvout;
        }
    }
    else
    {
        e1out = embed_graph[e1].link[0] == v ? 0 : 1;
    }
    if (e3 == e4)
    {
        ASSERT(embed_graph[e3].link[0] == embed_graph[e3].link[1]);
        e3out = 1 ^ vin;
        e4in = vin;
    }
    else
    {
        e4in = embed_graph[e4].link[0] == vv ? 0 : 1;
        e3out = embed_graph[e3].link[0] == vv ? 0 : 1;
    }

    IF_DEB(
           fprintf(stdout, "merge_simple_bicomp, before union of lists, e1\n");
           embedg_VES_print_edge(embed_graph, n, e1);
           fprintf(stdout, "merge_simple_bicomp, e3\n");
           embedg_VES_print_edge(embed_graph, n, e3);
           fprintf(stdout, "merge_simple_bicomp, e4\n");
           embedg_VES_print_edge(embed_graph, n, e4);
           )

    /*
      make e1 and e4 consecutive in the adjacency list
    */
    embed_graph[e1].link[e1out] = e4;
    embed_graph[e4].link[e4in] = e1;
    embed_graph[e3].link[e3out] = v;
    embed_graph[v].link[vin] = e3;

    IF_DEB(
           fprintf(stdout, "merge_simple_bicomp, after union of lists, e1\n");
           embedg_VES_print_edge(embed_graph, n, e1);
           fprintf(stdout, "merge_simple_bicomp, e3\n");
           embedg_VES_print_edge(embed_graph, n, e3);
           fprintf(stdout, "merge_simple_bicomp, e4\n");
           embedg_VES_print_edge(embed_graph, n, e4);
           )

    /*
      also, want to "disable" vv links, meaning then that
      vv is no longer a root of a bicomp
    */
    embed_graph[vv].link[0] = embed_graph[vv].link[1] = vv;

    IF_DEB_ADJL(
           fprintf(stdout, "merge_simple_bicomp, adj. list for %d (after)\n", vv);
           embedg_VES_print_adj_list(embed_graph, n, vv,
                                                TRUE);
           fprintf(stdout, "\n");
           embedg_VES_print_adj_list(embed_graph, n, vv,
                                                FALSE);
           fprintf(stdout, "\n");

           fprintf(stdout, "merge_simple_bicomp, adj. list for %d (after)\n", v);
           embedg_VES_print_adj_list(embed_graph, n, v,
                                                TRUE);
           fprintf(stdout, "\n");
           embedg_VES_print_adj_list(embed_graph, n, v,
                                                FALSE);
           )

    ASSERT(embedg_VES_is_adj_list_consistent(embed_graph, n, v));

    /*
      finally, give an orientation to the (formerly) root edge [vv, c]
      to keep traversal consistent (when recovering embedding)
    */
    if (vin == vvout)
        /*
          flip: set the sign of the root edge to clockwise

          note: a bicomp is merged only once, so there is no need to
          "flip" the root_edge's sign: it is set once at initialisation
          and then changed here if need be.
        */
    {
        embed_graph[root_edge].sign = CLOCKW;

        IF_VERB(
            fprintf(stdout, "merge_simple_bicomp, flip for %d, sign is now %d for %d of type %d\n",
                    c, embed_graph[root_edge].sign, root_edge, embed_graph[root_edge].type);
            embedg_VES_print_edge(embed_graph, n, root_edge);
            )
    }
}



void 
embedg_VES_merge_pertinent_bicomps (t_ver_edge *embed_graph, int n,
	int vv, int vvout, int v, int vin)
    /*
      the bicomps to be merged are pertinent: on top (and before)
      performing a simple merge, there are several things to do
      related to the merging to pertinent bicomps
    */
{
    /*
      a note of caution:
      it is (very) likely that after a bicomp merge the resulting
      bicomp is not biconnected (and hence traversal of the external face
      of the bicomp via embedg_VES_get_succ_on_ext_face is non-sensical)

      remembering that a PERTINENT bicomp merge is ALWAYS followed
      by a back edge embedding we see that the end result is then a bicomp
      where again traversal of the external face
      via embedg_VES_get_succ_on_ext_face will make sense
    */
    t_dlcl    *pertinent_list, *head, *rep_in_parent_list, *parent_list;
    int       c;

    /*
      find c such that [v^c, c] is the root edge of the bicomp
      rooted at vv = v^c
    */
    c = vv - n;
    ASSERT(embed_graph[c].DFS_parent == v);

    /*
      two things to do first:
      - remove vv from head of pertinent_bicomp_list of v
      - remove c from separated_DFS_child_list of v

      one may ask the point of this since the separated_DFS_child_list
      seems to mirror pertinent_bicomp_list: but this is not exactly so:
      + pertinent_bicomp_list is ordered according to the activity
        of the (virtual) vertices
      + separated_DFS_child_list is ordered according to the vertices'
        lowpoint values
      in effect, it could (almost?*) be said that these two lists
      are in reverse order (the *almost bit would warrant some thinking here)
    */
    
    /*
      remove vv from head of pertinent_bicomp_list of v
    */
    pertinent_list = head = embed_graph[v].pertinent_bicomp_list;
    ASSERT(!embedg_dlcl_is_empty(pertinent_list));
    ASSERT(head->info == vv);

    IF_DEB(
           fprintf(stdout, "merge_pertinent_bicomp, start: merge\n");
           embedg_VES_print_virtual_vertex(embed_graph, n, vv);
           fprintf(stdout, ":%d & ", vvout);
           embedg_VES_print_vertex(n, v);
           fprintf(stdout, ":%d\n", vin);
           )

    IF_DEB(
           fprintf(stdout, "merge_pertinent_bicomp, pertinent bicomp_list of %d (before)\n", v);
           embedg_dlcl_print(embed_graph[v].pertinent_bicomp_list);
           )

    
    embed_graph[v].pertinent_bicomp_list =
        embedg_dlcl_delete_first(pertinent_list);

    IF_DEB(
           fprintf(stdout, "merge_pertinent_bicomp, pertinent bicomp_list of %d (after)\n", v);
           embedg_dlcl_print(embed_graph[v].pertinent_bicomp_list);
           )

    /*
      vv = v^c: remove c from separated_DFS_child_list of v
    */
    rep_in_parent_list = embed_graph[c].rep_in_parent_list;
    ASSERT(!embedg_dlcl_is_empty(rep_in_parent_list));

    parent_list = embed_graph[v].separated_DFS_child_list;
    ASSERT(!embedg_dlcl_is_empty(parent_list));
    embed_graph[v].separated_DFS_child_list = 
        embedg_dlcl_delete_rec(parent_list, rep_in_parent_list);

    /*
      that's it, it remains to merge, ie. union the adjacency list,
      and flipping the bicomp if necessary
    */
    embedg_VES_merge_simple_bicomps(embed_graph, n,
                                               vv, vvout, v, vin);
}




    
    
/*
 *  embed_edge.c
 */
 
/*
  What:
  *****
  
  Implementing:

  Embedding an edge so that it lies on the external face of a bicomp.
  We work here with the VES structure.


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_DEB_EMBED(x)    {}
#define IF_VERB(x)   {}
 
 
/* aproto: header embed_graph_protos.h */
 

void 
embedg_VES_embed_edge (t_ver_edge *embed_graph, int n, int *edge_pos,
	int edge_type, int vv, int vvout, int w, int win)
    /*
      embed the edge (vv, w) (vv a virtual vertex, w a vertex) between
      vv and the edge vvout
      and the edge win and w

      so that after the embedding, one exits vv via (vv, w) and
      enters w via the twin (w, vv)
    */
{
    int       temp, tempin, tempout;

    ASSERT(edge_type == BE || edge_type == SCE);
    ASSERT(embedg_VES_is_virtual_vertex(n, vv));
    ASSERT(embedg_VES_is_vertex(n, w));

    IF_DEB(
           fprintf(stdout, "embed_edge, (%d:%d)\n", vv, w);
           )
    
    /*
      first, set the edge [vv, w] with the appropriate info

      when [vv, w] is a back edge there is some more work to do
      (see the walkup procedure for the extra information we need
      to copy here
    */
    (*edge_pos)++;
    ASSERT(*edge_pos < 2*n + 2 * MAXE(n));
    embed_graph[*edge_pos].neighbour = w;
    embed_graph[*edge_pos].type = edge_type;
    embed_graph[*edge_pos].sign = CCLOCKW;
    if (edge_type == BE)
    {
        ASSERT(embed_graph[w].adjacent_to ==
               embed_graph[vv - n].DFS_parent);

        /*
          PLUS: originally when the back edge [w, vv] was
          created (in the dfs preprocessing stage), it carried in
          .in_adjl the index of this directed edge in the
          adjacency list

          but now, note that we are actually inserting the
          directed edge [vv, w] in vv's adjacency list,
          meaning that in_adjl and twin_in_adjl
          must be exchanged!
        */
        embed_graph[*edge_pos].in_adjl = embed_graph[w].twin_in_adjl;
        embed_graph[*edge_pos].twin_in_adjl = embed_graph[w].in_adjl;

        ASSERT(embed_graph[w].mult % 2 == 0);
        /*
          the original graph is always undirected:
          we store its number of undirected edges
        */
        embed_graph[*edge_pos].mult = embed_graph[w].mult / 2;
    }
        
    /*
      insert this edge between vertex record for vv
      and edge record vv.link[vvout]
    */
    temp = embed_graph[vv].link[vvout];
    
    if (embed_graph[temp].link[0] == embed_graph[temp].link[1])
        /*
          this needs special treatment to ensure consistency of
          orientation
        */
    {
        ASSERT(embed_graph[temp].link[0] == vv);
        tempin = 1 ^ vvout;
    }
    else
    {
        tempin = embed_graph[temp].link[0] == vv ? 0 : 1;
    }

    IF_DEB(
           fprintf(stdout, "embed_edge, edge out of vv\n");
           embedg_VES_print_edge(embed_graph, n, temp);
           )
    
    embed_graph[vv].link[vvout] = *edge_pos;
    embed_graph[temp].link[tempin] = *edge_pos;
    /*
      the links for *edge_pos must also be "consistent"
    */
    embed_graph[*edge_pos].link[vvout] = temp;
    embed_graph[*edge_pos].link[vvout ^ 1] = vv;

    /*
      now create/set the twin edge, the directed edge [w, vv]
    */
    (*edge_pos)++;
    ASSERT(*edge_pos < 2*n + 2 * MAXE(n));
    embed_graph[*edge_pos].neighbour = vv;
    embed_graph[*edge_pos].type = edge_type;
    embed_graph[*edge_pos].sign = CCLOCKW;
    if (edge_type == BE)
    {
        embed_graph[*edge_pos].in_adjl = embed_graph[w].in_adjl;
        embed_graph[*edge_pos].twin_in_adjl = embed_graph[w].twin_in_adjl;
        embed_graph[*edge_pos].mult = embed_graph[w].mult / 2;
    }
 
    /*
      and insert the twin edge between edge record w.link[win]
      and vertex record for w
    */
    temp = embed_graph[w].link[win];

    if (embed_graph[temp].link[0] == embed_graph[temp].link[1])
        /*
          again, special treatment to ensure consistency of orientation
        */
    {
        ASSERT(embed_graph[temp].link[0] == w);
        tempout = 1 ^ win;
    }
    else
    {
        tempout = embed_graph[temp].link[0] == w ? 0 : 1;
    }
    
    IF_DEB(
           fprintf(stdout, "embed_edge, edge in of w\n");
           embedg_VES_print_edge(embed_graph, n, temp);
           )
    
    embed_graph[w].link[win] = *edge_pos;
    embed_graph[temp].link[tempout] = *edge_pos;
    /*
      and consistent orientation
    */
    embed_graph[*edge_pos].link[win] = temp;
    embed_graph[*edge_pos].link[win ^ 1] = w;
}



void 
embedg_VES_add_edge (t_ver_edge *embed_graph, int n, int *edge_pos,
	int v, int w, boolean MARK, int mark)
    /*
      add the edge (v, w): this is DIFFERENT from
      embedg_VES_embed_edge in the sense
      that the present function will only be used
      when building the Kuratowski homeomorphs:

      that is, we are in a situation where the graph is NON planar

      consequently it doesn't matter much where in the adjacency
      lists of v & w the edge is added:
      let's say that we always add it at the beginning

      for our sanity's sake, we'll ensure that the resulting
      adjacency lists remain consistent!
      
      and we add the edge as a BE!
      PLUS we mark it with mark in MARK true
    */
{
    int       temp;

    ASSERT(embedg_VES_is_vertex(n, v) ||
           embedg_VES_is_virtual_vertex(n, v));
    ASSERT(embedg_VES_is_vertex(n, w) ||
           embedg_VES_is_virtual_vertex(n, w));

    IF_DEB(
           fprintf(stdout, "add_edge, (%d:%d)\n", v, w);
           )

    /*
      not sure this is the best place to do this: mark the endpoints
    */
    if (MARK)
    {
        embed_graph[v].visited = mark;
        embed_graph[w].visited = mark;
    }
        
    /*
      first, set the edge [v, w] with the appropriate info
    */
    (*edge_pos)++;
    ASSERT(*edge_pos < 2*n + 2 * MAXE(n));
    embed_graph[*edge_pos].neighbour = w;
    embed_graph[*edge_pos].type = BE;
    /*
      the edge's orientation will be the same as the vertex
    */
    embed_graph[*edge_pos].sign = embed_graph[v].sign;
    /*
      and mark the edge
    */
    if (MARK)
    {
        embed_graph[*edge_pos].visited = mark;
    }
        
    /*
      insert this edge between vertex record for v 
      and edge record v.link[1]
    */
    temp = embed_graph[v].link[1];

    IF_DEB(
           fprintf(stdout, "add_edge, edge out of v\n");
           embedg_VES_print_edge(embed_graph, n, temp);
           )
    
    embed_graph[v].link[1] = *edge_pos;
    embed_graph[temp].link[0] = *edge_pos;
    /*
      the links for *edge_pos must also be "consistent"
    */
    embed_graph[*edge_pos].link[1] = temp;
    embed_graph[*edge_pos].link[0] = v;

    /*
      now create/set the twin edge, the directed edge [w, v]
    */
    (*edge_pos)++;
    ASSERT(*edge_pos < 2*n + 2 * MAXE(n));
    embed_graph[*edge_pos].neighbour = v;
    embed_graph[*edge_pos].type = BE;
    embed_graph[*edge_pos].sign = embed_graph[w].sign;
    if (MARK)
    {
        embed_graph[*edge_pos].visited = mark;
    }
 
    /*
      insert this edge between vertex record for w
      and edge record w.link[1]
    */
    temp = embed_graph[w].link[1];
    
    IF_DEB(
           fprintf(stdout, "add_edge, edge out of w\n");
           embedg_VES_print_edge(embed_graph, n, temp);
           )
    
    embed_graph[w].link[1] = *edge_pos;
    embed_graph[temp].link[0] = *edge_pos;
    /*
      and consistent orientation
    */
    embed_graph[*edge_pos].link[1] = temp;
    embed_graph[*edge_pos].link[0] = w;
}


/*
 *  recover.c
 */
 
/*
  What:
  *****
  
  Implementing:

  From the VES data structure recover either the embedding ot
  the obstruction into the

  t_sparseg_ver_struct,
  t_sparseg_adjl_struct,
  t_sparseg_embed_struct
  
  data types.


  (This is no even quite true: for some obscure reason
  I recover the obstruction as a dlcl[] structure to be
  converted later.
  The obvious reason being that it is easier to check as such.
  Maybe I leave it as it is...)

  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_DEB_EMBED_MULT(x) {}
#define IF_DEB_EMBED_LOOPS(x) {}
#define IF_DEB_EMBED(x)    {}
#define IF_DEB_CHECK_EMBED(x)    {}
#define IF_DEB_FACES(x) {}
#define IF_VERB(x)   {}
#define IF_DEB_SCE(x) {}
#define IF_DEB_OBS(x) {}
#define IF_DEB_CHECK_OBS(x) {}
#define IF_CPU(x) {}
 


/* aproto: header embed_graph_protos.h */

/* aproto: beginstatic -- don't touch this!! */
static void embedg_recover_embedding_embed_mult
 (t_dlcl **, t_embed_sparse_rep *, int, int, int, int, int *, boolean *, int *);
static void embedg_recover_embedding_embed_loops
 (t_dlcl **, t_embed_sparse_rep *, int, int, int *, boolean *);
static t_dlcl **embedg_get_reduced_obs (t_dlcl **, int);
static boolean embedg_is_red_obs_K33 (t_dlcl **, int);
static boolean embedg_is_red_obs_K5 (t_dlcl **, int);
/* aproto: endstatic -- don't touch this!! */
 
 
#ifndef PLANAR_IN_MAGMA
#endif

void 
embedg_recover_embedding (
    t_ver_sparse_rep *V,
    t_adjl_sparse_rep *A,             /* input (original sparse graph) */
    t_ver_edge *embed_graph,
    int n,
    int nbr_e,
    t_dlcl **mult_edges,
    t_ver_sparse_rep **vertices,
    t_embed_sparse_rep **embedding
)
    /*
      recover the embedding
      to prepare for the final Magma type for sparse & embedded graph

      we assume that all vertices/edges have been given their
      orientation
 
      at this stage we also embed the multiple edges and loops
      which were set aside in mult_edges by
      sparseg_adjl_dfs_preprocessing:

      as it turns out the last bit is pretty hairy!
    */
{
    /*
      the idea is to return an array of vertices and an array
      representing the embedding
      (careful: need to weedout the SCE)

      vertices: (*vertices)[i].first_edge contains index
                to first edge in embedding

      embedding: a doubly linked circular list of edges,
                 for each record/edge e = (*embedding)[i]:
                 e.in_adjl: index in A of e
                 e.next:    next edge in CLOCKW
                            (as an index in the embedding)
                 e.prev:    previous edge in CLOCKW
                            (as an index in embedding)
                 e.inv:     inverse edge (as an index in embedding)
                 e.mark:    a mark for this edge

      let's say that this new array is a slimmed down version of embed_graph

      one issue to address:
      - for edge e, find its index in A: this should be found
        in either the embed_graph[v] record of the mult_edges[v] record
    */
    int          index_embed, v, mult, w, v_w_in_embed, new_first_edge;
    boolean      set_next;

    IF_DEB(
           fprintf(stdout, "in recover emb.\n");
           sparseg_dlcl_print(mult_edges, n);
           );
    
    *vertices = (t_ver_sparse_rep *)
        mem_malloc(sizeof(t_ver_sparse_rep) * n);
    *embedding = (t_embed_sparse_rep *)
        mem_malloc(sizeof(t_embed_sparse_rep) * 2 * nbr_e);

    index_embed = 0;
    set_next = TRUE;
    for (v = 0; v < n; v++)
    {
        int       v_l, orient, in, out, e, cur_e, next_e;

        /*
          we take v's label
        */
        v_l = embed_graph[v].label;

        /*
          first let's deal with the isolated vertex case: those
          that refer to self
        */
        if (embed_graph[v].link[0] == v)
        {
            int       temp_index_embed;
            
            ASSERT(embed_graph[v].link[1] == v);

            /*
              there may be [v, v] loops for this vertex, must check this
            */
            temp_index_embed = index_embed - 1;
            /*
              temp_index_embed is pre-increased below
            */
            embedg_recover_embedding_embed_loops(mult_edges, *embedding,
                                                 nbr_e, v,
                                                 &temp_index_embed,
                                                 &set_next);
            
            if (temp_index_embed > index_embed - 1)
                /*
                  must fix beginning and end of adjacency list:
                */
            {
                (*vertices)[v_l].first_edge = index_embed;
                (*embedding)[temp_index_embed].next =
                    (*vertices)[v_l].first_edge;
                (*embedding)[(*vertices)[v_l].first_edge].prev =
                    temp_index_embed;

                index_embed = temp_index_embed;
                index_embed += 1;
            }
            else
            {
                (*vertices)[v_l].first_edge = NIL;
            }
            continue;
        }
        
        /*
          get v's orientation, and from this decide the way in which
          v's adjacency list will be traversed
          (recall that the list is supposed to be consistent, so no bad
          surprises)
        */
        orient = embed_graph[v].sign;
        in = orient == CCLOCKW ? 0 : 1;
        out = 1 ^ in;

        e = embed_graph[v].link[out];
        while (embedg_VES_is_short_cut_edge(embed_graph, n, e))
        {
            e = embed_graph[e].link[out];
        }
        ASSERT(embedg_VES_is_edge(n, e)
               && !embedg_VES_is_short_cut_edge(embed_graph, n, e));
        /*
          strictly speaking there should be no SCEs left at this stage...
          
          if there are SCEs in v's list, it must be the case that
          the list also contains tree or back edges...
        */
        
        (*vertices)[v_l].first_edge = index_embed;

        IF_DEB_EMBED(
                     fprintf(stdout, "recov. embed. DFI %d vertex %d at %d (edges) and %d (embedding)\n",
                             v, v_l, index_e, (*vertices)[v_l].first_edge);
                     )

        cur_e = e;
        while (TRUE)
        {
            next_e = embed_graph[cur_e].link[out];
            while (embedg_VES_is_short_cut_edge(embed_graph, n, next_e))
            {
                next_e = embed_graph[next_e].link[out];
            }
            ASSERT(!embedg_VES_is_short_cut_edge(embed_graph, n, next_e));

            if (next_e == v)
                /*
                  end of adjacency list
                */
            {
                break;
            }
            
            ASSERT(embedg_VES_is_edge(n, next_e));

            (*embedding)[index_embed].in_adjl = embed_graph[cur_e].in_adjl; 
            (*embedding)[index_embed].next = index_embed + 1; /* next in adj.
                                                                list */
            (*embedding)[index_embed].mark = NIL;  /* mark */
            
            /*
              cur_e's twin is trickier:
              we'll use twin's label field to store cur_e's index in
              the embedding

              if cur_e's label != NIL this means that cur_e's twin
              is already stored in edges/embedding  and consequently
              that cur_e.label = index of its twin (in the embedding)

              note that it is safe to do so since an edge's label
              has no meaning
            */
            if (embed_graph[cur_e].label != NIL)
            {
                (*embedding)[index_embed].inv = embed_graph[cur_e].label;

                /*
                  but fix the twin by the same token
                */
                (*embedding)[embed_graph[cur_e].label].inv = index_embed;
                ASSERT((*embedding)[embed_graph[cur_e].label].in_adjl ==
                       embed_graph[cur_e].twin_in_adjl);
            }
            else
                /*
                  we store cur_e's index in the embedding in twin's label
                */
            {
                int      twin;
                
                twin = embedg_VES_get_twin_edge(embed_graph, n, cur_e);
                embed_graph[twin].label = index_embed;
            }

            /*
              so the only thing we couldn't update yet is
              (*embedding)[index_embed].prev, cur_e previous edge in the list

              but we can do this for next_e
            */
            (*embedding)[index_embed + 1].prev = index_embed;

            /*
              we check if there are any multiple edges or loops
              to embed
            */
            w = embed_graph[cur_e].neighbour;
            mult = embed_graph[cur_e].mult - 1;
            /*
              one was for the TE or BE edge
            */

            if (index_embed == (*vertices)[v_l].first_edge)
                /*
                  when looking for multiple edges/loops
                  we must temporarily "close" this ordered
                  list of vertices when in presence of the first
                  edge in the list:

                  not doing this would mean that
                  (*embedding)[(*vertices)[v_l].first_edge].prev
                  contains some irrelevant value which may cause
                  (major) trouble when embedding inverses of
                  multiple edges...
                */
            {
                (*embedding)[(*vertices)[v_l].first_edge].prev = index_embed;
            }
            
            embedg_recover_embedding_embed_mult(mult_edges, *embedding,
                                                nbr_e, v, w, mult,
                                                &index_embed, &set_next,
                                                &new_first_edge);
            embedg_recover_embedding_embed_loops(mult_edges, *embedding,
                                                 nbr_e, v, &index_embed,
                                                 &set_next);
            set_next = TRUE;

            /*
              yes, it may be the case that (*vertices)[v_l].first_edge
              change while in embedg_recover_embedding_embed_mult
              -- see that function for more
            */
            (*vertices)[v_l].first_edge = new_first_edge == NIL ?
                (*vertices)[v_l].first_edge : new_first_edge;
            
            /*
              that's all, we proceed to read a new edge in the list
            */
            index_embed += 1;
            cur_e = next_e;
        }

        /*
          now next_e = v so that cur_e is the last edge in v's adjacency list
          we must deal with this case separately
        */

        /*
          fix cur_e in embedding (and its twin)
        */
        (*embedding)[index_embed].in_adjl = embed_graph[cur_e].in_adjl;

        /*
          we temporarily set next of cur_e in to index_embed + 1
        */
        (*embedding)[index_embed].next = index_embed + 1;
        (*embedding)[index_embed].mark = NIL;  /* mark */

        /*
          fix cur_e's twin
        */
        if (embed_graph[cur_e].label != NIL)
        {
            (*embedding)[index_embed].inv = embed_graph[cur_e].label;
            (*embedding)[embed_graph[cur_e].label].inv = index_embed;
            ASSERT((*embedding)[embed_graph[cur_e].label].in_adjl ==
                   embed_graph[cur_e].twin_in_adjl);
        }
        else
        {
            int      twin;
            
            twin = embedg_VES_get_twin_edge(embed_graph, n, cur_e);
            embed_graph[twin].label = index_embed;
        }

        /*
          we temporarily set the next record's prev field:
          but we can do that only if we haven't processed
          all the edges yet
        */
        if (index_embed < 2 * nbr_e - 1)
        {
            (*embedding)[index_embed + 1].prev = index_embed;

            /*
              again, check if there are any multiple edges/loops
              to embed
            */
            w = embed_graph[cur_e].neighbour;
            mult = embed_graph[cur_e].mult - 1;
            /*
              one was for the TE or BE edge
            */
            v_w_in_embed = index_embed;

            if (index_embed == (*vertices)[v_l].first_edge)
                /*
                  same comment as above
                */
            {
                (*embedding)[(*vertices)[v_l].first_edge].prev = index_embed;
            }
            
            embedg_recover_embedding_embed_mult(mult_edges, *embedding,
                                                nbr_e, v, w, mult,
                                                &index_embed, &set_next,
                                                &new_first_edge);
            embedg_recover_embedding_embed_loops(mult_edges, *embedding,
                                                 nbr_e, v, &index_embed,
                                                 &set_next);

            /*
              same comment as above
            */
             (*vertices)[v_l].first_edge = new_first_edge == NIL ?
                (*vertices)[v_l].first_edge : new_first_edge;
         }

        /*
          to finish off, we must set:
          
          cur_e's next field:
          next of cur_e in the list is ... vertices[v_l].first_edge
          
          cur_e's next's previous field...
        */
        if (set_next)
            /*
              set_next (poorly named) is used to indicate which
              edges must be updated to "close off" the list:

              if set_next is TRUE, we are in the standard case
              where the last edge in the ordered adj. list
              is at index_embed

              if set_next is FALSE, the last edge in the ordered adj. list
              is at v_w_in_embed: because it could have happened
              (in embedg_recover_embedding_embed_mult only)
              that the edges have been "wedged" between
              v_w_in_embed.prev and v_w_in_embed,
              leaving v_w_in_embed the last in the list
            */
        {
            (*embedding)[index_embed].next = (*vertices)[v_l].first_edge;
            (*embedding)[(*vertices)[v_l].first_edge].prev = index_embed;
        }
        else
        {
            (*embedding)[v_w_in_embed].next = (*vertices)[v_l].first_edge;
            (*embedding)[(*vertices)[v_l].first_edge].prev = v_w_in_embed;
        }
        set_next = TRUE;
            
        /*
          a simple check
        */
        ASSERT(embedg_dlcl_is_empty(mult_edges[v]));
        
        /*
          we can process another vertex
        */
        index_embed += 1;
    }
    /*
      when this is done there are a few things that must hold
    */
    ASSERT(index_embed == 2 * nbr_e);
}


static void 
embedg_recover_embedding_embed_mult (t_dlcl **mult_edges,
	t_embed_sparse_rep *embedding, int nbr_e, int v, int w,
	int mult, int *index_embed, boolean *set_next, int *first_edge)
    /*
      see if the directed edge [v, w] is multiple: if so embed it
      in embedding

      moreover if there are any [v, v] loops do that too
    */
{
    /*
      we take care of multiple edges: for tree edges and back
      edges their multiplicity is indicated by the
      embed_graph[cur_e].mult field (which records the number
      of undirected edges)
      
      for loops hovewer this information is stored in the mult
      field of the FIRST encountered neighbour v in v's neighbour
      list
    */
    t_dlcl      *p;
    int         v_w_in_embed, v_w_prev;
    boolean     do_twins, start, do_first_edge;

    IF_DEB_EMBED_MULT(
           fprintf(stdout, "in recover emb. mult, v %d w %d mult %d\n",
                   v, w, mult);
           )

    /*
      the current index_embed value is the edge [v, w]:
      I must record this value as it will be needed
      later
    */
    v_w_in_embed = *index_embed;
    start = TRUE;
    *set_next = TRUE;
    do_twins = FALSE;
    *first_edge = NIL;
    do_first_edge = FALSE;
    v_w_prev = NIL;
    while (mult > 0)
    {
        ASSERT(!embedg_dlcl_is_empty(mult_edges[v]));
        p = embedg_dlcl_find(mult_edges[v], w);
        /*
          note that using embedg_dlcl_find to always find
          the first in the list with p->info == w
          is ok here since any previous such records would
          have been deleted/removed from the list
        */
        ASSERT(p != NP);
        /*
          otherwise we couldn't have mult > 0 !
        */
        
        *index_embed += 1;
        
        /*
          once again I must use a similar sort of trick as in the
          main function to deal with the inverse edge:
          
          the inverse edge is to be found in mult_edges[w]:
          if p->twin_in_adjl (which was initialised to NIL
          and has NOT been set in the DFS preprocessing),
          if p->twin_in_adjl != NIL, then
          a. its inverse in mult_edges[w] has already been embedded
          in *embedding
          b. its index there is stored in p->twin_in_adjl
          precisely
         */
        if (p->twin_in_adjl != NIL)
        {
            if (! start)
                /*
                  if the first the multiple edges' inverse is already
                  stored, then this is true for ALL of them
                */
            {
                ASSERT(do_twins == TRUE);
            }
            do_twins = TRUE;
        }
        else
            /*
              similarly, if the first the multiple edges' inverse is
              not already stored, then this is true for ALL of them
            */
        {
            ASSERT(do_twins == FALSE);
        }

        embedding[*index_embed].in_adjl = p->in_adjl;
        embedding[*index_embed].mark = NIL;  

        /*
          as we will see do_twins has to be treated differently
        */
        if (!do_twins)
            /*
              this is pretty standard as works as the
              main recover function
            */
        {
            t_dlcl      *i_m_l, *i_p;
            
            embedding[*index_embed].next = *index_embed + 1;
            
            /*
              we store the current index in the embedding in
              the twin/inverse's twin_in_adjl field
            */
            i_p = i_m_l = mult_edges[w];
            ASSERT(!embedg_dlcl_is_empty(i_m_l));
            i_p = embedg_dlcl_find_with_NIL_twin_in_adjl(i_m_l, v);
            ASSERT(i_p != NP);
            ASSERT(i_p->twin_in_adjl == NIL);
            
            i_p->twin_in_adjl = *index_embed;
            
             /*
              to finish off this bit we set embedding[*index_embed + 1].prev
              
              but I can only set this prev field if I haven't reached
              the end of the embedding[] array: this is why we needed
              nbr_e (total number of edges to embed) as input
             */
            
            if (*index_embed < 2 * nbr_e - 1)
            {
                embedding[*index_embed + 1].prev = *index_embed;
            }
        }
        else
            /*
              how to insert the inverses of multiple edges already
              in the embedding:

              if one studies how the twin_in_adjl field has been
              set while dealing with the inverses of the
              present multiple edges one sees that
              the latter must be inserted in counter clockwise
              order (assuming that the inverses were inserted
              in clockwise order)

              this is necessariy to ensure a correct matching between
              the edge and its inverse
            */
        {
            
            embedding[*index_embed].inv = p->twin_in_adjl;
            
            /*
              fix the twin by the same token
            */
            embedding[p->twin_in_adjl].inv = *index_embed;

            /*
              general (reverse) insertion for these edges
            */
            embedding[*index_embed].prev = *index_embed + 1;
            embedding[*index_embed].next = *index_embed - 1;
 
            /*
              ok, that was the easy bit, things are a bit more complicated
              below...
            */
            if (start)
                /*
                  the edges are "wedged" between
                  embedding[v_w_in_embed].prev and v_w_in_embed,

                  hence the following
                */
            {
                v_w_prev = embedding[v_w_in_embed].prev;
                if (v_w_prev == v_w_in_embed)
                    /*
                      in this case the first edge in the adj. list
                      of the vertex whose first_edges is v_w_in_embed
                      will be changed
                      (because we insert in reverse order)
                    */
                {
                    do_first_edge = TRUE;
                }

                embedding[*index_embed].next = v_w_in_embed;
                embedding[v_w_in_embed].prev = *index_embed;
                
                ASSERT(embedding[embedding[*index_embed].inv].prev ==
                       embedding[v_w_in_embed].inv);
                ASSERT(embedding[embedding[v_w_in_embed].inv].next ==
                       embedding[*index_embed].inv);
            }
 
            if (mult == 1)
                /*
                  last inv. edge in this list to add
                */
            {
                ASSERT(v_w_prev != NIL);

                /*
                  must fix embedding[v_w_prev].next appropriately
                  (and embedding[*index_embed].prev)

                  this may be overwritten later on, but not necessarily so

                  the next_set flag will enable us to decide
                  which edge ends this adjacency list: see above
                */
                
                embedding[*index_embed].prev = v_w_prev;
                embedding[v_w_prev].next = *index_embed;
                *set_next = FALSE;

                ASSERT(embedding[embedding[*index_embed].inv].prev ==
                       embedding[*index_embed - 1].inv);
                ASSERT(embedding[embedding[*index_embed - 1].inv].next ==
                       embedding[*index_embed].inv);

                if (do_first_edge)
                    /*
                      the first edge is the last one added
                    */
                {
                    *first_edge = *index_embed;
                }

                embedding[v_w_in_embed].next = *index_embed + 1;
                if (*index_embed < 2 * nbr_e - 1)
                {
                   embedding[*index_embed + 1].prev = v_w_in_embed;
                }
            }

            ASSERT(embedding[embedding[*index_embed].inv].prev ==
                   embedding[embedding[*index_embed].next].inv);
        }
        
        /*
          to finish off this bit we delete the p record from m_l
          and set embedding[*index_embed + 1].prev
        */
        mult_edges[v] = embedg_dlcl_delete_rec(mult_edges[v], p);

        mult--;
        start = FALSE;
    }
    /*
      conclusion: sevral days to get this working! *sigh*
    */
}




static void 
embedg_recover_embedding_embed_loops (t_dlcl **mult_edges,
	t_embed_sparse_rep *embedding, int nbr_e, int v,
	int *index_embed, boolean *set_next)
    /*
      embed the [v, v] loops 
    */
{
    /*
      the loops' multiplicity is stored in the mult
      field of the FIRST encountered neighbour v in v's neighbour
      list
    */
    t_dlcl      *p;
    int         nbr_loops;
    
    /*
      have a look if there are any [v. v] loops
    */
    p = embedg_dlcl_find(mult_edges[v], v);
    if (p == NP)
    {
        return;
    }

    /*
      when there are loops to add to the adjaceny list,
      edge insertion resume in the "normal" clockwaise saya, way:
      so we reset set_next to true
    */
    *set_next = TRUE;
    
    nbr_loops = p->mult;
    ASSERT(nbr_loops % 2 == 0);
    /*
      we counted directed edges
    */
    nbr_loops /= 2;

    IF_DEB_EMBED_LOOPS(
           fprintf(stdout, "in recover emb. loops, nbr_loops [v, v] %d\n",
                   nbr_loops);
           )
        
    while (nbr_loops > 0)
        /*
          a loop requires to embed two directed edges
        */
    {
        p = embedg_dlcl_find(mult_edges[v], v);
        ASSERT(p != NP);

        *index_embed += 1;
        
        embedding[*index_embed].in_adjl = p->in_adjl;
        embedding[*index_embed].next = *index_embed + 1; 
        embedding[*index_embed].mark = NIL; 
        embedding[*index_embed].inv = *index_embed + 1;
        embedding[*index_embed + 1].prev = *index_embed;
        
        mult_edges[v] = embedg_dlcl_delete_rec(mult_edges[v], p);
        
        IF_DEB_EMBED_LOOPS(
           fprintf(stdout, "in recover emb. loops, mid\n");
           embedg_dlcl_print(mult_edges[v]);
           );

        /*
          now do the "inverse" loop
        */
        p = embedg_dlcl_find(mult_edges[v], v);
        ASSERT(p != NP);
        
        *index_embed += 1;
        
        embedding[*index_embed].in_adjl = p->in_adjl;
        embedding[*index_embed].next = *index_embed + 1; 
        embedding[*index_embed].mark = NIL;
        embedding[*index_embed].inv = *index_embed - 1;

        if (*index_embed < 2 * nbr_e - 1)
        {
            embedding[*index_embed + 1].prev = *index_embed;
        }
        mult_edges[v] = embedg_dlcl_delete_rec(mult_edges[v], p);
        
        nbr_loops--;

        IF_DEB_EMBED_LOOPS(
           fprintf(stdout, "in recover emb. loops, end\n");
           embedg_dlcl_print(mult_edges[v]);
           );
    }
}




void 
embedg_recov_embed_walk_proper_face (int n, int e, t_adjl_sparse_rep *A,
	t_embed_sparse_rep *embedding, boolean MARK, int mark)
    /*
      do a proper face walk in the recovered embedding starting
      at index e in the embedding
    */
{
    int          cur, next;
    
    IF_DEB_FACES(
                 fprintf(stdout, "recov. emb. proper face walk\n");
                 fprintf(stdout, "[-, %d] ",
                         A[embedding[e].in_adjl].end_vertex);
                 )
        
    cur = e;
    next = NIL;
    while (next != e)
        /*
          to get the next in a proper face traversal:
          get the previous of the cur's inverse
        */
    {
        int     inv;

        inv = embedding[cur].inv;
        next = embedding[inv].prev;

        ASSERT(embedding[next].mark != mark);
        
        if (MARK)
        {
            embedding[next].mark = mark;
        }
        
        cur = next;
        IF_DEB_FACES(
                     fprintf(stdout, "[-, %d] ",
                             A[embedding[cur].in_adjl].end_vertex);
                     )
    }
    IF_DEB_FACES(
                 fprintf(stdout, "\n");
                 )
}
      


boolean 
embedg_check_recov_embedding (int n, int nbr_e, int nbr_comp,
	t_ver_sparse_rep *vertices, t_adjl_sparse_rep *A,
	t_embed_sparse_rep *embedding)
    /*
      check if the recovered embedding is a valid embedding
      SHOULD ONLY be use after creation, that is, after having
      recovered the embedding from the VES structure
      (because of the mark MIN_EMBED_MARK we use)
    */
{
    int          v, e, f;

    f = 0;
    /*
      do all the edges in embedding:
      careful: we have 2 * nbr_e to visit (the edge and its inverse!)
    */
    for (e = 0; e < 2 * nbr_e; e++)
    {
        /*
          we check if the current edge is marked: if not, we
          traverse a proper face bordered by this edge
        */
        if (embedding[e].mark != MIN_EMBED_MARK)
            /*
              we --hopefully-- perform this check only after creation
              where mark == NIL
            */
        {
            embedg_recov_embed_walk_proper_face(n, e, A, embedding,
                                                TRUE, MIN_EMBED_MARK);
            f++;
        }
    }

    /*
      must also count a face for each isolated vertex
    */
    for (v = 0; v < n; v++)
    {
        if (vertices[v].first_edge == NIL)
            f++;
    }
    
    IF_DEB_CHECK_EMBED(
                       fprintf(stdout, "recovered embedding, n: %d\t e: %d\t C: %d\t f: %d\n",
                               n, nbr_e, nbr_comp, f);
                       )
        
    return f == 2 * nbr_comp + nbr_e - n ? TRUE : FALSE;
}


t_dlcl **
embedg_recover_obstruction (t_ver_edge *embed_graph, int n, minor m, int *nbr_e)
    /*
      recover the obstruction as a t_dlcl * structure:
      and return the number of edges: lets say we agree on returning
      the number of undirected edges
      -- I don't know yet which way to do, directed or undirected???

      so far in the algorithm we only dealt with DFIs,
      but now, we retrieve the obstruction not wrt DFIs but
      wrt the vertices' labels
    */
{
    /*
      so I am looking, in embed_graph, for the vertices and edges
      marked MARK_MINORS(n)
    */

    int          v;
    t_dlcl       **obs;

    obs = (t_dlcl **) mem_malloc(sizeof(t_dlcl *) * n);
    for (v = 0; v < n; v++)
        obs[v] = NP;

    *nbr_e = 0;
    for (v = 0; v < 2*n; v++)
        /*
          must check real vertices as well as virtual vertices
        */
    {
        int      e;

        if (embed_graph[v].link[0] == v)
            /*
              isolated vertex case
            */
        {
            ASSERT(embed_graph[v].link[1] == v);
            continue;
        }

        e = embed_graph[v].link[0];
        while (e != v)
        {
            ASSERT(embedg_VES_is_edge(n, e));
            if (embed_graph[e].visited == MARK_MINORS(n))
            {
                int        cur_v, neigh;

                /*
                  virtual vertices may still hang around
                */
                /*
                  let's get the "actual" v:
                  note that the statement below is safe since if v were
                  not a valid virtual vertex (ie [v - n].DFS_parent = n)
                  it would have an empty
                  adjacency list and we wouldn't be there anyway
                */
                cur_v = embedg_VES_get_ver(embed_graph, n, v);
        
                neigh = embedg_VES_get_ver(embed_graph, n,
                                           embed_graph[e].neighbour);

                /*
                  again, cur_v and neigh are DFIs,
                  we want vertex labels at this stage
                */
                cur_v = embed_graph[cur_v].label;
                neigh = embed_graph[neigh].label;
                sparseg_dlcl_append_to_neigh_list(obs, n, cur_v, neigh,
                                                  embed_graph[e].in_adjl);
                (*nbr_e)++;
            }
            e = embed_graph[e].link[0];
        }
    }
    
    IF_DEB_OBS(
               fprintf(stdout, "recovering the obstruction\n");
               sparseg_dlcl_print(obs, n);
    );

    ASSERT(*nbr_e % 2 == 0);
    *nbr_e /= 2;
    
    return obs;
}


static t_dlcl **
embedg_get_reduced_obs (t_dlcl **obs, int n)
    /*
      reduce the obstruction by removing all degree 2 vertices
      (so that they become isolated vertices)
    */
{
    t_dlcl       **reduced;
    int          v;

    reduced =  (t_dlcl **) mem_malloc(sizeof(t_dlcl *) * n);
    for (v = 0; v < n; v++)
    {
        reduced[v] = embedg_dlcl_copy(obs[v]);
    }
        
    for (v = 0; v < n; v++)
    {
        t_dlcl   *n_l, *n_l_b, *p, *new_n_v, *n_l_x, *b_in_n_x;
        int      a, b, n_x;

        n_l = reduced[v];
        while (!embedg_dlcl_is_empty(n_l)
                && embedg_dlcl_list_last(n_l) == embedg_dlcl_list_next(n_l))
            /*
              pick out which  vertices have deg 2
            */
        {
            a = n_l->info;
            b = embedg_dlcl_list_next(n_l)->info;
            /*
              we remove the edge (v, b), or rather, we identify v and b:
              b will then be an isolated vertex
              
              fix v's neighbour list: all of b's neighbours
              are now v's neighbours
            */
            reduced[v] = n_l =
                embedg_dlcl_delete_rec(n_l, embedg_dlcl_list_last(n_l));
            
            p = n_l_b = reduced[b];
            ASSERT(!embedg_dlcl_is_empty(n_l_b));
            n_x = p->info;
            if (n_x != v)
            {
                new_n_v = embedg_dlcl_rec_new(n_x);
                reduced[v] = n_l = embedg_dlcl_cat(n_l, new_n_v);
                
                /*
                  and in n_x neighbour list, we must replace b by v
                */
                n_l_x = reduced[n_x];
                b_in_n_x = embedg_dlcl_find(n_l_x, b);
                b_in_n_x->info = v;
            }
            /*
              and do this for all of b's neighbours
            */
            p = embedg_dlcl_list_next(p);
            while (p != n_l_b)
            {
                n_x = p->info;
                if (n_x != v)
                {
                    new_n_v = embedg_dlcl_rec_new(n_x);
                    reduced[v] = n_l = embedg_dlcl_cat(n_l, new_n_v);
                    n_l_x = reduced[n_x];
                    b_in_n_x = embedg_dlcl_find(n_l_x, b);
                    b_in_n_x->info = v;
                }
                p = embedg_dlcl_list_next(p);
            }
            embedg_dlcl_delete(reduced[b]);
            reduced[b] = NP;
        }
    }

    IF_DEB_CHECK_OBS(
                     fprintf(stdout, "reducing the obstruction\n");
                     sparseg_dlcl_print(reduced, n);
                     )

    /*
      now check no degree 2 vertices are left
    */
    for (v = 0; v < n; v++)
    {
        t_dlcl   *n_l;
        
        n_l = reduced[v];
        if (!embedg_dlcl_is_empty(n_l))
        {
            ASSERT(embedg_dlcl_list_last(n_l) != embedg_dlcl_list_next(n_l));
        }
    }

    return reduced;
}

static boolean 
embedg_is_red_obs_K33 (t_dlcl **reduced, int n)
    /*
      check if the (reduced) obstruction is indeed K33 
    */
{
    int          v, order, vs[6], i, b1[3];
    
    /*
      check that order == 6 and that the obstruction is cubic
    */
    order = 0;
    for (v = 0; v < n; v++)
    {
        if (!embedg_dlcl_is_empty(reduced[v]))
        {
            if (order == 6)
            {
                return FALSE;
            }
            order++;
            vs[order - 1] = v;
            
            if (embedg_dlcl_length(reduced[v]) != 3)
            {
                return FALSE;
            }
        }
    }
    if (order != 6)
    {
        return FALSE;
    }
    
    /*
      check if bipartite
    */
    v = vs[0];
    ASSERT(!embedg_dlcl_is_empty(reduced[v]));
    b1[0] = reduced[v]->info;
    b1[1] = embedg_dlcl_list_next(reduced[v])->info;
    b1[2] = embedg_dlcl_list_prev(reduced[v])->info;
    
    for (i = 1; i < 6; i++)
    {
        t_dlcl      *n_v;

        v = vs[i];
        n_v = reduced[v];
        ASSERT(!embedg_dlcl_is_empty(n_v));
        if (n_v->info == b1[0]
            || embedg_dlcl_list_next(n_v)->info == b1[0]
            || embedg_dlcl_list_prev(n_v)->info == b1[0])
        {
            if ((n_v->info != b1[1]
                 && embedg_dlcl_list_next(n_v)->info != b1[1]
                 && embedg_dlcl_list_prev(n_v)->info != b1[1])
                &&
                (n_v->info != b1[2]
                 && embedg_dlcl_list_next(n_v)->info != b1[2]
                 && embedg_dlcl_list_prev(n_v)->info != b1[2]))
            {
                return FALSE;
            }
        }
        else
        {
            if ((n_v->info == b1[1]
                 || embedg_dlcl_list_next(n_v)->info == b1[1]
                 || embedg_dlcl_list_prev(n_v)->info == b1[1])
                ||
                (n_v->info == b1[2]
                 || embedg_dlcl_list_next(n_v)->info == b1[2]
                 || embedg_dlcl_list_prev(n_v)->info == b1[2]))
            {
                return FALSE;
            }
        }
    }

    return TRUE;
}


static boolean 
embedg_is_red_obs_K5 (t_dlcl **reduced, int n)
    /*
      check if the (reduced) obstruction is indeed K5
    */
{
    int          v, order;
    
    /*
      check that order == 5 and that the obstruction is quadric
    */
    order = 0;
    for (v = 0; v < n; v++)
    {
        if (!embedg_dlcl_is_empty(reduced[v]))
        {
            if (order == 5)
            {
                return FALSE;
            }
            order++;
            
            if (embedg_dlcl_length(reduced[v]) != 4)
            {
                return FALSE;
            }
        }
    }

    return TRUE;
}


boolean 
embedg_check_recov_obs (t_dlcl **obs, int n, minor m)
    /*
      check if the recovered obstruction is one of K33 or K5
    */
{
    t_dlcl      **reduced;
    boolean     ans;

    reduced = embedg_get_reduced_obs(obs, n);
    if (m != MINOR_E5)
    {
        ans = embedg_is_red_obs_K33(reduced, n);
    }
    else
    {
        ans = embedg_is_red_obs_K5(reduced, n);
    }

    sparseg_dlcl_delete(reduced, n);
    return ans;
}
/*
 *  obstruction.c
 */
 
/*
  What:
  *****
  
  Implementing:

  The graph is not planar: we recover the obstruction from the VES structure
  and check it as well.
  (Some of these checks will disappear later)
  


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
#define IF_DEB_OBS(x) {}
#define IF_DEB_CHECK_OBS(x) {}
#define IF_CPU(x) {}
#define IF_DEB_MINOR(x) {}
 
 
/* aproto: header embed_graph_protos.h */

void 
embedg_obstruction (
    t_ver_sparse_rep *V,
    t_adjl_sparse_rep *A,       /* the input graph as a sparse graph */
    t_dlcl **dfs_tree,      /* a sparse graph rep. for the dfs tree
                                      -- vertices are as DFIs
                                      -- and children are ordered wrt
                                      lowpoint value
                                  */
    t_dlcl **back_edges,    /* for each vertex v, a dlcl
                                      of the back edges [v, x] incident to v
                                      where x is a DESCENDANT of v
                                      (vertices are given as DFIs)
                                  */
    t_ver_edge *embed_graph,    /* output of tester */
    int n,               /* size of the graph */
    int *edge_pos,       /* pos. in embed_graph for addition
                                     of the next edge */
    int v,
    int w_in,         /* the unembedded directed back edge
                                     [w_in, v]
                                   */
    t_ver_sparse_rep **OV,      /* the obstruction as an adjacency list */
    t_adjl_sparse_rep **OA,
    int *nbr_e_obs      /* obstruction's #edges */
)
                                    
    /*
      the graph is non planar: we must mark & recover the K33 or K5
      homeomorph
    */
{
    int          *ver_orient;
    minor        m;
    t_dlcl       **obs;
    
    /*
      this is magma code - must be removed
    */
    float      sttime, time_to_now;
 
 IF_CPU(
    sttime = time_current_user();
       )

    /*
      we will NOT remove the short-cut edges at this stage:
      we'll have to perform another walkdown in embedg_iso_is_minor_A
      so
      1. saves time when looking for ext. active vertices
      2. more importantly this enables us to ascertain that the number of
         edges in embed_graph (even after completing whichever obstruction
         applying in this case) will NEVER be > 3*n - 5!!!
      3. SCEs are then removed in embedg_iso_is_minor_A
         (obligatory path for every possible case)
    */

    /*
      we must compute each vertex's orientation (wrt flipped bicomps)
      and set the edges' orientation:
      
      the other day I was wondering why this was necessary in this
      instance (because after all we won't get an embedding):
      orientation is required bacause later in the piece we
      do a proper face traversal (I guess for Minor C testing)
    */
    ver_orient = embedg_vertices_orientation(embed_graph, n);
    embedg_VES_set_orientation(embed_graph, n, ver_orient);
    mem_free(ver_orient);

    m = embedg_mark_obstruction(dfs_tree, back_edges,
                                embed_graph, n, edge_pos, v, w_in);

    /*
      get the obstruction
    */
    obs = embedg_recover_obstruction(embed_graph, n, m, nbr_e_obs);

    /*
      and check it
    */
    if (!embedg_check_recov_obs(obs, n, m))
    {
        sparseg_dlcl_delete(obs, n);
        DIE();
    }
    
    sparseg_dlcl_to_sparseg(obs, n, *nbr_e_obs, OV, OA);
    sparseg_dlcl_delete(obs, n);

    /*
      just for the sake of it, chcek if the obstruction is
      a subgraph of the input graph
    */
    if (!sparseg_adjl_sub(*OV, n, *OA, V, n, A))
    {
        DIE();
    }
    
    IF_DEB_OBS(
               sparseg_adjl_print(*V, n, *A, FALSE);
               )
    
    IF_CPU(
           fprintf(stdout, "CPU for obstruction recovering %f\n",
                   (time_current_user() - sttime));
           )
}







minor 
embedg_mark_obstruction (
    t_dlcl **dfs_tree,      /* a sparse graph rep. for the dfs tree
                                      -- vertices are as DFIs
                                      -- and children are ordered wrt
                                      lowpoint value
                                  */
    t_dlcl **back_edges,    /* for each vertex v, a dlcl
                                      of the back edges [v, x] incident to v
                                      where x is a DESCENDANT of v
                                      (vertices are given as DFIs)
                                  */
    t_ver_edge *embed_graph,    /* output of tester */
    int n,               /* size of the graph */
    int *edge_pos,       /* pos. in embed_graph for addition
                                     of the next edge */
    int v,
    int w_in         /* the unembedded directed back edge
                                     [w_in, v]
                                   */
)
    /*
      the graph is non planar: we must mark & recover the K33 or K5
      homeomorph
    */
{
    int          c, vr, x, y, w;
    int          *path_v, *path_e, nbr_v, entry_in_path_e;
    boolean      px_attached_high, py_attached_high, is_minor_D;
    minor        m;
    
 
  IF_CPU(
    float      sttime; float time_to_now;

    sttime = time_current_user();
        )
    

    /*
      find c such that v^c is the root of the biconnected
      component on which the walkdown failed
    */
    c = embedg_iso_get_c_of_v(embed_graph, n, v, w_in);

    /*
      now: decide which minor we are dealing with and mark the
      appropriate one (vertices/edges marked as MARK_MINOR(n)
      in embed_graph)
    */
    if (embedg_iso_is_minor_A(embed_graph, n, edge_pos, v, c, &vr))
    {
        embedg_mark_minor_A(dfs_tree, back_edges,
                            embed_graph, n, edge_pos, v, c, vr);
        
        IF_DEB_MINOR(
                     fprintf(stdout, "Minor A\n");
                     )

        return MINOR_A;
    }

    /*
      get the externally active vertices x & y and the pertinent w
      on the external face of the bicomp rooted by v^c

      and determine if minor B
    */
    if (embedg_iso_is_minor_B(embed_graph, n, edge_pos, v, c,
                                   &x, &y, &w))
    {
        embedg_mark_minor_B(dfs_tree, back_edges,
                            embed_graph, n, edge_pos, v, c,
                            x, y, w);
        IF_DEB_MINOR(
                     fprintf(stdout, "Minor B\n");
                     )

        IF_CPU(
               fprintf(stdout, "CPU for obstruction isolation %f\n",
                       time_current_user() - sttime);
               )
            
        return MINOR_B;
    }

    /*
      the remaining cases: must get the highest x-y path

      it will be containing in path_v (vertices), path_e (edges)
    */
    embedg_iso_get_highest_x_y_path(embed_graph, n, MARK_EXT_FACE(n),
                                    MARK_EXT_FACE_L(n),
                                    MARK_EXT_FACE_R(n),
                                    v, c, x, y, w,
                                    &path_v, &path_e,
                                    &nbr_v, &entry_in_path_e,
                                    &px_attached_high,
                                    &py_attached_high,
                                    &is_minor_D);

    /*
      we are in the minor C case if either one of p_x or p_y
      is attached high
    */
    if (px_attached_high || py_attached_high)
    {
        embedg_mark_minor_C(dfs_tree, back_edges, embed_graph, n, edge_pos,
                            v, c, x, y, w,
                            path_v, path_e, nbr_v,
                            px_attached_high, py_attached_high);
        IF_DEB_MINOR(
                     fprintf(stdout, "Minor C\n");
                     )

        mem_free(path_v);
        mem_free(path_e);
        
        IF_CPU(
               fprintf(stdout, "CPU for obstruction isolation %f\n",
                       time_current_user() - sttime);
               )
            
        return MINOR_C;
    }

    if (is_minor_D)
    {
        embedg_mark_minor_D(dfs_tree, back_edges, embed_graph, n, edge_pos,
                            v, c, x, y, w,
                            path_v, path_e, nbr_v, entry_in_path_e);
        IF_DEB_MINOR(
                     fprintf(stdout, "Minor D\n");
                     )

        mem_free(path_v);
        mem_free(path_e);
        
        IF_CPU(
               fprintf(stdout, "CPU for obstruction isolation %f\n",
                       time_current_user() - sttime);
               )
            
        return MINOR_D;
    }

    /*
      finally, the minor E case
    */
    m = embedg_mark_minor_E(dfs_tree, back_edges, embed_graph, n, edge_pos,
                            v, c, x, y, w,
                            path_v, path_e, nbr_v);
    switch (m)
    {
    case MINOR_E1:
        IF_DEB_MINOR(
                     fprintf(stdout, "Minor E1\n");
                     )
        break;
    case MINOR_E2:
        IF_DEB_MINOR(
                     fprintf(stdout, "Minor E2\n");
                     )
        break;
    case MINOR_E3:
        IF_DEB_MINOR(
                     fprintf(stdout, "Minor E3\n");
                     )
        break;
    case MINOR_E4:
        IF_DEB_MINOR(
                     fprintf(stdout, "Minor E4\n");
                     )
        break;
    case MINOR_E5:
        IF_DEB_MINOR(
                     fprintf(stdout, "Minor E5\n");
                     )
        break;
    case MINOR_A: 
    case MINOR_B: 
    case MINOR_C: 
    case MINOR_D:
    case MINOR_E:
    case NBR_MINORS:
        break;
    }

    mem_free(path_v);
    mem_free(path_e);
    
    IF_CPU(
           fprintf(stdout, "CPU (scaled) for obstruction isolation %f\n",
                   (time_current_user() - sttime) / e);
           )
        
    return m;
}
/*
 *  isolator.c
 */
 
/*
  What:
  *****
  
  Implementing:

  The graph is non planar: we isolate the obstruction.
  
 
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 
 
  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca
 
 
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********
 
  Paulette Lieby (Magma), Brendan McKay (ANU)
 
  Started October 2001
*/
 
 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
#define IF_DEB_TREE(x)    {}
#define IF_DEB_EDGES(x) {}
#define IF_CPU(x) {}
/* #define IF_DEB_MINOR(x) {x}  -- Not Used  */
 
 
/* aproto: header embed_graph_protos.h */
 
#ifndef PLANAR_IN_MAGMA
#endif
 



int 
embedg_iso_get_c_of_v (t_ver_edge *embed_graph, int n, int v, int w)
    /*
      the edge [v, w] (w a descendant of v) remains unembedded
      after the walkdown returns

      find c such that v^c is the root of the biconnected
      component on which the walkdown failed
    */
{
    /*
      how to do this??? easy! follow the DFS tree path as given
      by the field DFS_parent
    */

    int           u;

    u = embed_graph[w].DFS_parent;
    while (embed_graph[u].DFS_parent != v)
    {
        u = embed_graph[u].DFS_parent;
    }
    /*
      this is guaranteed to succeed given the structure of the DFS tree
      and the fact that there exists a  back edge [w, v]
    */
    
    return u;
}


boolean 
embedg_iso_is_minor_A (t_ver_edge *embed_graph, int n,
	int *edge_pos, int v, int c, int *vr)
    /*
      determines if the obstruction is a minor A
    */
{
    /*
      to do this we again call the walkdown routine with v^c as input,
      the walkdown routine will fail (since there will be an
      un-embedded back edge incident to v and to a vertex
      in the subtree rooted by v^c)

      the obstruction is a minor A if the merge queue returned by the
      walkdown is non-empty, if this is the case we return
      the bicomp last appended to the queue
    */
    int             vv;
    t_merge_queue   q;

    vv = c + n;
    
    q = embedg_walkdown(embed_graph, n, edge_pos, vv);
    /*
      we MUST remove the SCEs here: this is the only place where it
      will be done when looking for and recovering an obstruction

      this is safe since this very function applies to ALL cases!
    */
    embedg_remove_SCE(embed_graph, n, *edge_pos);

    if (!embedg_merge_queue_empty(q))
        /*
          the bicomp of interest is the last in the queue
        */
    {
        int         r, rin, vrout;
        
        embedg_merge_queue_prune(&q, &r, &rin, vr, &vrout);
        embedg_merge_queue_delete(q);
        return TRUE;
    }
    else
    {
        embedg_merge_queue_delete(q);
        return FALSE;
    }
}


void 
embedg_iso_get_x_y_w (t_ver_edge *embed_graph, int n, int v, int r,
	int c, int mark, int mark_l, int mark_r, int *x, int *y, int *w)
    /*
      the obstruction is one of minor B, C, D, E.
      
      get the externally active vertices x & y along the
      external face  paths starting at r^c

      get a pertinent vertex w along the lower external
      face path between x and y

      external activity and pertinence are wrt v
      
      all the vertices on the external face r^c...x...w
      and r^c...y...w will be marked (the visited field)
    */
{
    int          vr, vrin, x_y[4];
    int          s, sin, cur, curin;

    vr = c + n;

    /*
      find x and y first:

      note that we mark the vertices on the external face r^c...x
      and r^c...y

      more on that below
    */
    embed_graph[vr].visited = mark;
    for (vrin = 0; vrin <= 1; vrin++)
    {
        int      m;

        m = vrin == 0 ? mark_l : mark_r;
        embedg_VES_get_succ_ext_active_on_ext_face(embed_graph, n, v,
                                                        vr, vrin,
                                                        TRUE, m,
                                                        &s, &sin);
        x_y[vrin] = s;
        x_y[vrin + 2] = sin;
        /*
          note the bizarre way I store the active vertex
          and the direction out of which to continue a walk
          on the lower external face as described above
        */
    }
    *x = x_y[0];
    *y = x_y[1];
    
    /*
      next get the pertinent w on the lower external face from x to y
    */
    cur = x_y[0];
    curin = x_y[2];
    embedg_VES_get_succ_pertinent_on_ext_face(embed_graph, n, v,
                                                   cur, curin,
                                                   TRUE, mark_l, w, &sin);

    /*
      now all the vertices  r^c...x...w and r^c...y have been marked,
      it remains to mark the vertices on the y...w external face path

      (will need to be able to distinguish the external face  later on)

      Note the way the external face is marked (needed when recovering
      the highest x-y path):
      mark_l for the path v^c...x...w
      mark_r for the path v^c...y
      mark for the lower external face y...w
    */
    cur = x_y[1];
    curin = x_y[3];
    s = n;
    while (s != *w)
    {
        embedg_VES_get_succ_pertinent_on_ext_face(embed_graph, n, v,
                                                       cur, curin,
                                                       TRUE, mark, &s, &sin);
        cur = s;
        curin = sin;
    }

    IF_DEB(
           fprintf(stdout, "get x, y & w: the external face\n");
           fprintf(stdout, "%d\t", vr);
           cur = vr;
           curin = 0;
           while (s != vr)
           {
               embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                                    cur, curin,
                                                    FALSE, 0, &s, &sin);
               cur = s;
               curin = sin;
               fprintf(stdout, "%d\t", s);
           }
           fprintf(stdout, "\n");
           )
}




boolean 
embedg_iso_is_minor_B (t_ver_edge *embed_graph, int n, int *edge_pos,
	int v, int c, int *x, int *y, int *w)
    /*
      determines if the obstruction is a minor B and return x, y
      (ext. active) and w (pertinent)
    */
{
    /*
      get x & y the ext. active vertices on the (external face)
      path out of v^c,
      and w the pertinent vertex on the lower external face x-y

      PLUS mark the whole external face with MARK_EXT_FACE(n)
    */
    embedg_iso_get_x_y_w(embed_graph, n, v, v, c,
                              MARK_EXT_FACE(n),
                              MARK_EXT_FACE_L(n), MARK_EXT_FACE_R(n),
                              x, y, w);

    if (embedg_dlcl_is_empty(embed_graph[*w].pertinent_bicomp_list))
        /*
          w has no pertinent child bicomp: not a minor B
        */
        return FALSE;
    else
    {
        t_dlcl      *pert_l;
        int         l;

        pert_l = embed_graph[*w].pertinent_bicomp_list;
        l = embedg_dlcl_list_last(pert_l)->info;
        /*
          if w has an ext. active pertinent child bicomp then minor B

          note that we need to know if w has an ext. active AND pertinent
          bicomp child: so it is NOT good enough to test
          w's separated_DFS_child_list as is done in
          embedg_VES_is_ver_ext_active!!!!!!!!!

          PLUS: l is actually a VIRTUAL vertex: to check its lowpoint
          I must take its DFS child l - n !!!!!!!!
        */
        ASSERT(embedg_VES_is_virtual_vertex(n, l));
        l = l - n;
        return embed_graph[l].lowpoint < v ? TRUE : FALSE;
    }
}

void 
embedg_iso_get_highest_x_y_path (
    t_ver_edge *embed_graph,
    int n,
    int mark,
    int mark_l,
    int mark_r,
    int v,
    int c,
    int x,
    int y,
    int w,
    int **path_v,     /* stack of vertices in x-y path */
    int **path_e,     /* stack of egdes in x-y path */
    int *nbr_v,         /* number of vertices in path_v */
    int *entry_in_path_e, /* the in direction for the FIRST edge in
                                      path_e: needed later on *sigh*
                                   */
    boolean *px_attached_high,
    boolean *py_attached_high,
    boolean *is_minor_D
)
    /*
      the obstruction is one of minor C, D, E.

      we want to recover the highest x-y path:
      the obstructing path attached to the external faces v^c - x - w
      and v^c - y - w

      while doing all this we also determine if the case is a minor C
      or a minor D
    */
{
    /*
      the path is obtained by walking the proper face starting at v
      where ALL the edges incident to v^c BUT the ones bordering
      the external face have been removed

      I won't I don't think remove these edges, but instead I'll be
      implementing an "avoidance" walk
    */

    int          vv, s, sin, p_x, p_y, cur_v, cur_vin;
    int          e, ein, s_e, s_ein;
    boolean      avoid_vv;
    
    /*
      must start the walk at edge embed_graph[v^c].link[1 ^ 0],
      (vvin = 0 is in direction of x, see embedg_iso_get_x_y_w)
    */
    vv = c + n;
    e = embed_graph[vv].link[1];
    ein = 0;     /* because of adjacency list consistency */

    *path_v = (int *) mem_malloc(sizeof(int) * n);
    *path_e = (int *) mem_malloc(sizeof(int) * n);
    (*nbr_v) = -1;
    
    /*
      recall that in embedg_iso_get_x_y_w we did mark
      (with mark, mark_l, mark_r)
      ALL the vertices lying on the external face walk starting
      & ending at v^c: we will use this fact to enable us
      to decide if a vertex is on the external face
      (as opposed to being on the internal face)
    */

    s = embed_graph[e].neighbour;
    ASSERT(embed_graph[s].visited == mark_l);
    /*
       this must be the case since s lies on the external face
       starting at v^c in x's direction
      -- we push s onto the stack
    */
    (*path_v)[++(*nbr_v)] = s;

    /*
      start the proper face walk which "avoids" v^c since the
      internal edges incident to v^c are supposed to have
      been removed

      please read on
    */
    avoid_vv = FALSE;
    while (TRUE)
    {
        boolean      av;
        
        av = 
            embedg_VES_get_succ_on_proper_face_with_avoidance(
                                                          embed_graph, n,
                                                          e, ein, vv,
                                                          FALSE, 0,
                                                          &s, &s_e, &s_ein);
        avoid_vv = av == TRUE ? av : avoid_vv;
        if (embed_graph[s].visited == mark_l)
            /*
              means that s is still on the external face:
              empty the path's stack and push s
            */
        {
            (*nbr_v) = -1;
            (*path_v)[++(*nbr_v)] = s;
            e = s_e;
            ein = s_ein;
        }
        else if (*nbr_v == 0)
            /*
              s is the first encountered vertex after
              path_v[0] which does not
              lie on the external face v^c...c...w

              given the way we pushed things on the vertex stack, path_v[0]
              will be the point of attachement of the x-y path
              on the v^c...x...w external face

              path_e[0] will contain nothing: a dummy
              
              path_e[1] will be the first edge in the x-y path
              (and entry_in_path will give the in-direction to this edge)
              
              oh yes!, we break the loop at this point if
              the vertex s lies on the v^c...y...w external face
            */
        {
            ASSERT(embed_graph[(*path_v)[0]].visited == mark_l);
            /*
              the first vertex on the path must be on the
              v^c...x...w external face
            */
            (*path_v)[++(*nbr_v)] = s;
            /*
              and now we also push the edge on the edge stack

              I'll need this later to initiate a proper face walk
              starting at the first vertex/edge in the x-y path,
              which is the same as starting from s_e
            */
            (*path_e)[*nbr_v] = s_e;
            *entry_in_path_e = s_ein;
            e = s_e;
            ein = s_ein;

            /*
              since we are at the start of the path, we must not
              forget to reset avoid_vv
            */
            avoid_vv = FALSE;
            
            if (embed_graph[s].visited == mark_r
                || embed_graph[s].visited == mark)
                /*
                  we have reached the v^c...y...w external face:
                  we can stop here
                */
            {
                break;
            }

            /*
              if not finished yet,
              we also mark s (and path_v[0]) as visited:
              later on we'll need to recognise which of the vertices
              in path have already been encountered
              (in case of encountering a cut-vertex due to the
              "removal" of the "internal" edges incidnet ot v^c)

              note that we mark s as visited iff s if not already
              on the v^c..y..w external face
            */

            ASSERT(embedg_VES_is_vertex(n, (*path_v)[0]));
            ASSERT(embedg_VES_is_vertex(n, s));

            embed_graph[s].visited = MARK_X_Y_PATH(n);
        }
        else  if (embed_graph[s].visited == MARK_X_Y_PATH(n))
            /*
              this means that s is a cut vertex on the internal
              face walk: pop all the vertices from path
              until s's last occurrence in path
            */
        {
            ASSERT((*nbr_v) >= 0);
            while ((*path_v)[(*nbr_v)] != s)
            {
                (*nbr_v)--;
                ASSERT((*nbr_v) >= 0);
                /*
                  note that s should be somewhere in path!
                */
            }
            /*
              note also that popping from path_v also implies
              popping from path_e
            */
            e = s_e;
            ein = s_ein;
        }
        else
            /*
              we push s and s_e on their respective stacks
            */
        {
            (*path_v)[++(*nbr_v)] = s;
            (*path_e)[*nbr_v] = s_e;
            e = s_e;
            ein = s_ein;
            
            if (embed_graph[s].visited == mark_r
                || embed_graph[s].visited == mark)
                /*
                  again, s lies on the v^c...y...w external face:
                  we end the walk: path_v now contains the highest x-y path

                  note that there can be no conflict between
                  mark_r or mark and MARK_X_Y_PATH(n) since
                  we mark with MARK_X_Y_PATH iff the vertex
                  is NOT marked with mark_r/mark!
                */
            {
                break;
            }
            else
                /*
                  we must mark this vertex as MARK_X_Y_PATH since we aren't
                  finished yet
                */
            {
                embed_graph[s].visited = MARK_X_Y_PATH(n);
            }
        }
    }

    /*
      there is only one thing remaining to do: see if p_x or
      p_y are attached high
      (ie closer to v^c than x or y resp.)

      we walk the external face starting at v^c in y's direction
      (again see embedg_iso_get_x_y_w)
    */
    *px_attached_high = TRUE;
    p_x = (*path_v)[0];
    /*
      p_y denotes the attachement point of the x-y path
      on the v^c...y...w external face
    */

    s = n;
    cur_v = vv;
    cur_vin = 0;
    while (s != p_x)
    {
        embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                             cur_v, cur_vin,
                                             FALSE, 0, &s, &sin);
        if (s == x)
        {
            *px_attached_high = FALSE;
            break;
        }
        cur_v = s;
        cur_vin = sin;
    }

    *py_attached_high = TRUE;
    p_y = (*path_v)[*nbr_v];
    /*
      p_y denotes the attachement point of the x-y path
      on the v^c...y...w external face
    */

    s = n;
    cur_v = vv;
    cur_vin = 1;
    while (s != p_y)
    {
        embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                             cur_v, cur_vin,
                                             FALSE, 0, &s, &sin);
        if (s == y)
        {
            *py_attached_high = FALSE;
            break;
        }
        cur_v = s;
        cur_vin = sin;
    }

    /*
      now we are in the minor C case if either p_x or p_y are
      attached high

      the minor D case:
      this happens when there is a path v^c - z where z lies
      on the x-y path

      that is, when
      
      either v^c has been effectively "avoided" within the
      embedg_VES_get_succ_on_proper_face_with_avoidance function
      BUT ONLY if this "avoidance" happened AFTER having
      encountered the very first vertex on the x-y path!

      or when a cut vertex has been encountered on the x-y path:
      separable components on this walk can only occur
      if one walks the face while skipping the edges incident to v^c

      in any case this means that checking the return from 
      the embedg_VES_get_succ_on_proper_face_with_avoidance function
      is enough: this is the purpose of avoid_vv.
    */

    *is_minor_D = !(*px_attached_high || *py_attached_high) && avoid_vv;


    IF_DEB(
           int    i;
           
           fprintf(stdout, "x-y path\t");
           for (i = 0; i <= *nbr_v; i++)
               fprintf(stdout, "%d\t", (*path_v)[i]);
           fprintf(stdout, "\n");
           )
}


/*
 *  embedg_misc.c
 */
 
/*
  What:
  *****
  
  Implementing:

  Some high level routinse for the VES structure.
  See VES_misc.c.


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
#define IF_DEB_TREE(x)    {}
 
 
 
/* aproto: header embed_graph_protos.h */
 
#ifndef PLANAR_IN_MAGMA
#endif
 



void 
embedg_VES_delete (t_ver_edge *embed_graph, int n)
{
    int          i;

    for (i = 0; i < n; i++)
    {
        embedg_dlcl_delete(embed_graph[i].separated_DFS_child_list);
        /*
          embedg_dlcl_delete(embed_graph[i].rep_in_parent_list);

          NO!!! this points to something in separated_DFS_child_list
        */
        embedg_dlcl_delete(embed_graph[i].pertinent_bicomp_list);
    }
    mem_free(embed_graph);
}    



void 
embedg_VES_print (t_ver_edge *embed_graph, int n)
{
    int          i;

    fprintf(stdout, "vertices\n");
    for (i = 0; i < n; i++)
    {
        t_ver_edge   rec;

        rec = embed_graph[i];

        fprintf(stdout, "\nDFI\t%d\tlabel\t%d\n", i, rec.label);
        fprintf(stdout, "DFS parent\t%d\tleast_a\t%d\tlowpoint\t%d\n",
                rec.DFS_parent, rec.least_ancestor, rec.lowpoint);
        fprintf(stdout, "separated_DFS_child_list\n");
        embedg_dlcl_print(rec.separated_DFS_child_list);
    }

    fprintf(stdout, "\nvirtual vertices\n");
    for (i = n; i < 2*n; i++)
    {
        int          c;
        
        c = i - n;
        fprintf(stdout, "%d^%d\t", embed_graph[c].DFS_parent, c);
    }
    fprintf(stdout, "\n");
    
    embedg_VES_print_bigcomps(embed_graph, n);
}


void 
embedg_VES_print_bigcomps (t_ver_edge *embed_graph, int n)
    /*
      walking the external faces of all the bicomp; for testing only
    */
{
    int          i;
 
    fprintf(stdout, "bicomponents\n");
    /*
      to get to the bicomps, it makes sense to start at the
      virtual vertices????
    */
    for (i = n + 1; i < 2*n; i++)
        /*
          a note of caution: there is no virtual vertex at
          embed_graph[n] since that would mean a virtual vertex x^0
          which makes no sense (0 is the root of the dfs_tree)
        */
    {
        embedg_VES_walk_bicomp(embed_graph, n, i, 0);
    }
    fprintf(stdout, "\n");
}
/*
 *  planar_alg_init.c
 */
 
/*
  What:
  *****
  
  Implementing:

  Initialising the embed_graph aka VES data structure from the information
  collected from the DFS.

  The embed_graph/VES data structure is an array consisting of vertices,
  virtual vertices and edges;
  vertices, virtual vertices and edges share a common record structure;
  one of the particular features is that any vertex is linked
  together with its incident edges into a doubly circular linked list.

  See also VES_misc.c.

  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/
 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_DEB_DFS(x) {}
#define IF_VERB(x)   {}
#define IF_DEB_TREE(x)    {}
#define IF_CPU(x) {}
 
 
/* aproto: header embed_graph_protos.h */

/* aproto: beginstatic -- don't touch this!! */
static void embedg_init_insert_TE (t_ver_edge *, int, int *, t_dlcl *);
/* aproto: endstatic -- don't touch this!! */
 
#ifndef PLANAR_IN_MAGMA
#endif


t_ver_edge *
embedg_planar_alg_init (
    t_ver_sparse_rep *V,
    int n,
    t_adjl_sparse_rep *A,        /* input sparse graph */
    int *nbr_c,        /* size of the graph, #components*/
    int *edge_pos,        /* pos in the struct where the last edge
                                      has been inserted
                                   */
    t_dlcl ***dfs_tree,      /* a sparse graph rep. for the dfs tree
                                      -- vertices are as DFIs
                                   */
    t_dlcl ***back_edges,    /* for each vertex v, a dlcl
                                      of the back edges [v, x] incident to v
                                      where x is a DESCENDANT of v
                                      -- vertices are as DFIs
                                   */
    t_dlcl ***mult_edges    /* for each vertex v, a dlcl
                                      of the back edges [v, x] incident to v
                                      where x is a DESCENDANT of v
                                      -- vertices are as DFIs
                                   */
)
    /*
      initialising embed_graph, the fundamental data structure
      underpinning the tester and obstruction isolator

      from there on, a vertex is exclusively referred to by its DFI!!
      -- so forget about labels
    */
{
    int          *dfs_nbr;        /* dfs numbering for each vertex */
    int          *dfs_order;      /* vertices in dfs order */
    int          *lowpoint;       /* lowpoint value for each DFI */
    int          *dfs_parent;     /* for each DFI, its DFS ancestor 
                                     as a DFI (DFS index)
                                  */
    int          *least_a;        /* for each DFI, its least ancestor's DFI
                                     (via a back edge exclusively)
                                  */

    t_ver_edge   *embed_graph;
    int          i;


  IF_CPU(
    float      sttime;  float time_to_now;
 
    sttime = time_current_user();
       )

    ASSERT(n >= 1);

    /*
      DFS and lowpoint calculations + ordering 
    */
    sparseg_adjl_dfs_preprocessing(V, n, A, nbr_c,
                                   &dfs_nbr, &dfs_order, &lowpoint,
                                   dfs_tree, back_edges,
                                   &dfs_parent, &least_a, mult_edges);

    IF_CPU(
           fprintf(stdout, "CPU for DFS only %f\n",
                   (time_current_user() - sttime));
    sttime = time_current_user();
           )
 
    IF_DEB_DFS(     
            fprintf(stdout, "DFS indices\n");
            for (i = 0; i < n; i++)
            fprintf(stdout, "%d ", dfs_nbr[i]);
            fprintf(stdout, "\n");
            
            fprintf(stdout, "DFS order\n");
            for (i = 0; i < n; i++)
            fprintf(stdout, "%d ", dfs_order[i]);
            fprintf(stdout, "\n");
            
            fprintf(stdout, "lowpoint values\n");
            for (i = 0; i < n; i++)
            fprintf(stdout, "%d ", lowpoint[i]);
            fprintf(stdout, "\n");
            );
    
    IF_VERB(
            fprintf(stdout, "DFS parent\n");
            for (i = 0; i < n; i++)
                fprintf(stdout, "%d ", dfs_parent[i]);
            fprintf(stdout, "\n");
            );

    IF_VERB(
            fprintf(stdout, "least ancestors\n");
            for (i = 0; i < n; i++)
                fprintf(stdout, "%d ", least_a[i]);
            fprintf(stdout, "\n");
            );
    
    IF_VERB(
            for (i = 0; i < n; i++)
            {
                fprintf(stdout, "the list of children ordered by lowpoint for %d\n",
                        i);
                embedg_dlcl_print((*dfs_tree)[i]);
            }
            );

    IF_DEB_DFS(
            fprintf(stdout, "the tree edges\n");
            sparseg_dlcl_print(*dfs_tree, n);
            
            fprintf(stdout, "the back edges\n");
            sparseg_dlcl_print(*back_edges, n);
 
            fprintf(stdout, "multiple edges\n");
            sparseg_dlcl_print(*mult_edges, n);
            );
    
    /*
      create the data structure for the embedded graph:
      it will have (max) size 2*n + 2 * MAXE(n)

      we will see that that number of edges is sufficient
      even when later adding short-cut edges (see embedg_walkdown)
    */
    embed_graph = (t_ver_edge *) mem_malloc(sizeof(t_ver_edge)
                                            * (2*n + 2 * MAXE(n)));
    /*
      initialisation
    */
    for (i = 0; i < 2*n + 2 * MAXE(n); i++)
        /*
          some fields are initialised to n as n is actually
          an "invalid" value 
        */
    {
        t_ver_edge   rec;
  
        rec.label = NIL;
        rec.DFS_parent = n;
        rec.least_ancestor = n;
        rec.lowpoint = n;
        rec.separated_DFS_child_list = NP;
        rec.rep_in_parent_list = NP;
        rec.pertinent_bicomp_list = NP;
        rec.adjacent_to = n;
        rec.visited = n;
        rec.neighbour = n;
        rec.in_adjl = NIL;
        rec.twin_in_adjl = NIL;
        rec.mult = 0;
        rec.type = NIL;
        rec.sign = NILSIGN;
        /*
          make the links refer back to self
        */
        rec.link[0] = rec.link[1] = i;
 
        embed_graph[i] = rec;
    }
    
    /*
      embed_graph[0..n-1]: the n vertices
      ATTENTION: the vertices are stored according to their DFS numbering
    */
    for (i = 0; i < n; i++)
    {
        t_ver_edge   rec;

        rec = embed_graph[i];

        rec.label = dfs_order[i];
        rec.DFS_parent = dfs_parent[i];
        rec.least_ancestor = least_a[i];  
        rec.lowpoint = lowpoint[i];
        rec.separated_DFS_child_list = embedg_dlcl_copy((*dfs_tree)[i]);

        IF_VERB(
                fprintf(stdout, "the list of children ordered by lowpoint for DFI %d\n",
                        i);
                embedg_dlcl_print(rec.separated_DFS_child_list);
        );

        embed_graph[i] = rec;
    }
    
    /*
      one more thing to do for these vertices:
      fix the rep_in_parent_list field
    */
    for (i = 1; i < n; i++)
    {
        t_dlcl       *parent_list, *rep;
        int          parent;

        parent = embed_graph[i].DFS_parent; /* careful: this is a DFI  */
        /*
          recall that the vertices in embed_graph are accessed via their DFI
        */

        if (parent != n)
            /*
              when parent == n this means that i the root of a DFS tree
              in the disconnected graph
            */
        {
            parent_list = embed_graph[parent].separated_DFS_child_list;
            rep = embedg_dlcl_find(parent_list, i);
            ASSERT(rep != NP);
            embed_graph[i].rep_in_parent_list = rep;
        }
    }

    /*
      embed_graph[n..2*n-1]: the n virtual vertices
      do I need to do anything here?????

      no - I don't think so

      let's try to explain what virtual vertices are:
      let v^c be a virtual vertex:
      - it is at position  c + n in the array,
      - c is the DFS child of v,
      - v can be retrieved by taking embed_graph[c].DFS_parent,
      - v^c is said virtual as long as the bicomp rooted by v^c is not
        merged with the vertex v
      - once v is merged (identified?) with v^c, then v^c
        is of no relevance anymore

      below we will see that we embed all the tree edges as singleton
      bicomps (bicomponent): (0^1, 1), (1^2, 2) etc...:
      this is what virtual vertices are there for:
      to distinguish them from their "real" counterpart with
      which they will be ultimately merged

      the primary reason for this is:
      while testing for planarity virtual vertices are the roots of bicomps
    */

    /*
      now the edges:
      we actually embed the tree edges so that each tree edge
      forms a (singleton) biconnected component

      embedding an edge in effect means creating the
      doubly linked circular list of [virtual] vertices & the edges incident
      to it

      this list is built using the links 0 & 1 in embed_graph[i]
    */

    /*
      for each tree edge (v,u) we embed (v^u, u) (v^u is the virtual vertex)

      CAREFUL: when talking about vertex v say, 
      we mean the vertex with DFI v, and NOT the vertex with label v
      **************************************************************
    */
    *edge_pos = 2*n - 1;
    /*
      edge_pos will tell us where to insert the next edge in embed_graph[]
    */
    for (i = 0; i < n; i++)
    {
        t_dlcl     *te_l, *p;

        te_l = (*dfs_tree)[i];
        p = te_l;
 
        if (!embedg_dlcl_is_empty(p))
        {
            /*
              the test below is a bit stupid... well...
            */
            ASSERT(embed_graph[p->info].DFS_parent == i);

            embedg_init_insert_TE(embed_graph, n, edge_pos, p);
            p = embedg_dlcl_list_next(p);
            while (p != te_l)
            {
                ASSERT(embed_graph[p->info].DFS_parent == i);
                embedg_init_insert_TE(embed_graph, n, edge_pos, p);

                p = embedg_dlcl_list_next(p);
            }
        }
    }
    
    mem_free(dfs_nbr);
    mem_free(dfs_order);
    mem_free(lowpoint);

    mem_free(dfs_parent);
    mem_free(least_a);
    
    IF_CPU(
           fprintf(stdout, "CPU for remainder of initialisation %f\n",
                   (time_current_user() - sttime));
           )

    return embed_graph;
}


static void
embedg_init_insert_TE (t_ver_edge *embed_graph, int n, int *edge_pos, t_dlcl *p)
    /*
      init and insert a tree edge in embed graph:
      
      the tree edge will form a singleton bicomponent (v^c, c)
      where c is p->info and v is c.DFS_parent
    */
{
    int             c, v;

    c = p->info;
    v = embed_graph[c].DFS_parent;
    ASSERT(v >= 0 && v < n);
    
    /*
      now (v, c) is a tree edge; embed the directed edge [v^c, c]
      
      -- and recall that v^c is a virtual vertex, at position  c + n
      in embed_graph, and that vertex c is at position c
    */
    
    /*
      first, set this edge with the appropriate info
    */
    (*edge_pos)++;
    ASSERT(*edge_pos < 2*n + 2 * MAXE(n));
    embed_graph[*edge_pos].neighbour = c;
    embed_graph[*edge_pos].in_adjl = p->in_adjl;
    embed_graph[*edge_pos].twin_in_adjl = p->twin_in_adjl;

    ASSERT(p->mult % 2 == 0);
    /*
      we want the number of undirected edges
    */
    embed_graph[*edge_pos].mult = p->mult / 2;
    embed_graph[*edge_pos].type = TE;
    embed_graph[*edge_pos].sign = CCLOCKW;
        
    /*
      link this with vertex v^c in a doubly linked circular list
    */
    embed_graph[c + n].link[0] =
        embed_graph[c + n].link[1] = *edge_pos;
    embed_graph[*edge_pos].link[0] =
        embed_graph[*edge_pos].link[1] = c + n;
        
    /*
      now create/set the twin edge, the directed edge [c, v^c] 
    */
    (*edge_pos)++;
    ASSERT(*edge_pos < 2*n + 2 * MAXE(n));
    embed_graph[*edge_pos].neighbour = c + n;
    embed_graph[*edge_pos].in_adjl = p->twin_in_adjl;
    embed_graph[*edge_pos].twin_in_adjl = p->in_adjl;
    embed_graph[*edge_pos].mult = p->mult / 2;
    embed_graph[*edge_pos].type = TE;
    embed_graph[*edge_pos].sign = CCLOCKW;
    
    /*
      and link it with vertex c in a doubly linked circular list
    */
    embed_graph[c].link[0] = embed_graph[c].link[1] = *edge_pos;
    embed_graph[*edge_pos].link[0] =
        embed_graph[*edge_pos].link[1] = c;
}
/*
 *  dfs_preprocessing.c
 */
 
/*
  What:
  *****
  
  Implementing:

  A DFS as an initialisation step for the planarity tester.
  This is an especially beefed up DFS that collects lots of
  marginal information:

  - a DFS tree as a list of DFS children for each vertex
  - the DFS children are sorted according to their lowpoint value
  - a back_edge structure as a list of descendants v for each
    vertex u such that [v, u] is a back edge
  - a multiple edges structure which stores multiple (directed) edges
    NOT in the DFS tree nor in the back_edge struc, and loops
    
  - the vertices in DFS order
  - the DFS index (DFI) for each vertex
  - the lowpoint value for each vertex
  - the number of components of the (possibly disconnected) graph
  - for each vertex, its DFS parent
  - for each vertex v, its least ancestor u such that [v, u]
    is a back edge

  ALL info above (except the vertices in DFS order) is given
  in terms of the vertices' DFIs and NOT their labels.
  

    
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 
 
  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca
 
 
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********
 
  Paulette Lieby (Magma), Brendan McKay (ANU)
 
  Started October 2001
*/


/*
  There are some dodgy things which need some thought; it would be nice
  to fix them so that the code get's cleaner:

  - we do store in_adj (and twin_in_adjl) for each directed edge:
  
    + this is ONLY needed at the time of recovering the embedding
      (see embedg_recover_embedding) and its sole use is to establish
      a link between an edge in the embedding structure
      t_embed_sparse_rep  *E and the corresponding edge
      in the t_adjl_sparse_rep   *A struc.

    + well, I cannot recall why I thought this correspondence
      was needed in the first place and it might well be the case
      that there is no use for it; in which case recovering the
      embedding is simplified
      (we would store the end-vertex in the embedding's edges instead
      of their index in the adjacency list)

  - there are some non-linear bits in the DFS below: when searching
    for an already existing tree/back/multiple edge.
    I couldn't fix this in less then one hour so I leave it as it is...
    for now.

    This shouldn't be a major issue, overall timings of the planarity
    tester do not show this non-linear "bump"...

  - also, this algorithm has been growing incrementally and I now
    realise that I am using some redundant data structures:
    for example visited[] and the vertex and could be dispensed with...
    ...more things to clean up...

    Paulette  07/02/02
*/
 

#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
#define IF_DEB_TREE(x)    {}
 
 
/* aproto: header embed_graph_protos.h */
 
#ifndef PLANAR_IN_MAGMA
#endif
 

void 
sparseg_adjl_dfs_preprocessing (
    t_ver_sparse_rep *V,
    int n,                /* size of the graph */
    t_adjl_sparse_rep *A,        /* input sparse graph */
    int *c,               /* nbr of components */
    int **dfs_nbr,        /* dfs numbering for each vertex */
    int **dfs_order,      /* vertices in dfs order */
    int **lowpoint,       /* lowpoint value for each DFI */
    t_dlcl ***dfs_tree,      /* a sparse graph rep. for the dfs tree:
                                      for each DFI, a list of its children's
                                      DFI ordered wrt their lowpoint values
                                   */
    t_dlcl ***back_edges,    /* for each DFI  v, a dlcl
                                      of the back edges [v, x] incident to v
                                      where x is a DESCENDANT of v */
    int **dfs_parent,     /* for each DFI its DFS ancestor */
    int **least_a,        /* for each DFI, its least ancestor's DFI
                                      via a back edge exclusively */
    t_dlcl ***mult_edges    /*  for each DFI  v, a dlcl
                                       of the multiple directed
                                       edges NOT included
                                       in either dfs_tree or back_edges
                                   */
)

    /*
      in ALL the returned info above BUT dfs_order[] we store
      the vertices' DFIs (DFS indices) and NOT their labels!

      -- shuffling between labels and vertices can then be done
         via dfs_nbr[] and dfs_order[]
    */
{
    int            pos_v_stack, pos_e_stack, dfs_n;
    int            *visited, *vertex_stack, *edge_stack, *lowpoint_order;
    int            *TE_in_adjl, *TE_twin_in_adjl, *TE_mult;
    int            v, lp, cur, cur_e, next;
    t_dlcl         **temp, *lowpoint_list, **new_dfs_tree;
    
    /*
      create the dfs tree as a sparse graph
    */
    *dfs_tree = (t_dlcl **) mem_malloc(sizeof(t_dlcl *) * n);
    /*
      the DFS numbering for the vertices
    */
    *dfs_nbr = (int *) mem_malloc(sizeof(int) * n);
    /*
      the vertices as ordered by their DFS index
    */
    *dfs_order = (int *) mem_malloc(sizeof(int) * n);
    /*
      the lowpoint value for each DFI
    */
    *lowpoint = (int *) mem_malloc(sizeof(int) * n);

    /*
      the (directed) back edges
    */
    *back_edges = (t_dlcl **) mem_malloc(sizeof(t_dlcl *) * n);


    /*
      the DFS parent for each DFI
    */
    *dfs_parent = (int *) mem_malloc(sizeof(int) * n);
    /*
      the least ancestor (via a back edge exlusively) for each DFI
    */
    *least_a = (int *) mem_malloc(sizeof(int) * n);
    
    /*
      the (directed) multiple edges
    */
    *mult_edges = (t_dlcl **) mem_malloc(sizeof(t_dlcl *) * n);

    /*
      the vertices visited while DFS
    */
    visited = (int *) mem_malloc(sizeof(int) * n);
    /*
      stack of vertices: last current vertex
    */
    vertex_stack = (int *) mem_malloc(sizeof(int) * n);
    /*
      stack of (tree) edges: last added tree edge
    */
    edge_stack = (int *) mem_malloc(sizeof(int) * n);
    
    /*
      the following will be used in order to recreate the dfs_tree
      so that the DFS children of each DFI are ordered
      according to their lowpoint value
    */
    lowpoint_order = (int *) mem_malloc(sizeof(int) * n);
    temp = (t_dlcl **) mem_malloc(sizeof(t_dlcl *) * n);
    new_dfs_tree = (t_dlcl **) mem_malloc(sizeof(t_dlcl *) * n);

    /*
      finally, three more holding arrays: a trick to remember which
      tree edges we are talking about:

      when constructing dfs_tree, back_edges, mult_edges
      - we NEED to record the index in A (the adjacency list)
        of some of the edges and their twins/inverses
        we are currently storing in either of these structures
      - we also need to record the number of multiple (directed)
        edges we encounter  when the graph is not simple
      
      this is easy to do when storing back edges and multiple edges,
      and tree edges also: but this lattest set of neighbour lists (dfs_tree)
      is subsequently reordered so that DFS children are ordered
      wrt lowpoint values;
      - consequently the info about position in adjacency list
      and edge multiplicity are lost in the ordering process

      the two following arrays will remember the info we'll need later
      - more about this below
    */
    TE_in_adjl = (int *) mem_malloc(sizeof(int) * n);
    TE_twin_in_adjl = (int *) mem_malloc(sizeof(int) * n);
    TE_mult = (int *) mem_malloc(sizeof(int) * n);
    
    
    /*
      initialization of the data structures
    */
    for (v = 0; v < n; v++)
    {
        (*dfs_tree)[v] = (*back_edges)[v] = (*mult_edges)[v] = NP;
        visited[v] = TE_mult[v] = 0;
        (*dfs_parent)[v] = (*least_a)[v] = n;
        temp[v] = new_dfs_tree[v] = NP;
        TE_in_adjl[v] = TE_twin_in_adjl[v] = NIL;
        /*
          note that in the 3rd last statement n is considered
          as an "invalid" value;
          will be if importance in the overall algorithm
        */
    }
    
    /*
      the DFS tree is rooted at vertex 0
    */
    dfs_n = -1;
    pos_v_stack = -1;
    pos_e_stack = -1;
    *c = 0;
    for (v = 0; v < n; v++)
    {
        if (visited[v])
            /*
              we come only at this level when looking for
              a new subtree (when graph is disconnected)
            */
        {
            continue;
        }
        else
        {
            (*c)++;
        }
        
        cur = v;
        visited[cur] = 1;
        (*dfs_nbr)[cur] = ++dfs_n;
        (*lowpoint)[(*dfs_nbr)[cur]] = dfs_n;
        (*dfs_order)[dfs_n] = cur;

        cur_e = V[cur].first_edge == NIL ? NIL : V[cur].first_edge;
        while (TRUE)
        {
            if (cur_e != NIL)
            {
                t_dlcl        *existing_e;
                
                next = A[cur_e].end_vertex;
                if (!visited[next])
                    /*
                      adding tree edges (careful: directed edges)
                      
                      AND tree edges are stored as 
                      [dfs_nbr[u], dfs_nbr[cv]]
                      instead of [u, cv]: that is we store the edges
                      according to the vertices' DFIs
                    */
                {
                    IF_DEB_TREE(
                                io_printf("add tree edge %d\t%d\n",
                                          cur+1, next+1);
                                );
                    
                    (*dfs_nbr)[next] = ++dfs_n;
                    (*lowpoint)[(*dfs_nbr)[next]] = dfs_n;
                    (*dfs_order)[dfs_n] = next;
                    
                    sparseg_dlcl_append_to_neigh_list(*dfs_tree, n,
                                                      (*dfs_nbr)[cur],
                                                      (*dfs_nbr)[next],
                                                      NIL);
                    TE_in_adjl[(*dfs_nbr)[next]] = cur_e;
                    TE_mult[(*dfs_nbr)[next]]++;
                    
                    /*
                      we push cur and the edge (cur, cur_e) on their
                      respective stacks
                    */
                    vertex_stack[++pos_v_stack] = cur;
                    edge_stack[++pos_e_stack] = cur_e;
                    
                    /*
                      and mark next as visited
                    */
                    visited[next] = 1;

                    /*
                      update dfs_parent (always deal with DFIs rembember!)
                    */
                    (*dfs_parent)[(*dfs_nbr)[next]] = (*dfs_nbr)[cur];
                    
                    /*
                      the DFS goes one level deeper
                    */
                    cur = next;
                    cur_e = V[cur].first_edge == NIL ?
                        NIL : V[cur].first_edge;
                }
                /*
                  the next three tests deal with multiple edges
                  and loops: apart from storing these (DIRECTED) edges
                  in mult_edges, we also need to update
                  the multipliciaty information about these edges
                */
                else if (sparseg_dlcl_is_adjacent(*dfs_tree, n,
                                                  (*dfs_nbr)[cur],
                                                  (*dfs_nbr)[next],
                                                  &existing_e))
                    /*
                      [cur, next] is a tree edge
                    */
                {
                    sparseg_dlcl_append_to_neigh_list(*mult_edges, n,
                                                      (*dfs_nbr)[cur],
                                                      (*dfs_nbr)[next],
                                                      cur_e);
                    TE_mult[(*dfs_nbr)[next]]++;
                    
                    cur_e = A[cur_e].next; /* next in cur's adjacency list */
                }
                else if (sparseg_dlcl_is_adjacent(*back_edges, n,
                                                  (*dfs_nbr)[next],
                                                  (*dfs_nbr)[cur],
                                                  &existing_e))
                    /*
                      [cur, next] is a back edge
                    */
                {
                    sparseg_dlcl_append_to_neigh_list(*mult_edges, n,
                                                      (*dfs_nbr)[cur],
                                                      (*dfs_nbr)[next],
                                                      cur_e);
                    (existing_e->mult)++;
                    
                    cur_e = A[cur_e].next; /* next in cur's adjacency list */
                }
                else if (next == cur)
                    /*
                      the case of a loop
                    */
                {
                    if (sparseg_dlcl_is_adjacent(*mult_edges, n,
                                                  (*dfs_nbr)[next],
                                                  (*dfs_nbr)[cur],
                                                  &existing_e))
                        /*
                          in this case we must update the multiplicity
                          of this edge: note that the elt. in cur's
                          neighbours list that gets updated is the first
                          in the list

                          dodgy??? certainly, but can't think
                          of a better way to do this

                          eventually it will happen that even myself
                          won't understand what I am doing..........
                        */
                    {
                        (existing_e->mult)++;
                    }
                    sparseg_dlcl_append_to_neigh_list(*mult_edges, n,
                                                      (*dfs_nbr)[cur],
                                                      (*dfs_nbr)[next],
                                                      cur_e);
                    
                    cur_e = A[cur_e].next; /* next in cur's adjacency list */
                }
                else if (sparseg_dlcl_is_adjacent(*dfs_tree, n,
                                                  (*dfs_nbr)[next],
                                                  (*dfs_nbr)[cur],
                                                  &existing_e))
                    /*
                      [next, cur] is a tree edge:
                      that is, [cur, next] is [next, cur]'s twin/inverse:

                      1. if it is the first time one encounters
                         [cur, next] (as it would always be the case
                         for a simple graph) then all I need to do
                         is to update the tree edge's multiplicity,
                         and the twin info in TE_[]

                      2. if [cur, next] is actually a multiple edge,
                         then I'll need to store it in mult_edges;
                         and I update the tree edge's multiplicity too.
                         No twin info will be required here.
                         Why? see how recover.c embeds the multiple
                         edges in the planar embedding.

                      3. how do I know it is the first time I encounter
                         [cur, next]?:
                         when TE_twin_in_adjl = NIL

                      4. finally, note that the present counting scheme
                         implies that the mult field always holds
                         the number of directed edges:
                         ie, if [a, b] is a tree edge, [a, b].mult = 2
                         because we would have counted [a, b] and [b, a]

                         this applies to tree edges, back edges, and loops
                    */
                {
                    ASSERT(TE_in_adjl[(*dfs_nbr)[cur]] != NIL);
                    if (TE_twin_in_adjl[(*dfs_nbr)[cur]] == NIL)
                    {
                        TE_twin_in_adjl[(*dfs_nbr)[cur]] = cur_e;
                    }
                    else
                    {
                        sparseg_dlcl_append_to_neigh_list(*mult_edges, n,
                                                          (*dfs_nbr)[cur],
                                                          (*dfs_nbr)[next],
                                                          cur_e);
                    }
                     
                    TE_mult[(*dfs_nbr)[cur]]++;

                    cur_e = A[cur_e].next; /* next in cur's adjacency list */
                }
                else if (sparseg_dlcl_is_adjacent(*back_edges, n,
                                                  (*dfs_nbr)[cur],
                                                  (*dfs_nbr)[next],
                                                  &existing_e))
                    /*
                      [next, cur] is a back edge: [cur, next] is its inverse:
                      we proceed as for the tree edge case above
                    */
                {
                    ASSERT(existing_e->in_adjl != NIL);
                    if (existing_e->twin_in_adjl == NIL)
                    {
                        existing_e->twin_in_adjl = cur_e;
                    }
                    else
                    {
                        sparseg_dlcl_append_to_neigh_list(*mult_edges, n,
                                                          (*dfs_nbr)[cur],
                                                          (*dfs_nbr)[next],
                                                          cur_e);
                    }
                    
                    (existing_e->mult)++;
                    
                    cur_e = A[cur_e].next; /* next in cur's adjacency list */
                }
                /*
                  the next bit concludes the DFS: it deals with the case
                  where a back edge needs to be added 
                */
                else 
                    /*
                      that is, next is visited and neither
                      the tree edge [next, cur] nor
                      the back edge [next, cur] exist:
                      
                      this implies that [cur, next] is a back edge
                      that must be added to the back_edges structure
                      (with dfs_nbr(next) < dfs_nbr(cur))
                    */
                {
                    IF_DEB_TREE(
                                io_printf("add back edge %d\t%d\n",
                                          cur+1, next+1);
                                );
                    
                    ASSERT(visited[next]);
                    ASSERT((*dfs_nbr)[cur] > (*dfs_nbr)[next]);
                    
                    sparseg_dlcl_append_to_neigh_list(*back_edges, n,
                                                      (*dfs_nbr)[next],
                                                      (*dfs_nbr)[cur],
                                                      cur_e);

                    /*
                      update cur's lowpoint
                    */
                    (*lowpoint)[(*dfs_nbr)[cur]] =
                        (*dfs_nbr)[next] < (*lowpoint)[(*dfs_nbr)[cur]] ?
                        (*dfs_nbr)[next] : (*lowpoint)[(*dfs_nbr)[cur]];
                    
                    /*
                      update least_a (of cur)
                      (always deal with DFIs remember!)
                    */
                    (*least_a)[(*dfs_nbr)[cur]] =
                        (*dfs_nbr)[next] < (*least_a)[(*dfs_nbr)[cur]] ?
                            (*dfs_nbr)[next] : (*least_a)[(*dfs_nbr)[cur]];

                    /*
                      get the next edge in cur's adjacency list
                    */
                    cur_e = A[cur_e].next;
                }
            }
            
            if (cur_e == NIL)
                /*
                  we are either at a leaf or have finished scanning
                  cur's adjacency list: backtrack
                */
            {
                if (pos_v_stack == -1)      /* no previous vertex */
                {
                    /*
                      no edge left on the stack: DFS ends for
                      this subtree:
                      we visit the next vertex
                    */
                    ASSERT(pos_e_stack == -1);
                    break;
                }
                else
                {
                    int      prev_e;
                    /*
                      Otherwise backtrack and pop cur from the stack
                      as well as the last tree edge  added to the tree.
                      We use next to get a new lowpoint value for cur:
                      This value will be min(lowpoint(cur), lowpoint(next)).
                    */
                    cur = vertex_stack[pos_v_stack--];
                    prev_e = edge_stack[pos_e_stack--];
                    next = A[prev_e].end_vertex;
                    (*lowpoint)[(*dfs_nbr)[cur]] =
                        (*lowpoint)[(*dfs_nbr)[cur]]
                        < (*lowpoint)[(*dfs_nbr)[next]] ?
                        (*lowpoint)[(*dfs_nbr)[cur]]
                        : (*lowpoint)[(*dfs_nbr)[next]];

                    cur_e = A[prev_e].next;
                }
                /*
                  we proceed with DFS
                */
            }
        }
    }
    mem_free(vertex_stack);
    mem_free(edge_stack);

    /*
      just for the sake of it, check that all vertices have
      been visited
    */
#ifdef ASSERTIONS
    for (v = 0; v < n; v++)
    {
        ASSERT(visited[v]);
    }
#endif
    mem_free(visited);

    /*
      we now order the DFIs wrt lowpoint values:
      use bucket sort (linear time)
    */
    /*
      for each lowpoint value, collect the DFIs (in a t_dlcl)
      with that lowpoint value
      (IMPORTANT: we want the DFIs since the aim is to rewrite dfs_tree
      which stores DFIs and not labels!)
    */
    for (v = 0; v < n; v++)
        /*
          v  is taken as a DFI here
        */
    {
        t_dlcl    *r;
 
        r = embedg_dlcl_rec_new(v);
        temp[(*lowpoint)[v]] = 
            embedg_dlcl_rec_append(temp[(*lowpoint)[v]], r);
    }
 
    /*
      concatenate these lists now
    */
    lowpoint_list = temp[0];
    for (lp = 1; lp < n; lp++)
    {
        lowpoint_list = embedg_dlcl_cat(lowpoint_list, temp[lp]);
    }
    ASSERT(embedg_dlcl_length(lowpoint_list) == n);
        
    lowpoint_order[0] = lowpoint_list->info;
    for (lp = 1; lp < n; lp++)
    {
        lowpoint_list = embedg_dlcl_list_next(lowpoint_list);
        lowpoint_order[lp] = lowpoint_list->info;
    }
    embedg_dlcl_delete(lowpoint_list);
    mem_free(temp);

    IF_DEB(
           fprintf(stdout, "dfs_preprocessing, lowpoint_order\n");
           for (lp = 0; lp < n; lp++)
           fprintf(stdout, "%d ", lowpoint_order[lp]);
           fprintf(stdout, "\n");
           fprintf(stdout, "dfs_preprocessing, lowpoint\n");
           for (lp = 0; lp < n; lp++)
           fprintf(stdout, "%d ", (*lowpoint)[lp]);
           fprintf(stdout, "\n");
           )
    
    /*
      we now use this order to rewrite dfs_tree such that
      the DFS children of each vertex are ordered wrt lowpoint values
    */
    for (lp = 0; lp < n; lp ++)
        /*
          for each DFI in lowpoint_order[] I know its DFS_parent
          from dfs_parent[] -- the rest is then trivial
        */
    {
        int       parent;
        
        v = lowpoint_order[lp];
        /*
          lowpoint_order stores DFIs as does dfs_parent, so the lot
          makes sense
        */
        parent = (*dfs_parent)[v];
        if (parent != n)
            /*
              v may be the root of a DFS tree
            */
        {
            t_dlcl   *temp;

            temp = embedg_dlcl_rec_new(v);

            /*
              this is where the TE_ holding arrays are useful *sigh*
            */
            ASSERT(TE_in_adjl[v] != NIL);
            temp->in_adjl = TE_in_adjl[v];

            ASSERT(TE_twin_in_adjl[v] != NIL);
            temp->twin_in_adjl = TE_twin_in_adjl[v];

            ASSERT(TE_mult[v] != 0 && TE_mult[v] % 2 == 0);
            temp->mult = TE_mult[v];
            
            new_dfs_tree[parent] =
                embedg_dlcl_rec_append(new_dfs_tree[parent], temp);
        }
    }
    mem_free(lowpoint_order);
    mem_free(TE_in_adjl);
    mem_free(TE_twin_in_adjl);
    mem_free(TE_mult);

    /*
      some checks are in order here
    */
#ifdef ASSERTIONS
    for (v = 0; v < n; v++)
    {
        ASSERT(embedg_dlcl_length((*dfs_tree)[v])
               == embedg_dlcl_length(new_dfs_tree[v]));

        IF_DEB(
               fprintf(stdout, "dfs_preprocessing    dfs_tree for %d\n", v);
               embedg_dlcl_print((*dfs_tree)[v]);
               fprintf(stdout, "dfs_preprocessing    new_dfs_tree for %d\n", v);
               embedg_dlcl_print(new_dfs_tree[v]);
    );
    }
#endif
    
    sparseg_dlcl_delete(*dfs_tree, n);
    *dfs_tree = new_dfs_tree;
}

/*
 *  embedding.c
 */
 
/*
  What:
  *****
  
  Implementing:

  The graph is planar: we recover the embedding from the VES structure
  and check it as well.
  (Some of these checks will disappear later)


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_DEB_EMBED(x)    {}
#define IF_DEB_CHECK_EMBED(x)    {}
#define IF_DEB_FACES(x) {}
#define IF_VERB(x)   {}
#define IF_DEB_SCE(x) {}
#define IF_CPU(x) {}
 

/* aproto: header embed_graph_protos.h */
 
 
#ifndef PLANAR_IN_MAGMA
#endif
 
void 
embedg_embedding (t_ver_sparse_rep *V, t_adjl_sparse_rep *A,
	t_ver_edge *embed_graph, int n, int e, int nbr_c,
	int edge_pos, t_dlcl **mult_edges, t_ver_sparse_rep **vertices,
	t_embed_sparse_rep **embedding)
    /*
      recovering the embedding for the (planar) graph

      - the embedding is returned in vertices and embedding, vertices
      indexes embedding, the ordered list of edges
      - edges in the embedding are given as their index in A, the graph's
      adajacency list
      - the nbr of edges in the embedding is given as nbr_e_embed:
      this may be different form the original number of edges when the graph
      iss not simple
    */
{
    int          *ver_orient, nbr_comp, nbr_e_embed;
    
    IF_CPU(
    float      sttime; float time_to_now;
 
    sttime = time_current_user();
    )

    IF_DEB(
           fprintf(stdout, "embedding, begin, which edges have been flipped\n");
           embedg_VES_print_flipped_edges(embed_graph, n, edge_pos);
           )

    IF_DEB(
           fprintf(stdout, "embedding, before removing SCE\n");
           embedg_VES_print_bigcomps(embed_graph, n);
           )

    /*
      several things to do:
      1. removing the short-cut edges
    */
    embedg_remove_SCE(embed_graph, n, edge_pos);

    IF_DEB(
           fprintf(stdout, "embedding, after removing SCE\n");
           embedg_VES_print_bigcomps(embed_graph, n);
           )
        
    /*
      2. computing each vertex's orientation (wrt flipped bicomps)
    */
    ver_orient = embedg_vertices_orientation(embed_graph, n);
    

    /*
      3. merging the remaining virtual vertices with their
         non-virtual counterpart
    */
    nbr_comp = embedg_merge_remaining_virtual(embed_graph, n);
    /*
      actually there is no need to return the nbr of components
      from the above function
      but let's do it for the sake of it and for possible checking
    */
    ASSERT(nbr_c == nbr_comp);

    IF_DEB(
           fprintf(stdout, "embedding, after merging of remaining vertices\n");
           )

    /*
      4. to be on the safe side: check that the embedding is a valid one

      for now, we DIE if not
    */

    if (!embedg_is_embed_valid(embed_graph, n, nbr_comp, edge_pos,
                               ver_orient, &nbr_e_embed))
    {
        mem_free(ver_orient);
        DIE();
    }
    mem_free(ver_orient);

    ASSERT(nbr_e_embed <= e);
    /*
      when the graph is not simple, multiple edges and loops are
      not in embed_graph[]: they will be added to the final
      embedding in embedg_recover_embedding below
    */
    
    /*
      5. recover the embedding in preparation for the Magma type,
      and check it as well
    */
    embedg_recover_embedding(V, A, embed_graph, n, e,
                             mult_edges, vertices, embedding);
    if (!embedg_check_recov_embedding(n, e, nbr_comp,
                                      *vertices, A, *embedding))
    {
        mem_free(*vertices);
        mem_free(*embedding);
        
        IF_CPU(
               fprintf(stdout, "CPU for embedding recovering %f\n",
                       time_current_user() - sttime);
               )
 
        DIE();
    }

    IF_DEB_EMBED(
                 fprintf(stdout, "embedding, original graph and embedding\n");
                 sparseg_adjl_print(V, n, A, FALSE);
                 fprintf(stdout, "\n");
                 sparseg_adjl_embed_print(*vertices, n, A, *embedding,
                                          FALSE);
                 )
    
    IF_CPU(
           fprintf(stdout, "CPU for embedding recovering %f\n",
                   time_current_user() - sttime);
           )
}


void 
embedg_remove_SCE (t_ver_edge *embed_graph, int n, int edge_pos)
    /*
      remove all the short-cut edges from the embedding
    */
{
    int          i, c;

    c = 0;
    for (i = 2*n; i <= edge_pos; i += 2)
        /*
          and edge and its twin occupy consecutive positions in embed_graph:
          need only to examine one out of two
          (removing an edge also entails removing its twin of course
        */
    {
        if (embedg_VES_is_short_cut_edge(embed_graph, n, i))
        {
            IF_DEB_SCE(
                       fprintf(stdout, "remove SCE\n");
                       embedg_VES_print_edge(embed_graph, n, i);
                       )

            embedg_VES_remove_edge(embed_graph, n, i);
            c++;
        }
    }

    IF_DEB_SCE(
               fprintf(stdout, "nbr of SCE edges removed %d\n", c);
               )
}


int *
embedg_vertices_orientation (t_ver_edge *embed_graph, int n)
    /*
      for each vertex return its orientation from the
      bicomps in embed_graph:
      perform a DFS of each bicomp
    */
{
    int          i, vv, prod_sign;
    int          *stack, *ver_orient, to_prev;

    /*
      the whole lot makes sense iff the adjacency lists are consistent:
      this is a very important issue and it might be the case
      that the ASSERT warrants replacement by a DIE
      (the check is linear - I think)
    */
    ASSERT(embedg_VES_are_adj_lists_consistent(embed_graph, n));

    ver_orient = (int *) mem_malloc(sizeof(int) * n);
    for (i = 0; i < n; i++)
    {
        ver_orient[i] = CCLOCKW;
    }

    /*
      create the stack for the DFS
    */
    stack = (int *) mem_malloc(sizeof(int) * 3*n);
    to_prev = -1;
    
    IF_DEB(
           fprintf(stdout, "vertex orientation, one line (of vert.) for each bicomp\n");
           )
                
    /*
      now visit all the bicomps, ie, all the virtual vertices
      in embed_graph
    */
    for (vv = n; vv < 2*n; vv++)
    {
        int      c, cur, cur_e;
        boolean  NEW_BICOMP;

        if (embed_graph[vv].link[0] == vv)
            /*
              means that vv is disabled and is not the root of a bicomp
            */
        {
            continue;
        }

        c = vv - n;
        IF_DEB(
               fprintf(stdout, "%d ", c);
               )
        /*
          orientation for c (vv is as yet unembedded) is CCLOCKW

          now find the orientation of all its DFS descendants
        */

        if (embed_graph[c].DFS_parent == n)
            /*
              this means that actually c is an isolated vertex:
              we initialise the sign to CCLOCKW
            */
        {
            prod_sign = CCLOCKW;
        }
        else
            /*
              we initialise the sign to CCLOCKW to the sign of c's parent
            */
        {
            prod_sign = ver_orient[embed_graph[c].DFS_parent];
        }

        /*
          we must not forget to set c's sign!!
          (won't be done below)
        */
        ver_orient[c] = prod_sign;
        
        NEW_BICOMP = FALSE;
        cur = c;
        cur_e = embed_graph[cur].link[0];
        ASSERT(embedg_VES_is_edge(n, cur_e));

        ASSERT(to_prev == -1);
        while (TRUE)
        {
            while (!embedg_VES_is_tree_edge(embed_graph, n, cur_e)
                   || !embedg_VES_is_vertex(n,
                                             embed_graph[cur_e].neighbour)
                   || embed_graph[cur_e].neighbour <= cur)
                /*
                  want to find a tree edge [cur, u]
                  where u is a descendant of cur
                */
            {
                cur_e = embed_graph[cur_e].link[0];
                
                while (cur_e == cur)
                    /*
                      back to the vertex where we started from:
                      no edge has been found:
                      cur is a leaf, backtrack
                    */
                {
                    if (to_prev == -1)
                    {
                        NEW_BICOMP = TRUE;
                        break;
                    }
                    prod_sign = stack[to_prev--];
                    cur_e = stack[to_prev--];
                    /*
                      must advance one more edge
                    */
                    cur_e = embed_graph[cur_e].link[0];
                    cur = stack[to_prev--];
                }
                if (NEW_BICOMP)
                {
                    break;
                }
            }

            if (NEW_BICOMP)
            {
                break;
            }
            else
                /*
                  now cur_e is the edge we were looking for, get its sign
                */
            {
                /*
                  push on stack the current vertex, the edge where we
                  stopped the DFS, AND the sign carried by that vertex
                  
                  and go down one level in the DFS
                */
                stack[++to_prev] = cur;
                stack[++to_prev] = cur_e;
                stack[++to_prev] = prod_sign;
                
                cur = embed_graph[cur_e].neighbour;
                prod_sign *= embed_graph[cur_e].sign;
                ver_orient[cur] = prod_sign;
                
                cur_e = embed_graph[cur].link[0];
                ASSERT(embedg_VES_is_edge(n, cur_e));
                
                IF_DEB(
                       fprintf(stdout, "%d with sign %d\n", cur, prod_sign);
                       )
            }
        }

        IF_DEB(
               fprintf(stdout, "\n");
               )
    }

    IF_DEB(
           fprintf(stdout, "vertex orientation\n");
           for (i = 0; i < n; i++)
           {
               fprintf(stdout, "%d ", ver_orient[i]);
           }
           fprintf(stdout, "\n");
           )

   mem_free(stack);
   return ver_orient;
}
    

int 
embedg_merge_remaining_virtual (t_ver_edge *embed_graph, int n)
    /*
      after the short-cut edges have been removed and the vertices'
      orientation computed, one finishes by merging all
      remaining virtual vertices with their virtual counterpart
      (without flip of course)

      and use this routine to return  the number of disconnected
      components of the graph
    */
{
    /*
      at this stage it is easy to see that all remaining
      virtual vertices are DFS roots (if the graph is not connected)
      or cut vertices
    */

    int          vv, nbr_comp;

    nbr_comp = 0;
    for (vv = n; vv < 2*n; vv++)
    {
        int      v, c;


        c = vv - n;
        v = embed_graph[c].DFS_parent;
        
        /*
          must fish out which virtual vertices are actual roots
          of DFS trees (esp. for the disconnected graph case):
          roots of DFS trees are those virtual vertices for which
          v = embed_graph[c].DFS_parent = n
        */
        if (v == n)
        {
            nbr_comp++;
            continue;
        }
        
        if (embed_graph[vv].link[0] == vv)
            /*
              means that vv is disabled and is not the root of a bicomp
            */
        {
            continue;
        }

        embedg_VES_merge_simple_bicomps(embed_graph, n,
                                                   vv, 1, v, 0);
        /*
          note:
          since v is a cut vertex in this intance the bicomp
          rooted by vv will be merged without flip;
          therefore we could have done
          embedg_VES_merge_simple_bicomps(embed_graph, n,
          vv, 0, v, 1)
          as well, the important thing being that vin != vvout
          (see embedg_VES_merge_simple_bicomps)
        */
    }

    return nbr_comp;
}


int 
embedg_nbr_faces (t_ver_edge *embed_graph, int n, int edge_pos,
	int *ver_orient, int *nbr_e_embed)
    /*
      count the number of faces and the number of edges of the embedding
    */
{
    int          v, e, f, total_e;

    IF_DEB_FACES(
                 int    v;
                 
                 fprintf(stdout, "nbr of faces, the vertices' adj. lists\n");
                 for (v = 0; v < n; v++)
                     embedg_VES_print_adj_list(embed_graph, n,
                                                          v, TRUE);
                 )
        
    /*
      the following is no more than a quick check -- certainly
      not very useful -- or could be done elsewhere
    */
    total_e = 0;
    for (e = 2*n; e <= edge_pos; e++)
    {
        if (!embedg_VES_is_short_cut_edge(embed_graph, n, e))
        {
            total_e++;
        }
    }
    ASSERT(total_e % 2 == 0);
    *nbr_e_embed = total_e / 2;

    /*
      I now set each edge's orientation

      QUESTION: do I really need to do this???
      so far, when doing a proper face traversal, the way in which
      the adjacency list of an edge must be traversed is given
      by the vertex's (in that list) orientation...
      So this seems sensible to me huh?
    */
    embedg_VES_set_orientation(embed_graph, n, ver_orient);

    /*
      I will be using the visited field to enable me to check
      if all edges have been traversed

      let's be smart (?!): so far the visited field has been used
      and set in the following circumstances:
      + initialisation: set to n
      + walkup: set to whatever DFI of interest

      so here we set it to MARK_EMBED(n)
    */
    f = 0;
    for (e = 2*n; e <= edge_pos; e++)
    {
        if (!embedg_VES_is_short_cut_edge(embed_graph, n, e)
            /*
              arrghh!!! I must also skip the SCE!!!
            */
            && embed_graph[e].visited != MARK_EMBED(n))
        {
            int          ein;
            
            IF_DEB_FACES(
                         fprintf(stdout, "nbr of faces, edges not visited\n");
                         embedg_VES_print_edge(embed_graph, n, e);
                         )

            ein = embed_graph[e].sign == CCLOCKW ? 0 : 1;
            /*
              the way I enter e in dependent on its sign:
              all the proper face traversal must obviously be done
              with the same orientation!
            */
            embedg_VES_walk_proper_face(embed_graph, n, e,
                                                       ein,
                                                       TRUE,
                                                       MARK_EMBED(n));
            f++;
        }
    }

    /*
      counting the faces by traversing all the edges does not
      account of the face defined by isolated vertices
      -- we do that now

      we only need to check which vertices refer to self, ie with
      no incident edges
    */
    for (v = 0; v < n; v++)
    {
        if (embed_graph[v].link[0] == v)
        {
            ASSERT(embed_graph[v].link[1] == v);
            f++;
        }
    }
    
    return f;
}


boolean 
embedg_is_embed_valid (t_ver_edge *embed_graph, int n, int nbr_comp,
	int edge_pos, int *ver_orient, int *nbr_e_embed)
    /*
      use Euler's formula to assertain that the embedding is a valid
      embedding:

      f = 2 * nbr_comp + nbr_e_embed - n

    */
{
    int         v, f;

    f = embedg_nbr_faces(embed_graph, n, edge_pos, ver_orient, nbr_e_embed);

    IF_DEB_CHECK_EMBED(
                       fprintf(stdout, "embedding, n: %d\t e: %d\t C: %d\t f: %d\n",
                               n, nbr_e, nbr_comp, f);
                       )
        
    return f == 2 * nbr_comp + *nbr_e_embed - n ? TRUE : FALSE;
}
/*
 *  ext_face_walk.c
 */
 
/*
  What:
  *****
  
  Implementing the external face walk of a bicomponent.
  The concept of an external face --in the context of the VES
  data structure-- makes only sense when talking
  about a bicomp.

  Recall that a vertex is linked together with the edges
  incident from it in a circular (doubly) linked list
  (this is the VES structure).

  One particular feature is that if a vertex v is on
  the external face of a component and if in the list
  we have edges e1, e2 such as e1  -> v -> e2
  then e1 and e2 border the external face.

  In other words, in the circular list of vertex v and edges,
  v is ALWAYS between the two edges bordering the external face

  Of course, when v is (maybe) pushed into the internal face
  (by embedding of some edge) then we don't care about this any more
  (for v that is).
  


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
 
 
/* aproto: header embed_graph_protos.h */
 

void 
embedg_VES_get_succ_on_ext_face (t_ver_edge *embed_graph, int n, int v,
	int vin, boolean MARK, int mark, int *s, int *sin)
    /*
      find the successor s of v (entered via vin) on the external face
      -- also return the direction in which s has been entered

      if MARK true mark the succ. vertex and the edges traversed
      with mark (the visited field)
    */
{
    int        e, twin;
    int        vout, ein, eout, tout;

    ASSERT(embedg_VES_is_vertex(n, v)
           || embedg_VES_is_virtual_vertex(n, v));

    IF_DEB(
           fprintf(stdout, "get_succ_on_ext_face, of %d:%d\n", v, vin);
           )
        
    /*
      find the direction out of the vertex, and get the edge
    */
    vout = vin == 0 ? 1 : 0;
    e = embed_graph[v].link[vout];
    if (embedg_VES_is_virtual_vertex(n, v) && e == v)
        /*
          this can happen if a virtual vertex has been "disabled"

          -- this should not never happen since we can only walk
          on the external face of a bicomp!
        */
    {
        *s = v;
        *sin = vin;
        return;
    }

    /*
      otherwise we must have an edge:
      note that it is entirely irrelevant if I walk SCEs:
      those are precisely there to "jump" over inactive vertices
    */
    ASSERT(embedg_VES_is_edge(n, e));

    /*
      get the twin edge
    */
    twin = embedg_VES_get_twin_edge(embed_graph, n, e);

    IF_DEB(
           fprintf(stdout, "get_succ_on_ext_face, edge [%d, %d]\n",
                   v, embed_graph[e].neighbour);
           fprintf(stdout, "get_succ_on_ext_face, twin edge [%d, %d]\n",
                   embed_graph[e].neighbour, embed_graph[twin].neighbour);
           )
    /*
      find which of twin's link links a vertex
    */
    tout = embedg_VES_is_vertex(n, embed_graph[twin].link[0])
        || embedg_VES_is_virtual_vertex(n,
                                                   embed_graph[twin].link[0])
        ?
        0 : 1;

    /*
      get this vertex: this is v's successor on the external face
    */
    *s = embed_graph[twin].link[tout];

    /*
      one more thing to do: find the direction in which s was entered
    */
    *sin = embed_graph[*s].link[0] == twin ? 0 : 1;

    IF_DEB(
           fprintf(stdout, "get_succ_on_ext_face, succ is %d:%d\n",
                   *s, *sin);
           )
    /*
      a special case: when the bicomp is a singleton bicomp
      (ie a single edge)
    */
    if (embed_graph[*s].link[0] == (embed_graph[*s].link[1]))
    {
        ASSERT(embed_graph[*s].link[0] = twin);
        *sin = vin;
    }

    /*
      finally, mark the vertex and edges if so requested
    */
    if (MARK)
    {
        embed_graph[*s].visited = mark;
        embed_graph[e].visited = mark;
        embed_graph[twin].visited = mark;
    }
}

void 
embedg_VES_get_succ_active_on_ext_face (t_ver_edge *embed_graph, int n,
	int v, int w, int win, boolean MARK, int mark, int *s, int *sin)
    /*
      find the ACTIVE (wrt v) successor s of w (entered via win)
      on the external face
      -- also return the direction in which s has been entered
 
      if MARK true mark the succ. vertex (and the edge)
      with mark (the visited field)
    */
{
    /*
      simply repeatedly calls embedg_VES_get_succ_on_ext_face
      until an active vertex is found
    */
    ASSERT(embedg_VES_is_vertex(n, w)
           || embedg_VES_is_virtual_vertex(n, w));

    embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                         w, win, MARK, mark, s, sin);
    while (embedg_VES_is_ver_inactive(embed_graph, n, v, *s))
    {
        embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                             *s, *sin, MARK, mark, s, sin);
    }
    ASSERT(!embedg_VES_is_ver_inactive(embed_graph, n, v, *s));
}

void 
embedg_VES_get_succ_ext_active_on_ext_face (t_ver_edge *embed_graph, int n,
	int v, int w, int win, boolean MARK, int mark, int *s, int *sin)
    /*
      find the externally active (wrt v) successor s of w (entered via win)
      on the external face
      -- also return the direction in which s has been entered
 
      if MARK true mark the succ. vertex (and the edge)
      with mark (the visited field)
    */
{
    ASSERT(embedg_VES_is_vertex(n, w)
           || embedg_VES_is_virtual_vertex(n, w));

    embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                         w, win, MARK, mark, s, sin);
    while (!embedg_VES_is_ver_ext_active(embed_graph, n, v, *s))
    {
        embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                             *s, *sin, MARK, mark, s, sin);
    }
    ASSERT(embedg_VES_is_ver_ext_active(embed_graph, n, v, *s));
}

void 
embedg_VES_get_succ_pertinent_on_ext_face (t_ver_edge *embed_graph, int n,
	int v, int w, int win, boolean MARK, int mark, int *s, int *sin)
    /*
      find the pertinent (wrt v) successor s of w (entered via win)
      on the external face
      -- also return the direction in which s has been entered
 
      if MARK true mark the succ. vertex (and the edge)
      with mark (the visited field)
    */
{
    ASSERT(embedg_VES_is_vertex(n, w)
           || embedg_VES_is_virtual_vertex(n, w));

    embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                         w, win, MARK, mark, s, sin);
    while (!embedg_VES_is_ver_pertinent(embed_graph, n, v, *s))
    {
        embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                             *s, *sin, MARK, mark, s, sin);
    }
    ASSERT(embedg_VES_is_ver_pertinent(embed_graph, n, v, *s));
}

/*
 *  mark_kur.c
 */
 
/*
  What:
  *****
  
  Implementing:

  Marking the Kuratowski obstruction (in the VES structure):
  this we do once we know which minor we are talking about
  (see isolator.c).
 
 
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 
 
  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca
 
 
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********
 
  Paulette Lieby (Magma), Brendan McKay (ANU)
 
  Started October 2001
*/
 
 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_VERB(x)   {}
#define IF_DEB_TREE(x)    {}
#define IF_DEB_EDGES(x) {}
#define IF_CPU(x) {}
 


/* aproto: header embed_graph_protos.h */

/* aproto: beginstatic -- don't touch this!! */
static void embedg_VES_walk_mark_part_ext_face (t_ver_edge *, int, int, int, int, int, int);
static void embedg_VES_walk_mark_ext_face (t_ver_edge *, int, int, int);
static void embedg_VES_walk_mark_part_proper_face (t_ver_edge *, int, int, int, int, int);
static boolean embedg_VES_is_part_ext_face_marked (t_ver_edge *, int, int, int, int, int, int);
static void embedg_get_u_x (t_ver_edge *, int, int, int, int *);
static int embedg_get_least_neigh (t_dlcl **, t_dlcl **, int, int, int);
static void embedg_add_mark_u_x (t_dlcl **, t_dlcl **, t_ver_edge *, int, int *, int, int, int *, int);
static void embedg_mark_tree_path (t_ver_edge *, int, int, int, int);
static void embedg_add_mark_v_w (t_dlcl **, t_dlcl **, t_ver_edge *, int, int *, int, int, int);
static void embedg_add_mark_v_w_for_B (t_dlcl **, t_dlcl **, t_ver_edge *, int, int *, int, int, int *, int);
static void embedg_mark_x_y_path (t_ver_edge *, int, int *, int *, int, int);
/* aproto: endstatic -- don't touch this!! */
 
#ifndef PLANAR_IN_MAGMA
#endif



static void 
embedg_VES_walk_mark_part_ext_face (t_ver_edge *embed_graph, int n,
	int v, int vin, int from, int to, int mark)
    /*
      walk & mark the external face:
      walk in the direction vin -> v -> vout and mark <from> <to>
    */
{
    int          cur, curin, next, nextin;

    embed_graph[from].visited = mark;
    embed_graph[to].visited = mark;

    IF_DEB(
           fprintf(stdout, "part. ext face marked\t");
           fprintf(stdout, "%d\t", from);
           )

    next = cur = v;
    curin = vin;
    while (next != from)
    {
        embedg_VES_get_succ_on_ext_face(embed_graph, n, cur, curin,
                                             FALSE, 0, &next, &nextin);
        cur = next;
        curin = nextin;
    }
    next = n;
    while (next != to)
    {
        embedg_VES_get_succ_on_ext_face(embed_graph, n, cur, curin,
                                             TRUE, mark, &next, &nextin);
        cur = next;
        curin = nextin;

        IF_DEB(
               fprintf(stdout, "%d\t", next);
               )
    }
    IF_DEB(
           fprintf(stdout, "\n");
           )
}
   
static void 
embedg_VES_walk_mark_ext_face (t_ver_edge *embed_graph, int n, int v, int mark)
    /*
      walk & mark the external face, starting & ending at vertex v
    */
{
    embedg_VES_walk_mark_part_ext_face(embed_graph, n, v, 0, v, v,
                                            mark);
}
   


static void 
embedg_VES_walk_mark_part_proper_face (t_ver_edge *embed_graph, int n,
	int from_e, int from_ein, int to, int mark)
    /*
      walk & mark a proper face starting at EDGE from_e and ending
      at VERTEX to

      walk in the direction from_ein -> from_e -> to and mark
      everything in between
    */
{
    int     s, cur_e, cur_ein, next_e, next_ein;
 
    next_e = s = n; /* this is an invalid value for an edge/vertex */
  
    cur_e = from_e;
    cur_ein = from_ein;
    while (s != to)
    {
        ASSERT(embedg_VES_is_edge(n, cur_e));
        ASSERT(!embedg_VES_is_short_cut_edge(embed_graph,
                                                        n, cur_e));
 
        embedg_VES_get_succ_on_proper_face(embed_graph, n,
                                                cur_e, cur_ein,
                                                TRUE, mark,
                                                &s, &next_e, &next_ein);
        cur_e = next_e;
        cur_ein = next_ein;
    }
}



static boolean 
embedg_VES_is_part_ext_face_marked (t_ver_edge *embed_graph, int n, int v,
	int vin, int from, int to, int mark)
    /*
      simple check to see if all the vertices on the external
      face walk starting at vin -> v -> vout <from> <to> are marked
      (with mark)
    */
{
    int          cur, curin, next, nextin;

    if (embed_graph[from].visited != mark || embed_graph[to].visited != mark)
        return FALSE;

    cur = v;
    curin = vin;
    next = n;
    while (next != from)
    {
        embedg_VES_get_succ_on_ext_face(embed_graph, n, cur, curin,
                                             FALSE, 0, &next, &nextin);
        cur = next;
        curin = nextin;
    }
    while (next != to)
    {
        embedg_VES_get_succ_on_ext_face(embed_graph, n, cur, curin,
                                             FALSE, 0, &next, &nextin);
        if (embed_graph[next].visited != mark)
            return FALSE;

        cur = next;
        curin = nextin;
    }
    
    return TRUE;
}


boolean 
embedg_VES_is_ext_face_marked (t_ver_edge *embed_graph, int n, int v, int mark)
    /*
      simple check to see if all the vertices on the external
      face walk starting/ending at v are marked (with mark)
    */
{
    return embedg_VES_is_part_ext_face_marked(embed_graph, n, v, 0,
                                                   v, v, mark);
}


static void 
embedg_get_u_x (t_ver_edge *embed_graph, int n, int v, int x, int *u_x)
    /*
      x is an externally active vertex (wrt v):
      we want u_x, the lowest point of "attachement" for
      the unembedded directed edge [x, u_x]
    */
{
    int          c;
    t_dlcl       *child_list;

    ASSERT(embedg_VES_is_ver_ext_active(embed_graph, n, v, x));
    if (embed_graph[x].least_ancestor < v)
        /*
          then there is a single unembedded back edge (u_x, x),
          u_x an ancestor of v
        */
    {
        *u_x = embed_graph[x].least_ancestor;
        return;
    }

    /*
      else there is a tree path x to d_x and an
      unembedded back edge (u_x, d_x)

      get the lowpoint of the first elt. in separated_DFS_child_list of x
    */
    child_list = embed_graph[x].separated_DFS_child_list;
    ASSERT(!embedg_dlcl_is_empty(child_list));
    c = child_list->info;
    *u_x = embed_graph[c].lowpoint;
}

static int 
embedg_get_least_neigh (t_dlcl **dfs_tree, t_dlcl **back_edges,
	int n, int v, int c)
    /*
      get the least neighbour of v >= c, ie a vertex in the sub tree
      rooted by c

      somehow this must always succeed
    */
{
    int          least_n;
    t_dlcl       *tree_l, *back_l, *p;

    /*
      neighbours are found in either dfs_tree[v] or back_edges[v]
    */

    tree_l = dfs_tree[v];
    back_l = back_edges[v];
    ASSERT(!embedg_dlcl_is_empty(tree_l) || !embedg_dlcl_is_empty(back_l));

    least_n = n;  /* ok, invalid value for any neighbour */
    p = tree_l;
    if (!embedg_dlcl_is_empty(p))
    {
        if (p->info >= c)
        {
            least_n = p->info < least_n ? p->info : least_n;
        }
        p = embedg_dlcl_list_next(p);
        while (p != tree_l)
        {
            if (p->info >= c)
            {
                least_n = p->info < least_n ? p->info : least_n;
            }
            p = embedg_dlcl_list_next(p);
        }
    }
    p = back_l;
    if (!embedg_dlcl_is_empty(p))
    {
        if (p->info >= c)
        {
            least_n = p->info < least_n ? p->info : least_n;
        }
        p = embedg_dlcl_list_next(p);
        while (p != back_l)
        {
            if (p->info >= c)
            {
                least_n = p->info < least_n ? p->info : least_n;
            }
            p = embedg_dlcl_list_next(p);
        }
    }

    ASSERT(least_n >= c);
    /*
      this is so because of the context where this function is called from
    */
    return least_n;
}
  
static void 
embedg_add_mark_u_x (t_dlcl **dfs_tree, t_dlcl **back_edges,
	t_ver_edge *embed_graph, int n, int *edge_pos, int v,
	int x, int *u_x, int mark)
    /*
      marking a Kuratowski homeomorph:
      
      marking and adding the unembedded dotted edge (u, x),
      x an ext. active vertex wrt v
    */
{
    int          c, d_x;
    t_dlcl       *child_list;

    ASSERT(embedg_VES_is_ver_ext_active(embed_graph, n, v, x));
    if (embed_graph[x].least_ancestor < v)
        /*
          then there is a single unembedded back edge (u_x, x),
          u_x an ancestor of v
        */
    {
        *u_x = embed_graph[x].least_ancestor;
        embed_graph[x].visited = mark;
        embed_graph[*u_x].visited = mark;
        embedg_VES_add_edge(embed_graph, n, edge_pos, *u_x, x,
                                       TRUE, mark);
        return;
    }

    /*
      else there is a tree path x to d_x and an
      unembedded back edge (u_x, d_x)

      get the lowpoint of the first elt. in separated_DFS_child_list of x
    */
    child_list = embed_graph[x].separated_DFS_child_list;
    ASSERT(!embedg_dlcl_is_empty(child_list));
    c = child_list->info;
    *u_x = embed_graph[c].lowpoint;

    /*
      search for the least neighbour of *u_x  >= c,
      that is in the subtree rooted by c
    */
    d_x = embedg_get_least_neigh(dfs_tree, back_edges, n, *u_x, c);
    ASSERT(d_x >= c);
    /*
      this must be true  since u_x is incident to a descendant of x
      (remember: x is externally active)
    */

    /*
      mark the DFS tree path from d_x to x
    */
    embedg_mark_tree_path(embed_graph, n, d_x, x, mark);
    /*
      add the unembedded (u_x, d_x) edge
    */
    embedg_VES_add_edge(embed_graph, n, edge_pos, *u_x, d_x,
                                   TRUE, mark);
}
  
static void 
embedg_mark_tree_path (t_ver_edge *embed_graph, int n, int d_x, int x, int mark)
    /*
      marking the DFS tree path d_x...x where x is an ancestor of d_x
    */
{
    int          cur_v, te, twe;

    ASSERT(d_x >= x);
    
    cur_v = d_x;

    while (cur_v != x)
    {
        embed_graph[cur_v].visited = mark;
        te = embed_graph[cur_v].link[0];
        ASSERT(embedg_VES_is_edge(n, te));
        while (!embedg_VES_is_tree_edge(embed_graph, n, te)
               || (embed_graph[te].neighbour > cur_v
                   && embed_graph[te].neighbour != cur_v + n))
            /*
              want to find a tree edge incident to an ancestor of d_x:
              given that d_x..x is a tree path, we MUST find such an edge!

              note also that I must take account of the fact that
              [te].neighbour could be a virtual vertex, in which case
              it can only be cur_v + n!
            */
        {
            te = embed_graph[te].link[0];
        }
        ASSERT(embedg_VES_is_tree_edge(embed_graph, n, te));
        ASSERT(embed_graph[te].neighbour == embed_graph[cur_v].DFS_parent
               || embed_graph[te].neighbour == cur_v + n);
        
        embed_graph[te].visited = mark;
        twe =  embedg_VES_get_twin_edge(embed_graph, n, te);
        embed_graph[twe].visited = mark;

        /*
          want only to deal with real vertices instead of virtual vertices
        */
        cur_v = embed_graph[te].neighbour < cur_v ?
            embed_graph[te].neighbour : embed_graph[cur_v].DFS_parent;
    }
    embed_graph[x].visited = MARK_MINORS(n);
}
  

static void 
embedg_add_mark_v_w (t_dlcl **dfs_tree, t_dlcl **back_edges,
	t_ver_edge *embed_graph, int n, int *edge_pos, int v, int w, int mark)
    /*
      marking a Kuratowski homeomorph:
      
      marking and adding the unembedded dotted edge (v, w),
      w is pertinent wrt v
    */
{
    int          vw, c, d_w;
    t_dlcl       *bicomp_list;
    
    if (embed_graph[w].adjacent_to == v)
        /*
          then there is a single unembedded back edge (v, w)
          w an ancestor of w
        */
    {
        embed_graph[v].visited = mark;
        embed_graph[w].visited = mark;
        embedg_VES_add_edge(embed_graph, n, edge_pos, v, w,
                                       TRUE, mark);
        return;
    }

    /*
      else there is a tree path w to d_w and an
      unembedded back edge (v, d_w)

      get the last elt in w's bicomp list
    */
    bicomp_list = embed_graph[w].pertinent_bicomp_list;
    ASSERT(!embedg_dlcl_is_empty(bicomp_list));
    vw = (embedg_dlcl_list_last(bicomp_list))->info;
    c = vw - n;

    /*
      search for the least neighbour of v >= c,
      that is in the subtree rooted by c
    */
    d_w = embedg_get_least_neigh(dfs_tree, back_edges, n, v, c);
    ASSERT(d_w >= c);
    /*
      this must be true since v is incident to a descendant of w
      (remember: w is pertinent)
    */

    /*
      mark the DFS tree path from d_w to w
    */
    embedg_mark_tree_path(embed_graph, n, d_w, w, mark);
    /*
      add the unembedded (d_w, v) edge
    */
    embedg_VES_add_edge(embed_graph, n, edge_pos, d_w, v,
                                   TRUE, mark);
}


static void 
embedg_add_mark_v_w_for_B (t_dlcl **dfs_tree, t_dlcl **back_edges,
	t_ver_edge *embed_graph, int n, int *edge_pos, int v, int w,
	int *u_z, int mark)
    /*
      marking a Kuratowski homeomorph:
      
      marking and adding the unembedded dotted edge (v, w) for minor B:
      w is pertinent wrt v
    */
{
    int          vz, z, d_z, d_w;
    t_dlcl       *bicomp_list;
    
    /*
      get the last elt in w's bicomp list
    */
    bicomp_list = embed_graph[w].pertinent_bicomp_list;
    ASSERT(!embedg_dlcl_is_empty(bicomp_list));
    vz = (embedg_dlcl_list_last(bicomp_list))->info;
    z = vz - n;

    /*
      get the lowpoint of z
    */
    *u_z = embed_graph[z].lowpoint;
 
    /*
      search for the least neighbour of *u_z  >= z,
      that is in the subtree rooted by c
    */
    d_z = embedg_get_least_neigh(dfs_tree, back_edges, n, *u_z, z);
    ASSERT(d_z >= z);
    /*
      this must be true since u_z is incident to z or a descendant of z
    */

    /*
      now do the same for neighbours of v
    */
    d_w = embedg_get_least_neigh(dfs_tree, back_edges, n, v, z);
    ASSERT(d_w >= z);
    /*
      this must be true since v is incident to a descendant of w
      (remember: w is pertinent)
    */
  
    /*
      mark the DFS tree path from d_w to w
    */
    embedg_mark_tree_path(embed_graph, n, d_w, w, mark);
    /*
      mark the DFS tree path from d_z to z
    */
    embedg_mark_tree_path(embed_graph, n, d_z, z, mark);
    /*
      add & mark the edges (u_z, d_z), (v, d_w)
    */
    embedg_VES_add_edge(embed_graph, n, edge_pos, *u_z, d_z,
                                   TRUE, mark);
    embedg_VES_add_edge(embed_graph, n, edge_pos, v, d_w,
                                   TRUE, mark);
}

static void 
embedg_mark_x_y_path (t_ver_edge *embed_graph, int n, int *path_v,
	int *path_e, int nbr_v, int mark)
{
    int          i;

    /*      
      have a look at embedg_iso_get_highest_x_y_path
      to see that path_e[0] is a dummy

      (note: path_v and path_e contain nbr_v + 1 elts!)
    */
    embed_graph[path_v[0]].visited = mark;
    for (i = 1; i <= nbr_v; i++)
    {
        int        e, twin;
        
        embed_graph[path_v[i]].visited = mark;
        e = path_e[i];
        twin = embedg_VES_get_twin_edge(embed_graph, n, e);
        embed_graph[e].visited =
            embed_graph[twin].visited = mark;
    }
}

void 
embedg_mark_minor_A (t_dlcl **dfs_tree, t_dlcl **back_edges,
	t_ver_edge *embed_graph, int n, int *edge_pos, int v, int c, int vr)
{
    int          r, r_c, x, y, w, u_x, u_y, u;

    ASSERT(embedg_VES_is_virtual_vertex(n, vr));
    r_c = vr - n;
    r = embed_graph[r_c].DFS_parent;

    /*
      find the ext. active x & y, and the pertinent w,
      and mark the external face of the bicomp rooted at vr
    */
    embedg_iso_get_x_y_w(embed_graph, n, v, r, r_c,
                              MARK_MINORS(n),
                              MARK_MINORS(n), MARK_MINORS(n), &x, &y, &w);

    /*
      mark the edges (u, x), (u, y), (v, w)
    */
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, x, &u_x,
                        MARK_MINORS(n));
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, y, &u_y,
                        MARK_MINORS(n));
    embedg_add_mark_v_w(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, w,
                        MARK_MINORS(n));

    /*
      mark the tree path from r to min(u_x, u_y)
    */
    u = u_x <= u_y ? u_x : u_y;
    embedg_mark_tree_path(embed_graph, n, r, u, MARK_MINORS(n));

    IF_DEB(
           fprintf(stdout, "mark minor A\n");
           fprintf(stdout, "v %d\t c %d\t r %d\t r_c %d\t x %d\t y %d\t w %d\t u_x %d\t u_y %d\n",
                   v, c, r, r_c, x, y, w, u_x, u_y);
           )
}

void 
embedg_mark_minor_B (t_dlcl **dfs_tree, t_dlcl **back_edges,
	t_ver_edge *embed_graph, int n, int *edge_pos, int v,
	int c, int x, int y, int w)
{
    int          vv, u_x, u_y, vz, u_z, u_max, u_min;

    vv = c + n;

    /*
      mark the external face of the bicomp rooted by v^c
    */
    embedg_VES_walk_mark_ext_face(embed_graph, n, vv, MARK_MINORS(n));
    ASSERT(embedg_VES_is_ext_face_marked(embed_graph, n, vv,
                                              MARK_MINORS(n)));

    /*
      mark the edges (u, x), (u, y)
    */
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, x, &u_x,
                        MARK_MINORS(n));
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, y, &u_y,
                        MARK_MINORS(n));

    /*
      mark the dotted edges (v, w), (v, u)
    */
    embedg_add_mark_v_w_for_B(dfs_tree, back_edges,
                              embed_graph, n, edge_pos, v, w,
                              &u_z, MARK_MINORS(n));

    /*
      mark the tree path from max(u_x, u_y, u_z) to min(u_x, u_y, u_z)
    */
    u_max = u_x > u_y ? u_x : u_y;
    u_max = u_max > u_z ? u_max : u_z;
    u_min = u_x < u_y ? u_x : u_y;
    u_min = u_min < u_z ? u_min : u_z;
    embedg_mark_tree_path(embed_graph, n, u_max, u_min, MARK_MINORS(n));

    IF_DEB(
           fprintf(stdout, "mark minor B\n");
           fprintf(stdout, "v %d\t c %d\t x %d\t y %d\t w %d\t u_x %d\t u_y %d\t u_z %d\n",
                   v, c, x, y, w, u_x, u_y, u_z);
           )
}

void 
embedg_mark_minor_C (t_dlcl **dfs_tree, t_dlcl **back_edges,
	t_ver_edge *embed_graph, int n, int *edge_pos, int v,
	int c, int x, int y, int w, int *path_v, int *path_e,
	int nbr_v, boolean px_attached_high, boolean py_attached_high)
{
    int          vv, p_x, p_y, u_x, u_y, u;

    vv = c + n;
    p_x = path_v[0];
    p_y = path_v[nbr_v];
    /*
      see embedg_iso_get_highest_x_y_path for the above
    */
    
    if (px_attached_high)
        /*
          mark the external face:
          - from v^c to p_y if py_attached_high
          - from v^c to y if !py_attached_high

          not too sure about that one....

          from v^c to p_y: so vvin = 0,
          in x's direction
        */
    {
        if (py_attached_high)
            embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 0,
                                                    vv, p_y, MARK_MINORS(n));
        else
            embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 0,
                                                    vv, y, MARK_MINORS(n));
    }
    else
        /*
          symmetric case:
          mark the external face from v^c to p_x: so vvin = 1,
          in y's direction
        */
    {
        if (px_attached_high)
            embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 1,
                                                    vv, p_x, MARK_MINORS(n));
        else
            embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 1,
                                                    vv, x, MARK_MINORS(n));
    }

    /*
      mark the edges (u, x), (u, y), (v, w)
    */
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, x, &u_x,
                        MARK_MINORS(n));
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, y, &u_y,
                        MARK_MINORS(n));
    embedg_add_mark_v_w(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, w,
                        MARK_MINORS(n));
 
    /*
      mark the tree path from v to min(u_x, u_y)
    */
    u = u_x <= u_y ? u_x : u_y;
    embedg_mark_tree_path(embed_graph, n, v, u, MARK_MINORS(n));

    /*
      finally, mark the x-y path, ie the vertices in path_v
      and the edges in path_e
    */
    embedg_mark_x_y_path(embed_graph, n, path_v, path_e, nbr_v,
                         MARK_MINORS(n));

    IF_DEB(
           fprintf(stdout, "mark minor C      p_x high %d\t p_y high %d\n",
                   px_attached_high, py_attached_high);
           fprintf(stdout, "v %d\t c %d\t x %d\t y %d\t w %d\t p_x %d\t p_y %d\t u_x %d\t u_y %d\n",
                   v, c, x, y, w, p_x, p_y, u_x, u_y);
           )
}

void 
embedg_mark_minor_D (t_dlcl **dfs_tree, t_dlcl **back_edges,
	t_ver_edge *embed_graph, int n, int *edge_pos, int v,
	int c, int x, int y, int w, int *path_v, int *path_e,
	int nbr_v, int entry_in_path_e)
{
    int          i, vv, p_x, p_y, u_x, u_y, u;

    vv = c + n;
    p_x = path_v[0];
    p_y = path_v[nbr_v];
    /*
      see embedg_iso_get_highest_x_y_path for the above
    */
    
    /*
      mark the lower external face from x to y: we can walk in
      either direction
    */
    embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 0,
                                            x, y, MARK_MINORS(n));

    /*
      mark the internal path which goes from the x-y path to v
      - since I haven't stored those vertices/edges I assume
      that a proper face walk should suffice

      BUT a walk that say starts at p_x and ends at vv,
      that is, a walk starting at path_e[1] entered from entry_in_path_e
      (recall that path_e[0] is a dummy)
    */
    embedg_VES_walk_mark_part_proper_face(embed_graph, n,
                                               path_e[1], entry_in_path_e,
                                               vv, MARK_MINORS(n));

    /*
      a note of caution here:
      ALWAYS mark external/internal faces before adding any other edges:
      since adding edges destroys the faces' consistency
      (adding edges makes no sense of face since we are in a non-planar
      situation)
    */
    /*
      mark the edges (u, x), (u, y), (v, w)
    */
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, x, &u_x,
                        MARK_MINORS(n));
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, y, &u_y,
                        MARK_MINORS(n));
    embedg_add_mark_v_w(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, w,
                        MARK_MINORS(n));
 
    /*
      mark the tree path from v to min(u_x, u_y)
    */
    u = u_x <= u_y ? u_x : u_y;
    embedg_mark_tree_path(embed_graph, n, v, u, MARK_MINORS(n));

    /*
      mark the x-y path, ie the vertices in path_v
      and the edges in path_e
    */
    embedg_mark_x_y_path(embed_graph, n, path_v, path_e, nbr_v,
                         MARK_MINORS(n));

    IF_DEB(
           fprintf(stdout, "mark minor D\n");
           fprintf(stdout, "v %d\t c %d\t x %d\t y %d\t w %d\t p_x %d\t p_y %d\t u_x %d\t u_y %d\n",
                   v, c, x, y, w, p_x, p_y, u_x, u_y);
           )
}




minor 
embedg_mark_minor_E (t_dlcl **dfs_tree, t_dlcl **back_edges,
	t_ver_edge *embed_graph, int n, int *edge_pos, int v,
	int c, int x, int y, int w, int *path_v, int *path_e, int nbr_v)
    /*
      while marking minor E return which of the minors we are dealing with
    */
{
    int          vv, p_x, p_y, u_x, u_y, u_w, u, u_max, u_min;

    vv = c + n;
    p_x = path_v[0];
    p_y = path_v[nbr_v];
    /*
      see embedg_iso_get_highest_x_y_path for the above
    */

    if (!embedg_VES_is_ver_ext_active(embed_graph, n, v, w))
        /*
          minor E1 case: we must find an ext. active z, distinct from w,
          on the external face p_x..w..p_y
        */
    {
        int       s, sin, cur, curin, z, u_z, u_xy;
        
        s = n;
        /*
          start searching at vv entered from 0 (in x's direction)
          -- we MUST reach p_x - hopefully! :)
        */
        cur = vv;
        curin = 0;
        while (s != p_x)
            /*
              first advance to p_x: we are sure of reaching it
            */
        {
            embedg_VES_get_succ_on_ext_face(embed_graph, n,
                                                 cur, curin,
                                                 FALSE, 0, &s, &sin);
            cur = s;
            curin = sin;
        }

        /*
          continue the walk on the external face:
          stop if either s is ext. active OR s == w

          we'll mark the lot later on
        */
        while (
               !(embedg_VES_is_ver_ext_active(embed_graph, n, v,
                                                         s)
                 && s != p_x)
               && s != w)
        {
            embedg_VES_get_succ_on_ext_face(embed_graph, n, cur, curin,
                                                 FALSE, 0, &s, &sin);
            cur = s;
            curin = sin;
        }
        /*
          now we must decide which symmetry we are in
        */
        if (embedg_VES_is_ver_ext_active(embed_graph, n, v, s))
            /*
              z is between x and w (recall that w is NOT ext. active)
            */
        {
            z = s;
            ASSERT(z != w);

            /*
              mark the external face from v^c to y in x's direction
            */
            embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 0,
                                                    vv, y, MARK_MINORS(n));
            /*
              add/mark dotted edge (u, y)
            */
            embedg_add_mark_u_x(dfs_tree, back_edges,
                                embed_graph, n, edge_pos,
                                v, y, &u_xy, MARK_MINORS(n));
        }
        else
            /*
              this is the symmetric case: must find z between w and p_y
            */
        {
            ASSERT(s == w);
            embedg_VES_get_succ_ext_active_on_ext_face(embed_graph, n, 
                                                            v, cur, curin,
                                                            FALSE, 0,
                                                            &s, &sin);
            /*
              and z is distinct from p_y!
            */
            z = s;
            ASSERT(z != p_y);

            /*
              mark the external face from v^c to x in y's direction
            */
            embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 1,
                                                    vv, x, MARK_MINORS(n));
            /*
              add/mark dotted edge (u, x)
            */
            embedg_add_mark_u_x(dfs_tree, back_edges,
                                embed_graph, n, edge_pos,
                                v, x, &u_xy, MARK_MINORS(n));
        }
        /*
          now the marked bits which are common to both cases:
          dotted edges (u, z), (v, w), the x-y path,
          the tree path (v, min(u_xy, u_z))
        */
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, z, &u_z, MARK_MINORS(n));
        embedg_add_mark_v_w(dfs_tree, back_edges,
                            embed_graph, n, edge_pos, v, w,
                            MARK_MINORS(n));

        embedg_mark_x_y_path(embed_graph, n, path_v, path_e, nbr_v,
                             MARK_MINORS(n));
        
        u = u_z <= u_xy ? u_z : u_xy;
        embedg_mark_tree_path(embed_graph, n, v, u, MARK_MINORS(n));

        IF_DEB(
               fprintf(stdout, "mark minor E1\n");
               fprintf(stdout, "v %d\t c %d\t x %d\t y %d\t z %d\t w %d\t p_x %d\t p_y %d\t u_xy %d\t u_z %d\n",
                       v, c, x, y, z, w, p_x, p_y, u_xy, u_z);
           )

        return MINOR_E1;
    }
    
    /*
      in all other cases we get u_x, u_y, u_w back
      from the ext. active vertices x, y, w resp.

      again, I CANNOT embed these edges now since that would destroy
      my external/internal faces
    */
    
    embedg_get_u_x(embed_graph, n, v, x, &u_x);
    embedg_get_u_x(embed_graph, n, v, y, &u_y);
    embedg_get_u_x(embed_graph, n, v, w, &u_w);

    if (u_w > u_x && u_w > u_y)
        /*
          minor E2 case:
          we mark the whole external face rooted by v^c
          and the tree path (v, min(u_x, u_y))
        */
    {
        embedg_VES_walk_mark_ext_face(embed_graph, n, vv,
                                           MARK_MINORS(n));
        u = u_x <= u_y ? u_x : u_y;
        embedg_mark_tree_path(embed_graph, n, v, u, MARK_MINORS(n));

        /*
          embed dotted edges (u, x), (u, y) & (u, w)
        */
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, x, &u_x, MARK_MINORS(n));
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, y, &u_y, MARK_MINORS(n));
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, w, &u_w, MARK_MINORS(n));
        
        IF_DEB(
               fprintf(stdout, "mark minor E2\n");
               fprintf(stdout, "v %d\t c %d\t x %d\t y %d\t w %d\t p_x %d\t p_y %d\t u_x %d\t u_y %d\t u_w %d\n",
                       v, c, x, y, w, p_x, p_y, u_x, u_y, u_w);
           )
            
        return MINOR_E2;
    }

    /*
      two more things common to all remaining cases:

      mark the dotted edge (v, w) (but we MUST do that later)
      
      mark the x-y path
    */
    embedg_mark_x_y_path(embed_graph, n, path_v, path_e, nbr_v,
                         MARK_MINORS(n));
     
    if (u_x < u_y && u_w < u_y)
        /*
          minor E3 case: one of the symmetric cases:
          the external face rooted at v_c from vv to x (in x's direction)
          the external face rooted at v_c from y to w (in y's direction)
          the (v, min(u_w, u_x)) tree path
        */
    {
        embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 0,
                                                vv, p_x, MARK_MINORS(n));
        embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 1,
                                                y, w, MARK_MINORS(n));
        
        u = u_x <= u_w ? u_x : u_w;
        embedg_mark_tree_path(embed_graph, n, v, u, MARK_MINORS(n));

        /*
          embed dotted edges (u, x), (u, y), (u, w), (v, w)
        */
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, x, &u_x, MARK_MINORS(n));
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, y, &u_y, MARK_MINORS(n));
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, w, &u_w, MARK_MINORS(n));
        embedg_add_mark_v_w(dfs_tree, back_edges,
                            embed_graph, n, edge_pos, v, w,
                            MARK_MINORS(n));

        IF_DEB(
               fprintf(stdout, "mark minor E3/a\n");
               fprintf(stdout, "v %d\t c %d\t x %d\t y %d\t w %d\t p_x %d\t p_y %d\t u_x %d\t u_y %d\t u_w %d\n",
                       v, c, x, y, w, p_x, p_y, u_x, u_y, u_w);
           )
            
        return MINOR_E3;
    }
    if (u_y < u_x && u_w < u_x)
        /*
          minor E3 case: the other symmetric case:
          the external face rooted at v_c from vv to y (in y's direction)
          the external face rooted at v_c from x to w (in x's direction)
          the (v, min(u_w, u_y)) tree path
        */
    {
        embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 1,
                                                vv, p_y, MARK_MINORS(n));
        embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 0,
                                                x, w, MARK_MINORS(n));

        u = u_y <= u_w ? u_y : u_w;
        embedg_mark_tree_path(embed_graph, n, v, u, MARK_MINORS(n));

        /*
          embed dotted edges (u, x), (u, y), (u, w), (v, w)
        */
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, x, &u_x, MARK_MINORS(n));
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, y, &u_y, MARK_MINORS(n));
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, w, &u_w, MARK_MINORS(n));
        embedg_add_mark_v_w(dfs_tree, back_edges,
                            embed_graph, n, edge_pos, v, w,
                            MARK_MINORS(n));

        IF_DEB(
               fprintf(stdout, "mark minor E3/b\n");
               fprintf(stdout, "v %d\t c %d\t x %d\t y %d\t w %d\t p_x %d\t p_y %d\t u_x %d\t u_y %d\t u_w %d\n",
                       v, c, x, y, w, p_x, p_y, u_x, u_y, u_w);
           )
            
        return MINOR_E3;
    }
    
    if (p_x != x)
        /*
          minor E4 case: one of the symmetric cases:
          the external face rooted at v_c from vv to w (in x's direction)
          the external face rooted at v_c from vv to p_y (in y's direction)
          the tree path from max(u_x, u_y, u_w) to min(u_x, u_y, u_w)
        */
    {
        embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 0,
                                                vv, w, MARK_MINORS(n));
        embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 1,
                                                vv, p_y, MARK_MINORS(n));
        
        u_max = u_x > u_y ? u_x : u_y;
        u_max = u_max > u_w ? u_max : u_w;
        u_min = u_x < u_y ? u_x : u_y;
        u_min = u_min < u_w ? u_min : u_w;
        embedg_mark_tree_path(embed_graph, n, u_max, u_min, MARK_MINORS(n));

        /*
          embed dotted edges (u, x), (u, y), (u, w), (v, w)
        */
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, x, &u_x, MARK_MINORS(n));
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, y, &u_y, MARK_MINORS(n));
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, w, &u_w, MARK_MINORS(n));
        embedg_add_mark_v_w(dfs_tree, back_edges,
                            embed_graph, n, edge_pos, v, w,
                            MARK_MINORS(n));

        IF_DEB(
               fprintf(stdout, "mark minor E4/a\n");
               fprintf(stdout, "v %d\t c %d\t x %d\t y %d\t w %d\t p_x %d\t p_y %d\t u_x %d\t u_y %d\t u_w %d\n",
                       v, c, x, y, w, p_x, p_y, u_x, u_y, u_w);
           )
            
        return MINOR_E4;
    }
    if (p_y != y)
        /*
          minor E4 case: the other symmetric case:
          the external face rooted at v_c from vv to w (in y's direction)
          the external face rooted at v_c from vv to x (in x's direction)
          (here p_x = x!)
          the tree path from max(u_x, u_y, u_w) to min(u_x, u_y, u_w)
        */
    {
        embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 1,
                                                vv, w, MARK_MINORS(n));
        embedg_VES_walk_mark_part_ext_face(embed_graph, n, vv, 0,
                                                vv, x, MARK_MINORS(n));
        
        u_max = u_x > u_y ? u_x : u_y;
        u_max = u_max > u_w ? u_max : u_w;
        u_min = u_x < u_y ? u_x : u_y;
        u_min = u_min < u_w ? u_min : u_w;
        embedg_mark_tree_path(embed_graph, n, u_max, u_min, MARK_MINORS(n));

        /*
          embed dotted edges (u, x), (u, y), (u, w), (v, w)
        */
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, x, &u_x, MARK_MINORS(n));
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, y, &u_y, MARK_MINORS(n));
        embedg_add_mark_u_x(dfs_tree, back_edges,
                            embed_graph, n, edge_pos,
                            v, w, &u_w, MARK_MINORS(n));
        embedg_add_mark_v_w(dfs_tree, back_edges,
                            embed_graph, n, edge_pos, v, w,
                            MARK_MINORS(n));

        IF_DEB(
               fprintf(stdout, "mark minor E$/b\n");
               fprintf(stdout, "v %d\t c %d\t x %d\t y %d\t w %d\t p_x %d\t p_y %d\t u_x %d\t u_y %d\t u_w %d\n",
                       v, c, x, y, w, p_x, p_y, u_x, u_y, u_w);
           )
            
        return MINOR_E4;
    }

    /*
      this is the last case for minor E: when the homeomorph is K5
      
      mark the whole external face rooted at v^c
      mark the tree path from v to min(u_x, u_y, u_w)
    */

    embedg_VES_walk_mark_ext_face(embed_graph, n, vv, MARK_MINORS(n));
    
    u = u_x < u_y ? u_x : u_y;
    u = u < u_w ? u : u_w;
    embedg_mark_tree_path(embed_graph, n, v, u, MARK_MINORS(n));

    /*
      embed dotted edges (u, x), (u, y), (u, w), (v, w)
    */
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos,
                        v, x, &u_x, MARK_MINORS(n));
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos,
                        v, y, &u_y, MARK_MINORS(n));
    embedg_add_mark_u_x(dfs_tree, back_edges,
                        embed_graph, n, edge_pos,
                        v, w, &u_w, MARK_MINORS(n));
    embedg_add_mark_v_w(dfs_tree, back_edges,
                        embed_graph, n, edge_pos, v, w,
                        MARK_MINORS(n));
    
    IF_DEB(
           fprintf(stdout, "mark minor E5\n");
           fprintf(stdout, "v %d\t c %d\t x %d\t y %d\t w %d\t p_x %d\t p_y %d\t u_x %d\t u_y %d\t u_w %d\n",
                   v, c, x, y, w, p_x, p_y, u_x, u_y, u_w);
           )
            
    return MINOR_E5;
}
/*
 *  proper_face_walk.c
 */
 
/*
  What:
  *****
  
  Implementing a proper face walk within the VES structure.
  This is obviously not the same as an external face walk,
  but is simply the standard face walk in a planar embedding.

  Not much to say, if only to emphasize that for our
  purposes here we assume:

  1. the short-cut edges have been removed from the VES structure
  2. each vertex/edge has been given its orientation
  2. the adjacency lists (vertex + plus its incident edges)
     are consistent:  (and this is IMPORTANT)
     that is, the way to traverse an adj. list (ie what
     constitute previous and next in the list which actually
     is a planar embedding at this stage) is indicated
     by the vertex/edge's orientation


     try to explain this better another time.... sorry...
  


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  from 

  Simplified O(n) Planarity Algorithms  (draft)
  ************************************
  
  John Boyer      JBoyer@PureEdge.com, jboyer@acm.org
  Wendy Myrvold   wendym@csr.uvic.ca


  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  authors:
  ********

  Paulette Lieby (Magma), Brendan McKay (ANU)

  Started October 2001
*/

 
#include "planarity.h"
 
#define IF_DEB(x)    {}
#define IF_DEB_PROPER_FACE(x)    {}
#define IF_VERB(x)   {}

 
 
/* aproto: header embed_graph_protos.h */




boolean 
embedg_VES_get_succ_on_proper_face_with_avoidance (t_ver_edge *embed_graph,
	int n, int e, int ein, int a, boolean MARK, int mark, int *s,
	int *next_e, int *next_ein)
    /*
      find the successor s of embed_graph[e].neighbour
      (entered via ein) on a proper face traversal
      which avoids (the vertex) a if a != n 

      also returns the edge next_e such that
      embed_graph[next_e].neighbour = s (to allow for continuation
      of the walk)

      assumes that short-cut edges have been removed and that each
      edge/vertex has been given its orientation

      and (more importantly) assumes that adjacency lists are consistent

      this function has been written especially to retrieve the highest
      x-y path for the isolator;
      (see embedg_iso_get_highest_x_y_path)
      but as I discovered later (when marking an internal face
      as for minor D) this function is general purpose

      PLUS: return true if the proper face walk has to skip an edge
      incident to a (ie had to "avoid" a)

      PLUS: mark s & next_e if so requested
    */
{
    int            eout;
    int            twin, twinout;
    boolean        avoid_a;

    ASSERT(embedg_VES_is_edge(n, e));
    ASSERT(!embedg_VES_is_short_cut_edge(embed_graph, n, e));

    IF_DEB(
           fprintf(stdout, "get_succ_on_proper_face, \n");
           )

    avoid_a = FALSE;
    /*
      find the direction out of the edge
    */
    eout = 1 ^ ein;

    /*
      get the twin edge
    */
    twin = embedg_VES_get_twin_edge(embed_graph, n, e);
 
    /*
      for each edge we must set the way to get to the next
      in the adjacency list:
      adjacency lists are traversed according to the vertex/edges
      orientation (one unique orientation per list of course)
    */
    if (embed_graph[e].sign != embed_graph[twin].sign)
        /*
          invert traversal
        */
    {
        twinout = 1 ^ eout;
    }
    else
        /*
          traversal is identical
        */
    {
        twinout = eout;
    }
 
    /*
      now, we want the edge previous to twin in twin's adjacency list,
      ie link[1 ^ twinout]
    */
    *next_e = embed_graph[twin].link[1 ^ twinout];
    /*
      next_e could be a vertex, I need an edge
    */
    if (embedg_VES_is_vertex(n, *next_e)
        || embedg_VES_is_virtual_vertex(n, *next_e))
        /*
          at this stage all virtual vertices should have
          been disabled BUT the vertices rooting the bicomps!!!
        */
    {
        *next_e = embed_graph[*next_e].link[1 ^ twinout];
    }
    ASSERT(embedg_VES_is_edge(n, *next_e));
    ASSERT(!embedg_VES_is_short_cut_edge(embed_graph, n, e));
    *s = embed_graph[*next_e].neighbour;

    if (*s == a)
        /*
          want to avoid this vertex, so must get yet previous
          edge in adjacency list
        */
    {
        avoid_a = TRUE;
        
        *next_e = embed_graph[*next_e].link[1 ^ twinout];
        if (embedg_VES_is_vertex(n, *next_e)
            || embedg_VES_is_virtual_vertex(n, *next_e))
        {
            *next_e = embed_graph[*next_e].link[1 ^ twinout];
        }
        ASSERT(embedg_VES_is_edge(n, *next_e));
        ASSERT(!embedg_VES_is_short_cut_edge(embed_graph, n, e));
    }
    *s = embed_graph[*next_e].neighbour;
    ASSERT(*s != a);

    /*
      finally (again, because lists  are consistent)
    */
    *next_ein = 1 ^ twinout;

    /*
      now mark s and next_e if required
    */
    if (MARK)
    {
        embed_graph[*s].visited = 
            embed_graph[*next_e].visited = mark;
        /*
          ouuh... must mark the twin as well....
          but ONLY when we mark the minors....
          that is poor design, can we do better????
          -- don't think so...

          (when we mark when counting the faces, we MUST only
          mark the edge and NOT its twin)
        */
        if (mark == MARK_MINORS(n))
        {
            twin =
                embedg_VES_get_twin_edge(embed_graph, n, *next_e);
            embed_graph[twin].visited = mark;
        }
    }

    return avoid_a;
}



void 
embedg_VES_get_succ_on_proper_face (t_ver_edge *embed_graph, int n, int e,
	int ein, int MARK, int mark, int *s, int *next_e, int *next_ein)
    /*
      same as above but without avoidance
    */
{
    boolean  avoid;
    
    avoid =
        embedg_VES_get_succ_on_proper_face_with_avoidance(
                                                       embed_graph, n,
                                                       e, ein, n,
                                                       MARK, mark,
                                                       s, next_e, next_ein);
    ASSERT(avoid == FALSE);
}


void 
embedg_VES_walk_proper_face (t_ver_edge *embed_graph, int n, int e,
	int ein, boolean MARK, int mark)
    /*
      traversing a proper face starting at edge e which has been entered
      via ein

      -- we mark the visited edges with mark if so requested

      assumes that short-cut edges have been removed and that each
      edge/vertex has been given its orientation
    */
{
    int     s, cur_e, cur_ein, next_e, next_ein;

    next_e = n; /* this is an invalid value for an edge */

    IF_DEB_PROPER_FACE(
                       fprintf(stdout, "proper face traversal\n");
                       )

    cur_e = e;
    cur_ein = ein;
    while (next_e != e)
    {
        ASSERT(embedg_VES_is_edge(n, cur_e));
        ASSERT(!embedg_VES_is_short_cut_edge(embed_graph,
                                                        n, cur_e));
        IF_DEB_PROPER_FACE(
                           embedg_VES_print_edge(embed_graph, n, cur_e);
                           )

        embedg_VES_get_succ_on_proper_face(embed_graph, n,
                                                cur_e, cur_ein,
                                                MARK, mark,
                                                &s, &next_e, &next_ein);
        cur_e = next_e;
        cur_ein = next_ein;
    }

    /*
      note that by doing so we would have marked e and the first of e's
      endpoints since by exiting the loop e = next_e and s is the
      actual starting vertex of the walk
    */
}





