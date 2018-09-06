/*
  data structures and stuff for the planarity algorithm

  Paulette Lieby, Brendan McKay
  October 2001
*/

#ifndef _PLANARITY_H_
#define _PLANARITY_H_



/* The following line must be uncommented for compiling into Magma. */
/* #define PLANAR_IN_MAGMA  */


#ifdef PLANAR_IN_MAGMA
#include "defs.h"
#include "system.h"   /* includes <stdio.h> <signal.h> "system_math.h"
                         <setjmp.h> <ctype.h> and more
                      */
#else
/* not PLANAR_IN_MAGMA */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#undef FALSE
#undef TRUE
#define FALSE 0
#define TRUE  1

#define NP NULL

#define ASSERT(x) assert(x)
#define DIE() exit(0)

#define mem_malloc malloc
#define mem_realloc realloc
#define mem_free free

#include "naututil.h"
#ifdef CPUDEFS
CPUDEFS
#endif
#define time_current_user() CPUTIME

#endif  /* not  PLANAR_IN_MAGMA */

#include "nauty.h" 





/*
  max number of edges (and directed edges) for the embed_graph
  data structure:
  1 more than for a (possibly) planar graph to allow search for obstructions

  1. if the graph is planar the embedding cannot possibly contain
  more than 3*n - 6 edges (including the short cut edges)
  2. if the graph is non planar, when retrieving and marking the
  obstruction, we introduce EXACTLY one edge crossing
*/
#define MAXE(n) ((n) > 1 ? 3*(n) - 5 : 0)
#define MAXDE(n) (6*(n) - 10)

#define NIL -1
#define CUTV -2     /* obviously used in diff. circ. than NILSIGN */
#define NILSIGN -2
#define CCLOCKW 1
#define CLOCKW -1
#define TE 1
#define BE 2
#define SCE 3

/*
  various "marks" for various purposes, ONLY for the t_ver_edge str

  note: do NOT use a mark in {0,..,n} since those values are
  used either while initialising or in the walkup procedure
*/
#define MARK_EMBED(n) ((n)+1)
#define MARK_EXT_FACE(n) ((n)+2)
#define MARK_EXT_FACE_L(n) ((n)+3)
#define MARK_EXT_FACE_R(n) ((n)+4)
#define MARK_X_Y_PATH(n) ((n)+5)
#define MARK_MINORS(n) ((n)+6)
#define MIN_EMBED_MARK    0   /* ONLY for the t_embed_sparse_rep str */
                        

typedef enum
{
    MINOR_A,
    MINOR_B,
    MINOR_C,
    MINOR_D,
    MINOR_E,
    MINOR_E1,
    MINOR_E2,
    MINOR_E3,
    MINOR_E4,
    MINOR_E5,
    NBR_MINORS
} minor;


/*
  a basic doubly linked circular list storing
  vertices/edges from the "big" 2*n + 2(3*n-5) array of vertices/edges

  only used internally in the planarity tester: especially
  where ordering of the vertices is important
*/
typedef struct dlcl {
    int            info;
    /*
      info is:
      - position in "big" array
    */
    int            in_adjl;      /* if relevant, the pos in the adjl. list
                                    of this edge
                                 */
    int            twin_in_adjl; /* if relevant, the pos in the adjl. list
                                    of the twin of this edge
                                 */
    int            mult;         /* if relevant, #occurences if this edge
                                    (when graph is not simple
                                 */
    struct dlcl    * right;
    struct dlcl    * left;
} t_dlcl;


/*
  a common structure for both (virtual) vertex & edge
*/
typedef struct ver_edge {
    /* vertex data */
    int          label;
    int          DFS_parent;
    int          least_ancestor;
    int          lowpoint;
    t_dlcl       * separated_DFS_child_list;
    t_dlcl       * rep_in_parent_list;
    t_dlcl       * pertinent_bicomp_list;
    int          adjacent_to;
    int          visited;
    /* edge data */
    int          neighbour;
    int          in_adjl;
    int          twin_in_adjl; 
    int          mult;
    int          type;
    int          sign;
    /* link the lot in a doubly linked circular list           */
    /* links indicate positions in the array of vertices/edges */
    int          link[2];
} t_ver_edge;


/*
  data structure for the merge queue
*/
typedef struct merge_queue {
    int          start, end;
    int          *b;
} t_merge_queue;



/*
  data structure for the sparse graph representation:
  the array of vertices
*/
typedef struct ver_sparse_rep {
    int          first_edge; /* can be index into an adj. list
                                or an embedding */
} t_ver_sparse_rep;

/*
  data structure for the sparse graph representation:
  a record in the adjacency list
*/
typedef struct adjl_sparse_rep {
    int          end_vertex; 
    int          next;   /* next in list as an index in the adj. list */
} t_adjl_sparse_rep;

/*
  data structure for the sparse graph representation:
  a record in the embedding
*/
typedef struct embed_sparse_rep {
    int          in_adjl;/* index of edge in adj. list */
    int          next;   /* next edge in embedding  */
    int          prev;   /* previous edge in embedding  */
    int          inv;    /* inverse edge */
    int          mark;   /* a spot for marking */
} t_embed_sparse_rep;

/*
  data structure for the sparse graph representation:
  a record an individual edge
*/
typedef struct edge_sparse_rep {
    int          ends[2];
} t_edge_sparse_rep;

/*
  data structure for the sparse graph representation:
  a record for a component
*/
typedef struct comp_sparse_rep {
    int          nbr_v;   /* nbr of vertices */
    int          *v;      /* the actual vertices */
} t_comp_sparse_rep;

typedef struct graph_sparse_rep {
        t_ver_sparse_rep    *V;
	int                 n;
	t_adjl_sparse_rep   *A;
	int                 size_A;
	int                 pos_A;
	int                 nbr_e;   /* ALWAYS # directed edges */
} t_graph_sparse_rep;

#ifdef __cplusplus
extern "C" {
#endif
/*
 * embed_graph_protos.h
 */

/* aproto: file embed_graph/sparseg_adjl_pred.c */
extern boolean sparseg_adjl_dir_edge_exists (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, int);
extern boolean sparseg_adjl_u_adj_v (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, int);
extern boolean sparseg_adjl_sub (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, t_ver_sparse_rep *, int, t_adjl_sparse_rep *);
extern boolean sparseg_adjl_eq (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, t_ver_sparse_rep *, int, t_adjl_sparse_rep *);
/* aproto: endfile */
/* aproto: file embed_graph/VES_misc.c */
extern boolean embedg_VES_is_vertex (int, int);
extern boolean embedg_VES_is_virtual_vertex (int, int);
extern boolean embedg_VES_is_edge (int, int);
extern boolean embedg_VES_is_tree_edge (t_ver_edge *, int, int);
extern boolean embedg_VES_is_back_edge (t_ver_edge *, int, int);
extern boolean embedg_VES_is_short_cut_edge (t_ver_edge *, int, int);
extern void embedg_VES_print_vertex (int, int);
extern void embedg_VES_print_virtual_vertex (t_ver_edge *, int, int);
extern void embedg_VES_print_any_vertex (t_ver_edge *, int, int);
extern void embedg_VES_print_any_rec (t_ver_edge *, int, int);
extern void embedg_VES_print_edge (t_ver_edge *, int, int);
extern void embedg_VES_print_flipped_edges (t_ver_edge *, int, int);
extern int embedg_VES_get_edge_from_ver (t_ver_edge *, int, int);
extern int embedg_VES_get_ver_from_edge (t_ver_edge *, int, int);
extern int embedg_VES_get_twin_edge (t_ver_edge *, int, int);
extern int embedg_VES_get_ver_from_virtual (t_ver_edge *, int, int);
extern int embedg_VES_get_ver (t_ver_edge *, int, int);
extern int embedg_VES_get_next_in_dlcl (t_ver_edge *, int, int, int);
extern void embedg_VES_walk_bicomp (t_ver_edge *, int, int, int);
extern void embedg_VES_print_adj_list (t_ver_edge *, int, int, boolean);
extern boolean embedg_VES_is_adj_list_consistent (t_ver_edge *, int, int);
extern boolean embedg_VES_are_adj_lists_consistent (t_ver_edge *, int);
extern void embedg_VES_remove_edge (t_ver_edge *, int, int);
extern void embedg_VES_set_orientation (t_ver_edge *, int, int *);
/* aproto: endfile */
/* aproto: file embed_graph/embed_edge.c */
extern void embedg_VES_embed_edge (t_ver_edge *, int, int *, int, int, int, int, int);
extern void embedg_VES_add_edge (t_ver_edge *, int, int *, int, int, boolean, int);
/* aproto: endfile */
/* aproto: file embed_graph/embedg_misc.c */
extern void embedg_VES_delete (t_ver_edge *, int);
extern void embedg_VES_print (t_ver_edge *, int);
extern void embedg_VES_print_bigcomps (t_ver_edge *, int);
/* aproto: endfile */
/* aproto: file embed_graph/ext_face_walk.c */
extern void embedg_VES_get_succ_on_ext_face (t_ver_edge *, int, int, int, boolean, int, int *, int *);
extern void embedg_VES_get_succ_active_on_ext_face (t_ver_edge *, int, int, int, int, boolean, int, int *, int *);
extern void embedg_VES_get_succ_ext_active_on_ext_face (t_ver_edge *, int, int, int, int, boolean, int, int *, int *);
extern void embedg_VES_get_succ_pertinent_on_ext_face (t_ver_edge *, int, int, int, int, boolean, int, int *, int *);
/* aproto: endfile */
/* aproto: file embed_graph/isolator.c */
extern int embedg_iso_get_c_of_v (t_ver_edge *, int, int, int);
extern boolean embedg_iso_is_minor_A (t_ver_edge *, int, int *, int, int, int *);
extern void embedg_iso_get_x_y_w (t_ver_edge *, int, int, int, int, int, int, int, int *, int *, int *);
extern boolean embedg_iso_is_minor_B (t_ver_edge *, int, int *, int, int, int *, int *, int *);
extern void embedg_iso_get_highest_x_y_path (t_ver_edge *, int, int, int, int, int, int, int, int, int, int **, int **, int *, int *, boolean *, boolean *, boolean *);
/* aproto: endfile */
/* aproto: file embed_graph/mark_kur.c */
extern boolean embedg_VES_is_ext_face_marked (t_ver_edge *, int, int, int);
extern void embedg_mark_minor_A (t_dlcl **, t_dlcl **, t_ver_edge *, int, int *, int, int, int);
extern void embedg_mark_minor_B (t_dlcl **, t_dlcl **, t_ver_edge *, int, int *, int, int, int, int, int);
extern void embedg_mark_minor_C (t_dlcl **, t_dlcl **, t_ver_edge *, int, int *, int, int, int, int, int, int *, int *, int, boolean, boolean);
extern void embedg_mark_minor_D (t_dlcl **, t_dlcl **, t_ver_edge *, int, int *, int, int, int, int, int, int *, int *, int, int);
extern minor embedg_mark_minor_E (t_dlcl **, t_dlcl **, t_ver_edge *, int, int *, int, int, int, int, int, int *, int *, int);
/* aproto: endfile */
/* aproto: file embed_graph/merge_bicomps.c */
extern void embedg_VES_merge_simple_bicomps (t_ver_edge *, int, int, int, int, int);
extern void embedg_VES_merge_pertinent_bicomps (t_ver_edge *, int, int, int, int, int);
/* aproto: endfile */
/* aproto: file embed_graph/merge_queue_misc.c */
extern t_merge_queue embedg_merge_queue_new (int);
extern void embedg_merge_queue_delete (t_merge_queue);
extern boolean embedg_merge_queue_empty (t_merge_queue);
extern void embedg_merge_queue_print (t_merge_queue);
extern void embedg_merge_queue_append (t_merge_queue *, t_ver_edge *, int, int, int, int, int);
extern void embedg_merge_queue_append_vertex (t_merge_queue *, t_ver_edge *, int, int, int);
extern void embedg_merge_queue_append_virtual_vertex (t_merge_queue *, t_ver_edge *, int, int, int);
extern void embedg_merge_queue_get (t_merge_queue *, int *, int *, int *, int *);
extern void embedg_merge_queue_prune (t_merge_queue *, int *, int *, int *, int *);
/* aproto: endfile */
/* aproto: file embed_graph/obstruction.c */
extern void embedg_obstruction (t_ver_sparse_rep *, t_adjl_sparse_rep *, t_dlcl **, t_dlcl **, t_ver_edge *, int, int *, int, int, t_ver_sparse_rep **, t_adjl_sparse_rep **, int *);
extern minor embedg_mark_obstruction (t_dlcl **, t_dlcl **, t_ver_edge *, int, int *, int, int);
/* aproto: endfile */
/* aproto: file embed_graph/post_dfs_preproc.c */
extern t_dlcl *sparseg_order_wrt_lowpoint (int, int *, int *, t_dlcl *);
extern int *sparseg_find_least_ancestor (int, t_dlcl **);
extern int *sparseg_find_dfs_parent (int, t_dlcl **);
/* aproto: endfile */
/* aproto: file embed_graph/proper_face_walk.c */
extern boolean embedg_VES_get_succ_on_proper_face_with_avoidance (t_ver_edge *, int, int, int, int, boolean, int, int *, int *, int *);
extern void embedg_VES_get_succ_on_proper_face (t_ver_edge *, int, int, int, int, int, int *, int *, int *);
extern void embedg_VES_walk_proper_face (t_ver_edge *, int, int, int, boolean, int);
/* aproto: endfile */
/* aproto: file embed_graph/vertex_activity.c */
extern boolean embedg_VES_is_ver_pertinent (t_ver_edge *, int, int, int);
extern boolean embedg_VES_is_ver_ext_active (t_ver_edge *, int, int, int);
extern boolean embedg_VES_is_ver_int_active (t_ver_edge *, int, int, int);
extern boolean embedg_VES_is_ver_inactive (t_ver_edge *, int, int, int);
/* aproto: endfile */
/* aproto: file embed_graph/walkdown.c */
extern t_merge_queue embedg_walkdown (t_ver_edge *, int, int *, int);
/* aproto: endfile */
/* aproto: file embed_graph/sparseg_dlcl_misc.c */
extern void sparseg_dlcl_delete (t_dlcl **, int);
extern void sparseg_dlcl_print (t_dlcl **, int);
extern boolean sparseg_dlcl_is_adjacent (t_dlcl **, int, int, int, t_dlcl **);
extern void sparseg_dlcl_append_to_neigh_list (t_dlcl **, int, int, int, int);
extern void sparseg_dlcl_to_sparseg (t_dlcl **, int, int, t_ver_sparse_rep **, t_adjl_sparse_rep **);
extern boolean sparseg_dlcl_sub (t_dlcl **, int, t_dlcl **, int);
/* aproto: endfile */
/* aproto: file embed_graph/walkup.c */
extern void embedg_walkup (t_ver_edge *, int, int, t_dlcl *);
/* aproto: endfile */
/* aproto: file embed_graph/dfs_preprocessing.c */
extern void sparseg_adjl_dfs_preprocessing (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int *, int **, int **, int **, t_dlcl ***, t_dlcl ***, int **, int **, t_dlcl ***);
/* aproto: endfile */
/* aproto: file embed_graph/planar_by_edge_addition.c */
extern boolean sparseg_adjl_is_planar (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int *, t_dlcl ***, t_dlcl ***, t_dlcl ***, t_ver_edge **, int *, int *, int *);
/* aproto: endfile */
/* aproto: file embed_graph/recover.c */
extern void embedg_recover_embedding (t_ver_sparse_rep *, t_adjl_sparse_rep *, t_ver_edge *, int, int, t_dlcl **, t_ver_sparse_rep **, t_embed_sparse_rep **);
extern void embedg_recov_embed_walk_proper_face (int, int, t_adjl_sparse_rep *, t_embed_sparse_rep *, boolean, int);
extern boolean embedg_check_recov_embedding (int, int, int, t_ver_sparse_rep *, t_adjl_sparse_rep *, t_embed_sparse_rep *);
extern t_dlcl **embedg_recover_obstruction (t_ver_edge *, int, minor, int *);
extern boolean embedg_check_recov_obs (t_dlcl **, int, minor);
/* aproto: endfile */
/* aproto: file embed_graph/planar_alg_init.c */
extern t_ver_edge *embedg_planar_alg_init (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int *, int *, t_dlcl ***, t_dlcl ***, t_dlcl ***);
/* aproto: endfile */
/* aproto: file embed_graph/dlcl_misc.c */
extern t_dlcl *embedg_dlcl_rec_new (int);
extern void embedg_dlcl_rec_print (t_dlcl *);
extern void embedg_dlcl_print (t_dlcl *);
extern t_dlcl *embedg_dlcl_rec_append (t_dlcl *, t_dlcl *);
extern t_dlcl *embedg_dlcl_rec_prepend (t_dlcl *, t_dlcl *);
extern t_dlcl *embedg_dlcl_cat (t_dlcl *, t_dlcl *);
extern t_dlcl *embedg_dlcl_find (t_dlcl *, int);
extern t_dlcl *embedg_dlcl_find_with_NIL_twin_in_adjl (t_dlcl *, int);
extern t_dlcl *embedg_dlcl_delete_first (t_dlcl *);
extern t_dlcl *embedg_dlcl_delete_rec (t_dlcl *, t_dlcl *);
extern boolean embedg_dlcl_is_empty (t_dlcl *);
extern t_dlcl *embedg_dlcl_list_next (t_dlcl *);
extern t_dlcl *embedg_dlcl_list_prev (t_dlcl *);
extern t_dlcl *embedg_dlcl_list_last (t_dlcl *);
extern void embedg_dlcl_delete (t_dlcl *);
extern t_dlcl *embedg_dlcl_copy (t_dlcl *);
extern int embedg_dlcl_length (t_dlcl *);
/* aproto: endfile */
/* aproto: file embed_graph/sparseg_adjl.c */
extern boolean sparseg_adjl_plan_and_iso (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, int *, t_ver_sparse_rep **, t_adjl_sparse_rep **, t_embed_sparse_rep **, int *);
extern int *sparseg_adjl_footprint (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int);
extern void sparseg_adjl_print (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, boolean);
extern void sparseg_adjl_embed_print (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, t_embed_sparse_rep *, boolean);
extern graph *sparseg_adjl_to_nauty_graph (t_ver_sparse_rep *, int, t_adjl_sparse_rep *);
extern t_edge_sparse_rep *sparseg_adjl_edges (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, boolean);
/* aproto: endfile */
/* aproto: file embed_graph/embedding.c */
extern void embedg_embedding (t_ver_sparse_rep *, t_adjl_sparse_rep *, t_ver_edge *, int, int, int, int, t_dlcl **, t_ver_sparse_rep **, t_embed_sparse_rep **);
extern void embedg_remove_SCE (t_ver_edge *, int, int);
extern int *embedg_vertices_orientation (t_ver_edge *, int);
extern int embedg_merge_remaining_virtual (t_ver_edge *, int);
extern int embedg_nbr_faces (t_ver_edge *, int, int, int *, int *);
extern boolean embedg_is_embed_valid (t_ver_edge *, int, int, int, int *, int *);
/* aproto: endfile */
/* aproto: file embed_graph/sparseg_adjl_modify.c */
extern boolean sparseg_adjl_add_edge (t_ver_sparse_rep *, int, t_adjl_sparse_rep **, int *, int *, int, int, boolean);
extern boolean sparseg_adjl_add_edge_no_extend (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, int *, int, int, boolean);
extern boolean sparseg_adjl_add_dir_edge (t_ver_sparse_rep *, int, t_adjl_sparse_rep **, int *, int *, int, int, boolean);
extern boolean sparseg_adjl_add_dir_edge_no_extend (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, int *, int, int, boolean);
extern boolean sparseg_adjl_remove_edge_no_red (t_ver_sparse_rep *, t_adjl_sparse_rep *, int, int);
extern boolean sparseg_adjl_remove_dir_edge_no_red (t_ver_sparse_rep *, t_adjl_sparse_rep *, int, int);
extern int sparseg_adjl_remove_all_dir_edge_no_red (t_ver_sparse_rep *, t_adjl_sparse_rep *, int, int);
extern void sparseg_adjl_add_vertices (t_ver_sparse_rep **, int, int);
extern void sparseg_adjl_add_vertices_no_extend (t_ver_sparse_rep *, int, int);
extern void sparseg_adjl_remove_vertex (t_ver_sparse_rep **, int, t_adjl_sparse_rep *, int, int, int *);
extern void sparseg_adjl_remove_vertex_no_red (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, int *);
extern void sparseg_adjl_relabel_vertex (t_adjl_sparse_rep *, int, int);
/* aproto: endfile */
/* aproto: file embed_graph/sparseg_adjl_misc.c */
extern void sparseg_adjl_assign_V (t_ver_sparse_rep *, t_ver_sparse_rep *, int);
extern t_ver_sparse_rep *sparseg_adjl_dup_V (t_ver_sparse_rep *, int);
extern void sparseg_adjl_assign_A (t_adjl_sparse_rep *, t_adjl_sparse_rep *, int);
extern t_adjl_sparse_rep *sparseg_adjl_dup_A (t_adjl_sparse_rep *, int);
extern void sparseg_embed_assign_E (t_embed_sparse_rep *, t_embed_sparse_rep *, int);
extern t_embed_sparse_rep *sparseg_embed_dup_E (t_embed_sparse_rep *, int);
extern void sparseg_adjl_underlying_undir (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, t_ver_sparse_rep **, t_adjl_sparse_rep **, int *, int *);
extern void sparseg_adjl_edge_union (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, t_ver_sparse_rep *, t_adjl_sparse_rep *, boolean, t_ver_sparse_rep **, t_adjl_sparse_rep **, int *, int *);
extern void sparseg_adjl_add_edges (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, int *, t_ver_sparse_rep *, t_adjl_sparse_rep *, boolean, int *);
/* aproto: endfile */
/* aproto: file embed_graph/degree.c */
extern int *sparseg_adjl_degree_seq (t_ver_sparse_rep *, int, t_adjl_sparse_rep *);
/* aproto: endfile */
/* aproto: file embed_graph/neighbours.c */
extern void sparseg_adjl_vertex_out_neighbours (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, int **, int *);
extern void sparseg_adjl_vertex_in_neighbours (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, int **, int *);
/* aproto: endfile */
/* aproto: file embed_graph/misc.c */
extern void sparseg_comp_delete (t_comp_sparse_rep *, int);
extern void sparseg_comp_print (t_comp_sparse_rep *, int);
/* aproto: endfile */
/* aproto: file embed_graph/dfs_things.c */
extern void sparseg_adjl_dfs_things (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, boolean, int *, t_ver_sparse_rep **, t_adjl_sparse_rep **, int *, int **, int **, int **, int **, int **, int *, int **, int *, t_comp_sparse_rep **);
extern void sparseg_adjl_dfs_add_vertex_to_comp (t_comp_sparse_rep **, int *, int);
extern void sparseg_adjl_dfs_create_new_comp (t_comp_sparse_rep **, int *);
extern void sparseg_adjl_from_comps_to_sparseg_adjl (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, t_comp_sparse_rep *, int, t_graph_sparse_rep **);
extern void sparseg_adjl_from_start_comp_to_sparseg_adjl (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int *, t_graph_sparse_rep **);
/* aproto: endfile */
/* aproto: file embed_graph/dfs_for_digraph.c */
extern void sparseg_adjl_dfs_digraph (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int, int, t_ver_sparse_rep **, t_adjl_sparse_rep **, int *, int **, int **, int **, int **, int **, int *, t_comp_sparse_rep **);
/* aproto: endfile */
/* aproto: file embed_graph/bfs.c */
extern void sparseg_adjl_bfs (t_ver_sparse_rep *, int, t_adjl_sparse_rep *, boolean, int, t_ver_sparse_rep **, t_adjl_sparse_rep **, int *);
/* aproto: endfile */
/* aproto: file embed_graph/faces.c */
extern void sparseg_adjl_walk_proper_face (int, t_adjl_sparse_rep *, t_embed_sparse_rep *, int, boolean, int, int *, t_edge_sparse_rep *);
extern void sparseg_adjl_get_face_edges (int, t_adjl_sparse_rep *, int, t_embed_sparse_rep *, int, t_edge_sparse_rep **, int *, boolean, int);
/* aproto: endfile */

#ifdef __cplusplus
}
#endif

#endif  /* _PLANARITY_H_ */
