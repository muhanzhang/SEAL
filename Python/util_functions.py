from __future__ import print_function
import numpy as np
import random
from tqdm import tqdm
import os, sys, pdb, math, time
#import cPickle as cp
import _pickle as cp  # python3 compatability
import networkx as nx
import argparse
import scipy.io as sio
import scipy.sparse as ssp
from sklearn import metrics
from gensim.models import Word2Vec
import warnings
warnings.simplefilter('ignore', ssp.SparseEfficiencyWarning)
cur_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append('%s/../../pytorch_DGCNN' % cur_dir)
sys.path.append('%s/software/node2vec/src' % cur_dir)
from util import GNNGraph, HGNNHypergraph,HGNNBihypergraph
import node2vec
import multiprocessing as mp
from collections import defaultdict
from scipy.sparse import csr_matrix,hstack, vstack
sys.path.append('/home2/e1-313-15477/hynetworkx/src')
from data_preparer import clean_train_hypergraph



class ForkedPdb(pdb.Pdb):
    """A Pdb subclass that may be used
    from a forked multiprocessing child

    """
    def interaction(self, *args, **kwargs):
        _stdin = sys.stdin
        try:
            sys.stdin = open('/dev/stdin')
            pdb.Pdb.interaction(self, *args, **kwargs)
        finally:
            sys.stdin = _stdin

def sample_neg(net, test_ratio=0.1, train_pos=None, test_pos=None, max_train_num=None):
    # get upper triangular matrix
    net_triu = ssp.triu(net, k=1)
    # sample positive links for train/test
    row, col, _ = ssp.find(net_triu)
    # sample positive links if not specified
    if train_pos is None or test_pos is None:
        perm = random.sample(range(len(row)), len(row))
        row, col = row[perm], col[perm]
        split = int(math.ceil(len(row) * (1 - test_ratio)))
        train_pos = (row[:split], col[:split])
        test_pos = (row[split:], col[split:])
    # if max_train_num is set, randomly sample train links
    if max_train_num is not None:
        perm = np.random.permutation(len(train_pos[0]))[:max_train_num]
        train_pos = (train_pos[0][perm], train_pos[1][perm])
    # sample negative links for train/test
    train_num, test_num = len(train_pos[0]), len(test_pos[0])
    neg = ([], [])
    n = net.shape[0]
    print('sampling negative links for train and test')
    while len(neg[0]) < train_num + test_num:
        i, j = random.randint(0, n-1), random.randint(0, n-1)
        if i < j and net[i, j] == 0:
            neg[0].append(i)
            neg[1].append(j)
        else:
            continue
    train_neg  = (neg[0][:train_num], neg[1][:train_num])
    test_neg = (neg[0][train_num:], neg[1][train_num:])
    return train_pos, train_neg, test_pos, test_neg

def sample_neg_bip(net, test_ratio=0.1, train_pos=None, test_pos=None, max_train_num=None):
    '''
    Note that net is NOT a symmetric matrix!
    '''
    # sample positive links for train/test
    row, col, _ = ssp.find(net)
    # sample positive links if not specified
    if train_pos is None or test_pos is None:
        perm = random.sample(range(len(row)), len(row))
        row, col = row[perm], col[perm]
        split = int(math.ceil(len(row) * (1 - test_ratio)))
        train_pos = (row[:split], col[:split])
        test_pos = (row[split:], col[split:])
    # if max_train_num is set, randomly sample train links
    if max_train_num is not None:
        perm = np.random.permutation(len(train_pos[0]))[:max_train_num]
        train_pos = (train_pos[0][perm], train_pos[1][perm])
    # sample negative links for train/test
    train_num, test_num = len(train_pos[0]), len(test_pos[0])
    neg = ([], [])
    m,n = net.shape
    print('sampling negative links for train and test')
    while len(neg[0]) < train_num + test_num:
        i, j = random.randint(0, m-1), random.randint(0, n-1)
        if net[i, j] == 0:
            neg[0].append(i)
            neg[1].append(j)
        else:
            continue
    train_neg  = (neg[0][:train_num], neg[1][:train_num])
    test_neg = (neg[0][train_num:], neg[1][train_num:])
    return train_pos, train_neg, test_pos, test_neg


def links2sub_bihypergraphs(S, S_, B, train_pos, train_neg, test_pos, test_neg, h=1, max_nodes_per_hop=None, 
                    node_information=None, edge_information=None):
    # automatically select h from {1, 2}
    if h == 'auto':
        # split train into val_train and val_test
        _, _, val_test_pos, val_test_neg = sample_neg(A, 0.1)
        val_A = A.copy()
        val_A[val_test_pos[0], val_test_pos[1]] = 0
        val_A[val_test_pos[1], val_test_pos[0]] = 0
        val_auc_CN = CN(val_A, val_test_pos, val_test_neg)
        val_auc_AA = AA(val_A, val_test_pos, val_test_neg)
        print('\033[91mValidation AUC of AA is {}, CN is {}\033[0m'.format(val_auc_AA, val_auc_CN))
        if val_auc_AA >= val_auc_CN:
            h = 2
            print('\033[91mChoose h=2\033[0m')
        else:
            h = 1
            print('\033[91mChoose h=1\033[0m')

    # extract enclosing subgraphs
    max_label_n = {'value': 0}
    max_label_n_ = {'value': 0}
#     pdb.set_trace()
    def helper(S, S_, B, links, g_label):
        '''
        g_list = []
        for i, j in tqdm(zip(links[0], links[1])):
            g, n_labels, n_features = subgraph_extraction_labeling((i, j), A, h, max_nodes_per_hop, node_information)
            max_n_label['value'] = max(max(n_labels), max_n_label['value'])
            g_list.append(GNNGraph(g, g_label, n_labels, n_features))
        return g_list
        '''
        # the new parallel extraction code
        start = time.time()
        pool = mp.Pool(mp.cpu_count())
        results = pool.map_async(parallel_worker, [((i, j), S, S_, B, h, max_nodes_per_hop, node_information, edge_information) for i, j in zip(links[0], links[1])])
        remaining = results._number_left
        pbar = tqdm(total=remaining)
        while True:
            pbar.update(remaining - results._number_left)
            if results.ready(): break
            remaining = results._number_left
            time.sleep(1)
        results = results.get()
        pool.close()
        pbar.close()
#         g_list = [GNNGraph(g, g_label, n_labels, n_features) for g, _, n_labels, n_features in results]
#         hyg_list = [HGNNHypergraph(S_g, g_label, n_labels, n_features, f_labels) for _, S_g, n_labels, f_labels, n_features in results]
#         pdb.set_trace()
        bihyg_list = [HGNNBihypergraph(sub_S, sub_S_, sub_B, g_label, labels_n, labels_n_, features_n, features_n_) for sub_S, sub_S_, sub_B, labels_n, labels_n_, features_n, features_n_ in results]
#         pdb.set_trace()
        max_label_n['value'] = max(max([max(labels_n) for _, _, _, labels_n, _, _, _ in results]) if len(results) > 0 else 0, max_label_n['value'])
        max_label_n_['value'] = max(max([max(labels_n_) for _, _, _, _, labels_n_, _, _ in results]) if len(results) > 0 else 0, max_label_n_['value'])
        
        end = time.time()
        print("Time eplased for subgraph extraction: {}s".format(end-start))
        return bihyg_list

    print('Enclosing subgraph extraction begins...')
#     train_graphs = helper(A, train_pos, 1) + helper(A, train_neg, 0)
#     test_graphs = helper(A, test_pos, 1) + helper(A, test_neg, 0)
    train_bihypergraphs = helper(S, S_, B, train_pos, 1) + helper(S, S_, B, train_neg, 0)
    test_bihypergraphs = helper(S, S_, B, test_pos, 1) + helper(S, S_, B, test_neg, 0)
    print(max_label_n, max_label_n_)
    return train_bihypergraphs, test_bihypergraphs, max_label_n['value'], max_label_n_['value']


def parallel_worker(x):
    return subgraph_extraction_labeling(*x)
    

def filter_hypergraph(S, nodes):
#     S_cleaned = clean_train_hypergraph(S, csr_matrix(([1, 1], ([nodes[0], nodes[1]], [nodes[1], nodes[0]])), shape=(S.shape[0], S.shape[0])))
    S_cleaned=S
    F_dict = incidence_to_hyperedges(S_cleaned, _type=dict)
    matching_F_dict = {j: f for j, f in F_dict.items() if f.issubset(set(nodes))}
    I = nodes
    J = list(matching_F_dict)
    try:
        S_g = S_cleaned[I, :][:, J]
    except IndexError:
        print('Caught exception', type(I), type(J), S_cleaned.shape)
        raise (IndexError)
    return S_g


def subgraph_extraction_labeling(ind, S, S_, B, h=1, max_nodes_per_hop=None,
                                 node_information=None, edge_information=None,
                                 node_information_=None):
    # extract the h-hop enclosing subgraph around link 'ind'
#     print(B.shape)
    dist = 0
    nodes = set([ind[0]])
    visited = set([ind[0]])
    fringe = set([ind[0]])
    
    nodes_ = set([ind[1]])
    visited_ = set([ind[1]])
    fringe_ = set([ind[1]])
    
    nodes_dist = [0, 0]
    nodes_dist_ = [0, 0]
#     pdb.set_trace()
    
    for dist in range(1, 3*h+1):
        
        fringe_, fringe = neighbors(fringe, B.T), neighbors(fringe_, B)
        fringe_ = fringe_ - visited_
        fringe = fringe - visited
        visited_ = visited_.union(fringe_)
        visited = visited.union(fringe)
        if max_nodes_per_hop is not None:
            if max_nodes_per_hop < len(fringe):
                fringe = random.sample(fringe, max_nodes_per_hop)
            if max_nodes_per_hop_ < len(fringe_):
                fringe_ = random.sample(fringe_, max_nodes_per_hop)
        if len(fringe_) == 0 or len(fringe) == 0:
            break
        nodes_ = nodes_.union(fringe_)
        nodes = nodes.union(fringe)
        
        nodes_dist_ += [dist] * len(fringe_)
        nodes_dist += [dist] * len(fringe)
    # move target nodes to top
    nodes.remove(ind[0])
    nodes_.remove(ind[1])
    nodes = [ind[0]] + list(nodes) 
    nodes_ = [ind[1]] + list(nodes_)
    
    sub_B = B[nodes, :][:, nodes_]
#     sub_S = filter_hypergraph(S.T, nodes).T
#     sub_S_ = filter_hypergraph(S_.T, nodes_).T
    sub_S = S[:, nodes]
    sub_S=sub_S[(sub_S.sum(axis=1)!=0).nonzero()[0], :]
    sub_S_ = S_[:, nodes_]
    sub_S_=sub_S_[(sub_S_.sum(axis=1)!=0).nonzero()[0], :]
#     print('ind: {}, sV: {}, sV_: {}, S: {}, S_: {}, B: {}'.format(ind, nodes, nodes_, sub_S.shape, sub_S_.shape, sub_B.shape))
    # apply node-labeling
#     print("Sub_s",sub_S.todense())
#     print("Sub_s_",sub_S_.todense())
    labels, labels_ = node_label(sub_S, sub_S_, sub_B)
    # get node features
    features = None
    features_ = None
    if node_information is not None:
        features = node_information[nodes]
        features_ = node_information_[nodes_]
#     if edge_information is not None:
#         features = edge_information[edges]
    # construct nx graph
#     g = nx.from_scipy_sparse_matrix(subgraph)
    # remove link between target nodes
#     if g.has_edge(0, 1):
#         g.remove_edge(0, 1)
    # TODO: remove edge 0, 1 from S_g also
   
    return sub_S, sub_S_, sub_B, labels.tolist(), labels_.tolist(), features, features_


def neighbors(fringe, A):
    # find all 1-hop neighbors of nodes in fringe from A
    res = set()
    for node in fringe:
        nei, _, _ = ssp.find(A[:, node])
        nei = set(nei)
        res = res.union(nei)
    return res


def S_to_bipartite(S):
    B = nx.bipartite.from_biadjacency_matrix(S)
    return B

def get_Z(n_rows, n_cols):
    return csr_matrix(([], ([], [])), shape=(n_rows, n_cols))

def stack_quadrant(tl, tr, bl, br):
    if (tl is None and tr is None) or (bl is None and br is None) or \
       (tl is None and bl is None) or (tr is None and br is None):
        print('Warning: Unstackable! Size of zero matrices not known.')
        return None
    if tl is None:
        tl = get_Z(tr.shape[0], bl.shape[1])
    if tr is None:
        tr = get_Z(tl.shape[0], br.shape[1])
    if bl is None:
        bl = get_Z(br.shape[0], tl.shape[1])
    if br is None:
        br = get_Z(bl.shape[0], tr.shape[1])
    l = vstack([tl, bl])
    r = vstack([tr, br])
    return hstack([l, r]).tocsr()

def untetrapartite(S, S_, B):
    # Return an (n'+m+m'+n) x (n'+m+m'+n) matrix, in that order
    n, m = S.shape
    n_, m_ = S_.shape
    rows_n_, rows_m, rows_m_, rows_n = map(list,(range(n_), range(n_, n_+m), range(n_+m, n_+m+m_), range(n_+m+m_, n_+m+m_+n)))
    cols_m_, cols_n, cols_n_, cols_m = map(list,(range(m_), range(m_, m_+n), range(m_+n, m_+n+n_), range(m_+n+n_, m_+n+n_+m)))
    tl = stack_quadrant(S_, None, B, S.T)
    A = stack_quadrant(tl, None, None, tl.T)
    A = A[:, cols_n_+cols_m+cols_m_+cols_n]
    return A, rows_n, rows_n_, rows_m, rows_m_

def node_label(sub_S, sub_S_, sub_B):
    # an implementation of the proposed double-radius node labeling (DRNL)
    subgraph, rows_n, rows_n_, rows_m, rows_m_ = untetrapartite(sub_S, sub_S_, sub_B)
#     print(subgraph.todense())
#     print('S: {}, S_: {}, B: {}\t n: {}, n_: {}, m: {}, m_: {}, subg: {}'.format(sub_S.shape, sub_S_.shape, sub_B.shape, rows_n, rows_n_, rows_m, rows_m_, subgraph.shape))
    K = subgraph.shape[0]
    u = rows_m[0]
    v = rows_m_[0]
    u, v = min([u, v]), max([u, v])
    nodes_wou = list(sorted(set(range(K)).difference({u})))
    nodes_wov = list(sorted(set(range(K)).difference({v})))
    subgraph_wou = subgraph[nodes_wou, :][:, nodes_wou]
    subgraph_wov = subgraph[nodes_wov, :][:, nodes_wov]
    dist_to_v = ssp.csgraph.shortest_path(subgraph_wou, directed=False, unweighted=True)
    dist_to_v = dist_to_v[list(sorted(set(range(dist_to_v.shape[0])).difference({v-1}))), v-1]
    dist_to_u = ssp.csgraph.shortest_path(subgraph_wov, directed=False, unweighted=True)
    dist_to_u = dist_to_u[list(sorted(set(range(dist_to_u.shape[0])).difference({u}))), u]
    d = (dist_to_u + dist_to_v).astype(int)
    d_over_2, d_mod_2 = np.divmod(d, 2)
    labels = 1 + np.minimum(dist_to_u, dist_to_v).astype(int) + d_over_2 * (d_over_2 + d_mod_2 - 1)
    
    labels = np.concatenate((labels[:u], np.array([1]), labels[u:v], np.array([1]), labels[v:]))
    labels[np.isinf(labels)] = 0
    labels[labels>1e6] = 0  # set inf labels to 0
    labels[labels<-1e6] = 0  # set -inf labels to 0
#     print(labels, rows_n, rows_n_)
    labels_ = labels[rows_n_]
    labels = labels[rows_n]
    return labels, labels_


def generate_node2vec_embeddings(A, emd_size=128, negative_injection=False, train_neg=None):
    if negative_injection:
        row, col = train_neg
        A = A.copy()
        A[row, col] = 1  # inject negative train
        A[col, row] = 1  # inject negative train
    nx_G = nx.from_scipy_sparse_matrix(A)
    G = node2vec.Graph(nx_G, is_directed=False, p=1, q=1)
    G.preprocess_transition_probs()
    walks = G.simulate_walks(num_walks=10, walk_length=80)
    walks = [map(str, walk) for walk in walks]
    model = Word2Vec(walks, size=emd_size, window=10, min_count=0, sg=1, 
            workers=8, iter=1)
    wv = model.wv
    embeddings = np.zeros([A.shape[0], emd_size], dtype='float32')
    sum_embeddings = 0
    empty_list = []
    for i in range(A.shape[0]):
        if str(i) in wv:
            embeddings[i] = wv.word_vec(str(i))
            sum_embeddings += embeddings[i]
        else:
            empty_list.append(i)
    mean_embedding = sum_embeddings / (A.shape[0] - len(empty_list))
    embeddings[empty_list] = mean_embedding
    return embeddings

def incidence_to_hyperedges(S, silent_mode=True, _type=set):
    I, J = S.nonzero()
    hyperedges = defaultdict(set)
    indices = list(zip(I, J))
    if not silent_mode:
        print('Converting incidence matrix to hyperedge {} for faster processing...'.format(_type))
    for i, j in (tqdm(indices) if not silent_mode else indices):
        hyperedges[j].add(i)
    if _type == set:
        return set(map(frozenset, hyperedges.values()))
    elif _type == list:
        return set(map(frozenset, hyperedges.values()))
    elif _type == dict:
        return {i: frozenset(f) for i, f in hyperedges.items()}
    return hyperedges


def hyperedges_to_incidence(hyperedges, nV):
    nF = len(hyperedges)
    hyperedges = list(set(hyperedges))
    I = []
    J = []
    for j, f in enumerate(hyperedges):
        I.extend(f)
        J.extend([j] * len(f))
    S = csr_matrix(([1] * len(I), (I, J)), shape=(nV, nF))
    return S


def AA(A, test_pos, test_neg):
    # Adamic-Adar score
    A_ = A / np.log(A.sum(axis=1))
    A_[np.isnan(A_)] = 0
    A_[np.isinf(A_)] = 0
    sim = A.dot(A_)
    return CalcAUC(sim, test_pos, test_neg)
    
        
def CN(A, test_pos, test_neg):
    # Common Neighbor score
    sim = A.dot(A)
    return CalcAUC(sim, test_pos, test_neg)


def CalcAUC(sim, test_pos, test_neg):
    pos_scores = np.asarray(sim[test_pos[0], test_pos[1]]).squeeze()
    neg_scores = np.asarray(sim[test_neg[0], test_neg[1]]).squeeze()
    scores = np.concatenate([pos_scores, neg_scores])
    labels = np.hstack([np.ones(len(pos_scores)), np.zeros(len(neg_scores))])
    fpr, tpr, _ = metrics.roc_curve(labels, scores, pos_label=1)
    auc = metrics.auc(fpr, tpr)
    return auc

import pandas as pd
import os
def load_bipartite_hypergraph(data_params):
    id_p_map = pd.read_csv(os.path.join(data_params['home_path'], data_params['r_label_file']), sep='\t', header=None)
    id_a_map = pd.read_csv(os.path.join(data_params['home_path'], data_params['u_label_file']), sep='\t', header=None)
    id_a_map = dict(zip(id_a_map[0], id_a_map[1]))
    id_k_map = pd.read_csv(os.path.join(data_params['home_path'], data_params['v_label_file']), sep='\t', header=None)
    id_k_map = dict(zip(id_k_map[0], id_k_map[1]))
    p_a_list_map = pd.read_csv(os.path.join(data_params['home_path'], data_params['r_u_list_file']), sep=':', header=None)
    p_k_list_map = pd.read_csv(os.path.join(data_params['home_path'], data_params['r_v_list_file']), sep=':', header=None)
    n_p, na, nk = len(id_p_map), len(id_a_map), len(id_k_map)
    pos_A = list(map(lambda x: list(map(int, x.split(','))), p_a_list_map[1]))
    pos_B = list(map(lambda x: list(map(int, x.split(','))), p_k_list_map[1]))    
    # I, J, V: row, col, value of author-hypergraph
    # I_, J_, V_: row, col, value of keyword-hypergraph
    # I_B, J_B, V_B: row, col, value of author_hyperedge-keyword_hyperedge link
    I=[]
    J=[]
    V=[]
    I_=[]
    J_=[]
    V_=[]

    I_B=[]
    J_B=[]
    V_B=[]
    U_set=set()
    V_set=set()
    u_map={}
    v_map={}
    j_u=-1
    j_v=-1
    for u,v in zip(pos_A,pos_B):
        u=frozenset(u)
        v=frozenset(v)

        if u not in U_set:
            j_u+=1
            U_set.add(u)
            u_map[u]=j_u
            I.extend(list(u))
            J.extend([j_u]*len(u))
            V.extend([1]*len(u))
        if v not in V_set:
            j_v+=1
            V_set.add(v)
            v_map[v]=j_v
            I_.extend(list(v))
            J_.extend([j_v]*len(v))
            V_.extend([1]*len(v))

        I_B.append(u_map[u])
        J_B.append(v_map[v])
        V_B.append(1)

    n=max(I)+1
    m=len(U_set)
    n_=max(I_)+1
    m_=len(V_set)
    S = csr_matrix((V, (I, J)), shape=(n, m))
    S_ = csr_matrix((V_, (I_, J_)), shape=(n_, m_))
    B = csr_matrix((V_B, (I_B, J_B)), shape=(m, m_))    
    return S,S_,B





