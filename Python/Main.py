import torch
import pdb
import numpy as np
import sys, copy, math, time, pdb
# import cPickle as pickle
import scipy.io as sio
import scipy.sparse as ssp
import os.path

import argparse
sys.path.append('%s/../../pytorch_DGCNN' % os.path.dirname(os.path.realpath(__file__)))
from main import *
from util_functions import *
import random
from joblib import Memory

data_path = '/storage2/home2/e1-313-15477/swyam/SEAL/Python'
cachedir = os.path.join(data_path, 'cache')
memory = Memory(cachedir, verbose=0)

def run_parser():
    parser = argparse.ArgumentParser(description='Link Prediction with SEAL')
    # general settings
    parser.add_argument('--data-name', default='USAir', help='network name')
    parser.add_argument('--train-name', default=None, help='train name')
    parser.add_argument('--test-name', default=None, help='test name')
    parser.add_argument('--max-train-num', type=int, default=100000, 
                        help='set maximum number of train links (to fit into memory)')
    parser.add_argument('--no-cuda', action='store_true', default=False,
                        help='disables CUDA training')
    parser.add_argument('--seed', type=int, default=1, metavar='S',
                        help='random seed (default: 1)')
    parser.add_argument('--test-ratio', type=float, default=0.5,
                        help='ratio of test links')
    # model settings
    parser.add_argument('--hop', default=1, metavar='S', 
                        help='enclosing subgraph hop number, \
                        options: 1, 2,..., "auto"')
    parser.add_argument('--max-nodes-per-hop', default=None, 
                        help='if > 0, upper bound the # nodes per hop by subsampling')
    parser.add_argument('--use-embedding', action='store_true', default=False,
                        help='whether to use node2vec node embeddings')
    parser.add_argument('--use-edge-embedding', action='store_true', default=False,
                        help='whether to use node2vec edge embeddings')
    parser.add_argument('--use-attribute', action='store_true', default=False,
                        help='whether to use node attributes')
    parser.add_argument('--use-edge-attribute', action='store_true', default=False,
                        help='whether to use edge attributes')

    args = parser.parse_args('')
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    print(torch.cuda.is_available())
    torch.manual_seed(args.seed)
    if args.cuda:
        torch.cuda.manual_seed(args.seed)
    print(args)

    random.seed(cmd_args.seed)
    np.random.seed(cmd_args.seed)
    torch.manual_seed(cmd_args.seed)
    if args.hop != 'auto':
        args.hop = int(args.hop)
    if args.max_nodes_per_hop is not None:
        args.max_nodes_per_hop = int(args.max_nodes_per_hop)


    '''Prepare data'''
    args.file_dir = os.path.dirname(os.path.realpath('__file__'))
    args.res_dir = os.path.join(args.file_dir, 'results/{}'.format(args.data_name))

return args


@memory.cache
def initialize_1(data_params):
    S,S_,B = load_bipartite_hypergraph(data_params)
    return S,S_,B


@memory.cache
def initialize_2(B,test_ratio,max_train_num):
    my_net = B
    train_pos, train_neg, test_pos, test_neg = sample_neg_bip(my_net, test_ratio, max_train_num)
    return my_net, train_pos, train_neg, test_pos, test_neg

@memory.cache
def initialize_3(S, S_,A, train_pos, train_neg, test_pos, test_neg,
                 node_information, edge_information,hop,max_nodes_per_hop ):
    train_hypergraphs, test_hypergraphs, max_n_label, max_f_label = links2sub_bihypergraphs(S, S_,A, train_pos, 
    train_neg, test_pos, test_neg, hop, max_nodes_per_hop, node_information, edge_information)
    return train_hypergraphs, test_hypergraphs, max_n_label, max_f_label



''' Reading BiHypergraph Data '''
# home_path = 'sample_data/'
home_path = 'main_data/'
data_params = {'home_path': home_path,
               'r_label_file': 'id_p_map.txt',
               'u_label_file': 'id_a_map.txt',
               'v_label_file': 'id_k_map.txt',
               'r_u_list_file': 'p_a_list_train.txt',
               'r_v_list_file': 'p_k_list_train.txt',
               'emb_pkl_file': 'nodevectors.pkl'}

test_ratio = args.test_ratio
max_train_num = args.max_train_num 
node_information = None
edge_information = None
hop = args.hop
max_nodes_per_hop = args.max_nodes_per_hop

S, S_, B = initialize_1(data_params)
my_net, train_pos, train_neg, test_pos, test_neg = initialize_2(B, test_ratio,max_train_num)
A = my_net.copy()  # the observed network
A[test_pos[0], test_pos[1]] = 0  # mask test links
train_hypergraphs, test_hypergraphs, max_n_label, max_f_label = initialize_3(S, S_, A, train_pos, train_neg, test_pos, test_neg, node_information, edge_information, hop, max_nodes_per_hop)

print('# train: %d, # test: %d' % (len(train_hypergraphs), len(test_hypergraphs)))

# DGCNN configurations
cmd_args.gm = 'DGCNN'
cmd_args.sortpooling_k = 0.6
cmd_args.latent_dim = [64, 32, 16, 1]
cmd_args.hidden = 128
cmd_args.out_dim = 0
cmd_args.dropout = True
# cmd_args.dropout_prob = 0.9
cmd_args.num_class = 2
cmd_args.mode = 'gpu' if args.cuda else 'cpu'
cmd_args.num_epochs = 40
cmd_args.learning_rate = 1e-4
cmd_args.batch_size = 50
cmd_args.printAUC = True
cmd_args.feat_dim = max_n_label + 1
cmd_args.edge_feat_dim = 0 #max_f_label + 1
cmd_args.attr_dim = 0
cmd_args.edge_attr_dim = 0

if node_information is not None:
    cmd_args.attr_dim = node_information.shape[1]
    
if edge_information is not None:
    cmd_args.edge_attr_dim = edge_information.shape[1]
    
if cmd_args.sortpooling_k <= 1:
    num_nodes_list = sorted([g.num_nodes for g in train_hypergraphs + test_hypergraphs])
    cmd_args.sortpooling_k = num_nodes_list[int(math.ceil(cmd_args.sortpooling_k * len(num_nodes_list))) - 1]
    cmd_args.sortpooling_k = max(10, cmd_args.sortpooling_k)
    print('k used in SortPooling is: ' + str(cmd_args.sortpooling_k))

#pdb.set_trace()
classifier = Classifier()
if cmd_args.mode == 'gpu':
    classifier = classifier.cuda()

optimizer = optim.Adam(classifier.parameters(), lr=cmd_args.learning_rate)

train_idxes = list(range(len(train_hypergraphs)))
best_loss = None
for epoch in range(cmd_args.num_epochs):
    random.shuffle(train_idxes)
    classifier.train()
    avg_loss = loop_dataset(train_hypergraphs, classifier, train_idxes, optimizer=optimizer)
    if not cmd_args.printAUC:
        avg_loss[2] = 0.0
    print('\033[92maverage training of epoch %d: loss %.5f acc %.5f auc %.5f\033[0m' % (epoch, avg_loss[0], avg_loss[1], avg_loss[2]))

    classifier.eval()
    pdb.set_trace()
    test_loss = loop_dataset(test_hypergraphs, classifier, list(range(len(test_hypergraphs))))
    if not cmd_args.printAUC:
        test_loss[2] = 0.0
    print('\033[93maverage test of epoch %d: loss %.5f acc %.5f auc %.5f\033[0m' % (epoch, test_loss[0], test_loss[1], test_loss[2]))

with open('acc_results.txt', 'a+') as f:
    f.write(str(test_loss[1]) + '\n')

if cmd_args.printAUC:
    with open('auc_results.txt', 'a+') as f:
        f.write(str(test_loss[2]) + '\n')
