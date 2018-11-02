import numpy as np
import scipy.io as scio
from scipy import sparse
import argparse


parse = argparse.ArgumentParser()
parse.add_argument("-g", "--graph", help="each row is a training link")
parse.add_argument("-u", "--users", help="users list, each row is a node")
parse.add_argument("-n", "--dimension", help="dimension of node profile")
parse.add_argument("-f", "--profile", help="profile information of nodes, corresponding to node list")
parse.add_argument("-o", "--output", help="output .mat file")
args = parse.parse_args()


class MatFile:
    def __init__(self):
        pass

    @staticmethod
    def _transform_graph(self) -> object:
        data = np.loadtxt(args.graph, delimiter=",", dtype=np.int64)
        net = sparse.dok_matrix((data.max() + 1, data.max() + 1))
        for line in data:
            net[line[0], line[1]] = 1.0
        return net

    @staticmethod
    def _transform_feature(self) -> object:
        users = np.loadtxt(args.users, dtype=np.int64)
        features = np.loadtxt(args.profile, dtype=np.int64)
        group = sparse.dok_matrix((users.max() + 1 , int(args.dimension)))
        for index, user in enumerate(users):
            for f_id, feature in enumerate(features[index]):
                group[user, f_id] = feature
        return group

    def create_mat_file(self, file_path):
        net = self._transform_graph()
        group = self._transform_feature()
        scio.savemat(file_path, {"net": net, "group": group})

    @staticmethod
    def load_mat_file(self, file_path):
        data = scio.loadmat(file_path)
        for key, value in data.items():
            print(key)
            print(value)

if __name__ == "__main__":
    mat_file = MatFile()
    mat_file.create_mat_file(args.output)
    # mat_file.load_mat_file(args.output)