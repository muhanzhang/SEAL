#!/bin/bash

cd ../../
git clone https://github.com/muhanzhang/pytorch_DGCNN
cd pytorch_DGCNN
unzip pytorch_structure2vec-master.zip
cd pytorch_structure2vec-master/s2v_lib/
make -j4
cd "$(dirname "$0")"
pip install --user numpy
pip install --user scipy
pip install --user networkx
pip install --user tqdm
pip install --user sklearn
pip install --user gensim
