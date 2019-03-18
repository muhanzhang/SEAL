#!/bin/bash

cd ../../
git clone https://github.com/muhanzhang/pytorch_DGCNN
cd pytorch_DGCNN
cd lib
make -j4
cd "$(dirname "$0")"
pip install --user numpy
pip install --user scipy
pip install --user networkx
pip install --user tqdm
pip install --user sklearn
pip install --user gensim
