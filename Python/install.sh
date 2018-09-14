#!/bin/bash

cd ../../
git clone https://github.com/muhanzhang/pytorch_DGCNN
cd pytorch_DGCNN-master
unzip pytorch_structure2vec-master.zip
cd pytorch_structure2vec-master/s2vlib/
make -j4
cd SEAL/Python


