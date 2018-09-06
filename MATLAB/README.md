SEAL -- learning from Subgraphs, Embeddings, and Attributes for Link prediction
===============================================================================

About
-----

MATLAB version of SEAL (learning from Subgraphs, Embeddings, and Attributes for Link prediction). This repo also contains implementations of other link prediction baselines (heuristics, latent feature methods, embedding methods, WLK and WLNM).

How to run
----------

Please download our [\[DGCNN software\]](https://github.com/muhanzhang/DGCNN) to the same level as the root SEAL folder (not this MATLAB folder). DGCNN is the default graph neural network in SEAL.

Install the required libraries of DGCNN. Then run "Main.m" in MATLAB to do the experiments. By default, it will run SEAL on the USAir dataset for 1 time. Modify variables such as _numOfExperiment_, _ratioTrain_, _dataname_, _method_ in the setting part of "Main.m" to change experimental settings.

Results will be saved in "data/result".

Requirements
------------

### Requirements for SEAL

Please follow its instruction to install DGCNN. To generate Torch-readable graphs, Torch library matio is required to be installed:

    1# OSX
    brew install homebrew/science/libmatio
    2# Ubuntu
    sudo apt-get install libmatio2
    luarocks install --local matio

To calculate AUC within DGCNN, install Torch library metrics:

    git clone https://github.com/hpenedones/metrics.git
    cd metrics
    luarocks --local make

MATLAB toolbox Bioinformatics is required to calculate graph shortest path distance. 

Please install the network embedding software [\[node2vec\]](https://github.com/aditya-grover/node2vec) into "software/".

If you cannot call Torch within MATLAB, please refer to this [\[README\]](https://github.com/muhanzhang/LinkPrediction) for fixings.

### Requirements for WLK

The baseline WLK (Weisfeiler-Lehman graph kernel) requires installing libsvm. We include libsvm-3.22 already in "software/" folder. Type:

    cd software/libsvm-3.22
    make

to compile libsvm on your unix machine.

### Requirements for WLNM

The baseline WLNM (Weisfeiler-Lehman Neural Machine) has the following requirements: liblinear, nauty, nnsparse, svm, metircs. Liblinear is used to save .mat to sparse libsvm format, so that Torch neural networks can read them. Type:

    cd software/liblinear-2.1
    make

to compile liblinear on your unix machine. 

A graph canonization software Nauty is required to break ties in WL labelings:

    cd software/nauty26r7
    ./configure
    make

Torch libaries nnsparse, svm are required in the neural network:

    luarocks install --local nnsparse
    luarocks install --local svm

### Requirements for embedding methods

Two network embedding software: node2vec and LINE, have been included in "software/". If they do not work, you may need to reinstall them from source. To run embedding-based link prediction, liblinear is required. After compiling liblinear, cd to "software/liblinear-2.1/matlab/" in MATLAB and type "make" to install its MATLAB interface. Please also change the names of "train.mexa64" and "predict.mexa64" to "liblinear_train.mexa64" and "liblinear_predict.mexa64" in "software/liblinear-2.1/matlab/", respectively.

Reference
---------

If you find the code useful, please cite our paper:

    @article{zhang2018link,
      title={Link Prediction Based on Graph Neural Networks},
      author={Zhang, Muhan and Chen, Yixin},
      journal={arXiv preprint arXiv:1802.09691},
      year={2018}
    }

Muhan Zhang, Washington University in St. Louis
muhan@wustl.edu
2/10/2018
