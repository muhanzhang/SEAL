Weighted Stochastic Blockmodel

This tool box is an implementation of our Bayesian variational algorithm 
for inferring community structure in weighted networks using the 
Weighted Stochastic Block Model (WSBM) based on our papers written by 
Christopher Aicher, Abigail Z. Jacobs, and Aaron Clauset. 
Our goal is for the methods to be widely accessible to the community.

Journal References
	Aicher, Christopher, Abigail Z. Jacobs, and Aaron Clauset. `Adapting 
    the Stochastic Block Model to Edge-Weighted Networks.' (2013).  
    arXiv:1305.5782 http://arxiv.org/abs/1305.5782.
	
	Aicher, Christopher, Abigail Z. Jacobs, and Aaron Clauset. “Learning 
	Latent Block Structure in Weighted Networks.” (2013)
	arXiv:1404.0431 http://arxiv.org/abs/1404.0431.


Download all the Code:
	Get the most up-to-date versions of the complete implementations. 
	At http://tuvalu.santafe.edu/~aaronc/wsbm/ 
	
	The full code package contains the core package along with all 
    additional files. The code was written for Matlab, with optional 
    MEX functions that can be installed for additional scalability.

	Type 'help <filename>' for specific information about any file in MATLAB.
	
The Core Package Overview: 
	The implementation of the Bayesian variational algorithm consists of 
    the following files:

    wsbm.m - The main function for inferring community structure.
	
	setup_distr.m - A helper function for selecting exponential family 
        distributions to use in the wsbm.m model. This function is 
        required function for wsbm.m to work properly.
	
    Adj2Edg.m, Edg2Adj.m - Helper functions for converting between 
        adjacency matrices and edge lists. These functions are required 
        for wsbm.m to work properly.
	
    /private - A folder containing MATLAB and MEX files that are used 
        privately by wsbm.m. This is a required folder for wsbm.m to 
        work properly. 
	
	InstallMEXFiles.m - A script file detailing how to install optional 
        (but highly recommended) MEX functions (in the private folder) 
        for wsbm.m to use. 
	
Additional Tools:
	There are several additional tools included in the full package, 
    which provide some demo scripts to help users learn the software, 
    some visualization scripts for displaying the posterior vertex-label 
    probabilities and matrices, and for generating synthetic networks.
	
	plotMu.m - Visualizing the Posterior Vertex-Label Probabilities	
	This function visualizes the posterior label probabilities each vertex.

	plotWSBM.m - Visualizing Adjacency Matrices
	This function visualizes the block structure of the network using the 
    labels inferred.

	generateEdges.m - Generate Synthetic Networks
	This function generates synthetic networks. For advanced users.
	
	varInfo.m - Calculates variation of information 
	nmi.m - Calculates normalized mutual information
	
	wsbmLooper.m - Wrapper Function
	This function is a wrapper for fitting and testing multiple models at 
	once. For advanced users.

	crossValidEdg.m - Data Partition Function
	This function splits an edge list into a training and test set for cross validation
	
Demo Script Files:
	The following script files are examples of how to use the package.
	
	WSBMDemo.m - Demo Script File
	This script file demonstrates examples of how to use wsbm.m	
	
	WSBMOptionsDemo.m - Demo Script File
	This script file demonstrates examples of how to change default wsbm.m options
	
	WSBMAdvancedDemo.m - Demo Script File
	This script file demonstrates examples of how to perform model selection and 
	cross valdiation with wsbm.m
	
	DegreeCorrectedDemo.m - Demo Script File
	This script file demonstrates examples of how to use degree correction in wsbm.m

	BIWSBMDemo.m - Demo Script File
	This script file demonstrates examples of how to use biwsbm.m for bipartite networks
	

A note about MATLAB compatibility
	The MATLAB functions were designed to be compatible with MATLAB v7.13
    (2011). They are not necessarily compatible with older versions of MATLAB.
	
    Also see InstallMEXFiles.m for details on how to install optional 
    (but highly recommended) MEX functions for wsbm.m to use.
	
A note about bugs
	The code provided here is provided as-is, with no warranty, with no 
    guarantees of technical support or maintenance, etc. If you experience 
    problems while using the code, please let Christopher Aicher know via 
    email (chrsitopher.aicher@colorado.edu).

Finally, if you use our code in an academic publication, it would be 
courteous of you to thank Christopher Aicher in your acknowledgements for 
providing you with implementations of the methods.

Updates
9 January 2014 - v1.0 -  First Version Posted
3 April 2014 - v1.1 - Documentation Updated and Helper Functions added
19 May 2014 - v1.2 - Additional Demos added and Bug-Fixes 

---------------------------------------------------------------------------
Copyright (C) 2013-2014 Christopher Aicher

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (LICENSE.txt).  
If not, see <http://www.gnu.org/licenses/>

