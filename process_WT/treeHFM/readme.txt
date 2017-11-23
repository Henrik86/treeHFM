The treeHFM (hidden Factor graph model) generalises Hidden Markov Models to tree structured data. The distinctive feature of treeHFM is that it learns a transition matrix for first order (sequential) and for second order (splitting) events. It can be applied to all discrete and continuous data that is structured as a binary tree. In the case of continuous observations, treeHFM has Gaussian distributions as emissions.

In order to run treeHFM, you have to compile the mex files first. 

Therefore, move into the treeHFM folder, set a C++ compiler (using ’mex -setup’) and run compile_mex. 

Copyright (c) 2016 Henrik Failmezger