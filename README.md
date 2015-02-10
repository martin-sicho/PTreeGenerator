# Introduction
A sort of primitive command line tool which lets the user input a multiple sequence alignment in FASTA format and then generates an unrooted phylogenetic tree from it. For more information on the functionality see the **Usage** section of this README.

This project was created as a part of a bioinformatics algorithms course at my university. It is nothing great, but does the job. If you are interested in the source code see the **Source Code** section of this README or jump straight to the [documentation](http://web.vscht.cz/~sichom/PTreeGen/docs).

# Installation
The program is written in the Python programming language so you will need a working Python interpreter. It should work fine with any 2.7.x series of the interpreter.

It also makes use of two bioinformatics toolkits:

 - [Biopython](http://biopython.org/wiki/Main_Page) (version 1.65)
 - [ETE Toolkit](http://etetoolkit.org/) (version 2.2)
 
Both can be easily obtained from the [Python Package Index](https://pypi.python.org/pypi). Refer to the links above on how to install these.

# Usage

There are two methodologies that can be used to build a phylogenetic tree:

 1. **distance based**: this approach uses the distances between biological sequences to determine the tree toplogy, at the moment this is handled by an implementation of the popular [Neigbor-Joining algorithm](http://en.wikipedia.org/wiki/Neighbor_joining).
 2. **character based:** uses the number of mutations needed to transform one sequence to another to build the tree, this approach is called the [Maximum Parsimony](http://en.wikipedia.org/wiki/Maximum_parsimony_%28phylogenetics%29).
 
There are also a few options for handling gaps and alignment cleanup such as removing badly conserved regions or all gapped regions entirely.

For a complete set of instructions run the program with the ```-h``` option.

# Source Code

The source code can be divided into two parts:

 2. **the library**: the whole [ptreegen](./src/ptreegen) package. It is basically a self contained library of various functional elements that interface with each other. It facilitates all of the core functionality from tree building to visualization.
 1. **the main module**: [this](./src/main.py) can be regarded as a frontend to the [ptreegen](./src/ptreegen) package. It defines command line arguments and feeds user input to an instance of the [Computation](./src/ptreegen/computation.py) class, which is the starting point for every tree building procedure and also stores the computed results.
 
There is a complete source code documentation available [here](http://web.vscht.cz/~sichom/PTreeGen/docs).

#Credits

I would like to thank [Francois-Jose Serra](https://github.com/fransua) for sparing me the pain of having to implement a consensus tree selection function for the Maximum Parsimony approach by including [it](https://github.com/fransua/utils/blob/master/pmodeltest/consensus.py) in his [repository](https://github.com/fransua/utils).
