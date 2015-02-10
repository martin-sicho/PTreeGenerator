# Introduction
A sort of primitive command line tool which lets the user input a multiple sequence alignment in FASTA format and then generates an unrooted phylogenetic tree from it. For more information on the functionality see the **[Usage](#usage)** section of this README.

This project was created as a part of a bioinformatics algorithms course at my university. It is nothing great, but does the job. If you are interested in the source code see the **[Source Code](#source-code)** section of this README or jump straight to the [documentation](http://web.vscht.cz/~sichom/PTreeGen/docs/namespaceptreegen.html).

# Installation
The program is written in the Python programming language so you will need a working Python interpreter. It should work fine with any 2.7.x series of the interpreter.

It also makes use of two bioinformatics toolkits:

 - [Biopython](http://biopython.org/wiki/Main_Page) (version 1.65)
 - [ETE Toolkit](http://etetoolkit.org/) (version 2.2)
 
Both can be easily obtained from the [Python Package Index](https://pypi.python.org/pypi). Refer to the links above on how to install these.

# Usage

There are two methodologies that can be used to build a phylogenetic tree:

 1. **distance based**: this approach uses the distances between biological sequences to determine the tree toplogy, at the moment this is handled by an implementation of the popular [Neigbor-Joining algorithm](http://en.wikipedia.org/wiki/Neighbor_joining).
 2. **character based:** uses the least evolutionary change principle (minimum number of transformations needed to transform one sequence to another) to build the tree, this approach is called the [Maximum Parsimony](http://en.wikipedia.org/wiki/Maximum_parsimony_%28phylogenetics%29).
 
There are also a few options for handling gaps and alignment cleanup such as removing badly conserved regions or all gapped regions entirely.

For a complete set of instructions and options, run the program with the ```-h``` option or read it here:

```
usage: ptreegen [-h] [--version] [-m METHOD] [-i PARS_TREE_COUNT] [-g]
                [-p GAP_PENALTY] [-c] [-u GAP_CUTOFF] [-r PAIR_CUTOFF]
                [-s SEQUENCE_TYPE] [-d DIST_MEASURE] [-f OUT_FORM]
                [-t TREE_TYPE]
                alignment_file

This is a simple tool for generating phylogenetic trees from multiple sequence
alignments. It implements two tree building approaches (Neigbor-Joining and
Maximum Parsimony).It can also do some simple visualizations and export the
tree in the Newick format.

positional arguments:
  alignment_file        Multiple sequence alignment in FASTA format.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -m METHOD, --method METHOD
                        Method used to build the tree. One of the following:
                        "NJ" for neigbor joining and "PARSIMONY" for a
                        parsimony method. Neigbor joining is the default
                        method.
  -i PARS_TREE_COUNT, --pars-tree-count PARS_TREE_COUNT
                        Number of trees used to build a consensus tree from
                        when using the parsimony method. Default value is
                        1000.
  -g, --no-gaps         Remove all gapped postions from the alignment before
                        tree building.
  -p GAP_PENALTY, --gap-penalty GAP_PENALTY
                        Gap penalty. Default value is 0.5.
  -c, --no-cleaning     Do not clean badly conserved regions from the
                        alignment before tree building.
  -u GAP_CUTOFF, --gap-cutoff GAP_CUTOFF
                        When cleaning the alignment, keep only columns with
                        non-gap frequency above this threshold. Default value
                        is 0.8.
  -r PAIR_CUTOFF, --pair-cutoff PAIR_CUTOFF
                        When cleaning the alignment, keep only columns where
                        the frequency of identical pairs is above this
                        threshold. Default value is 0.3.
  -s SEQUENCE_TYPE, --sequence-type SEQUENCE_TYPE
                        Type of the sequences in the alignment (proteins by
                        default): "AA" for proteins, "DNA" for DNA and "RNA"
                        for RNA.
  -d DIST_MEASURE, --dist-measure DIST_MEASURE
                        Distance function to be used to compute distance
                        between sequences (Jukes-Cantor by default):
                        "P_DISTANCE" for p-distance, "POISSON_CORRECTED" for
                        Poisson correction and "JUKES_CANTOR" for Jukes-
                        Cantor.
  -f OUT_FORM, --out-form OUT_FORM
                        The output formats for the resulting tree as a comma
                        separated list. Possible options: "PRINT" prints the
                        tree to command line, "NEWICK" saves the tree in
                        newick format to a file in the input directory,
                        "IMAGE_PNG" saves the tree as a PNG image in the input
                        directory, "IMAGE_SVG" saves the tree as a PNG image
                        in the input directory, "GUI" shows the tree in a
                        graphical viewer.
  -t TREE_TYPE, --tree-type TREE_TYPE
                        The type of the tree to be rendered. The default is
                        circular. Can be one of: "CIRC" for circular or "RECT"
                        for rectangular.

For more info see the GitHub page (https://github.com/martin-
sicho/PTreeGenerator) or contact 'sichom@vscht.cz'.
```

# Source Code

The source code can be divided into two parts:

 2. **the library**: the whole [ptreegen](./src/ptreegen) package. It is basically a self contained library of various functional elements that interface with each other. It facilitates all of the core functionality from tree building to visualization.
 1. **the main module**: [this](./src/main.py) can be regarded as a frontend to the [ptreegen](./src/ptreegen) package. It defines command line arguments and feeds user input to an instance of the [Computation](./src/ptreegen/computation.py) class, which is the starting point for every tree building procedure and also stores the computed results.
 
There is a complete source code documentation available [here](http://web.vscht.cz/~sichom/PTreeGen/docs/namespaceptreegen.html).

#Credits

I would like to thank [Francois-Jose Serra](https://github.com/fransua) for sparing me the pain of having to implement a consensus tree selection function for the Maximum Parsimony approach by including [it](https://github.com/fransua/utils/blob/master/pmodeltest/consensus.py) in his [repository](https://github.com/fransua/utils).
