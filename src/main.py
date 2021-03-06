## @package main
# The main module with the main method.
# Handles the comand line options
# and delegates actions to other modules.

from argparse import ArgumentParser

import ptreegen

##
# The main method. First to be executed
# when the program is run.
def main():
    description="This is a simple tool for generating phylogenetic trees from multiple sequence alignments. " \
                "It implements two tree building approaches (Neigbor-Joining and Maximum Parsimony)." \
                "It can also do some simple visualizations and export the tree in the Newick format."
    epilog = "For more info see the GitHub page (https://github.com/martin-sicho/PTreeGenerator) " \
             "or contact 'sichom@vscht.cz'."
    arg_parser = ArgumentParser(prog="ptreegen", description=description, epilog=epilog)
    arg_parser.add_argument('--version'
                            , action='version'
                            , version='1.0-alpha')
    arg_parser.add_argument('alignment_file'
                            , help='Multiple sequence alignment in FASTA format.')
    arg_parser.add_argument('-m'
                            ,'--method'
                            , action='store'
                            , default=ptreegen.TreeBuildAlgorithms.NJ
                            , help='Method used to build the tree. One of the following: "'
                            + ptreegen.TreeBuildAlgorithms.NJ + '" for neigbor joining and "'
                            + ptreegen.TreeBuildAlgorithms.PARSIMONY + '" for a parsimony method. '
                            + 'Neigbor joining is the default method.'
                            )
    arg_parser.add_argument('-i'
                            , '--pars-tree-count'
                            , action='store'
                            , default=1000
                            , help="Number of trees used to build a consensus tree from when using "
                                   + "the parsimony method. Default value is 1000.")
    arg_parser.add_argument('-g'
                            , '--no-gaps'
                            , action='store_true'
                            , help="Remove all gapped postions from the alignment before tree building.")
    arg_parser.add_argument('-p'
                            , '--gap-penalty'
                            , action='store'
                            , default=0.5
                            , help="Gap penalty. Default value is 0.5.")
    arg_parser.add_argument('-c'
                            , '--no-cleaning'
                            , action='store_true'
                            , help="Do not clean badly conserved regions from the alignment before tree building.")
    arg_parser.add_argument('-u'
                            , '--gap-cutoff'
                            , action='store'
                            , default=0.8
                            , help="When cleaning the alignment, keep only columns with "
                                   + "non-gap frequency above this threshold. Default value is 0.8.")
    arg_parser.add_argument('-r'
                            , '--pair-cutoff'
                            , action='store'
                            , default=0.3
                            , help="When cleaning the alignment, keep only columns where the frequency of identical pairs is above this threshold. Default value is 0.3.")
    arg_parser.add_argument('-s'
                            , '--sequence-type'
                            , action='store'
                            , default=ptreegen.SeqTypes.AA
                            , help='Type of the sequences in the alignment (proteins by default): "'
                                   + ptreegen.SeqTypes.AA + '" for proteins, "'
                                   + ptreegen.SeqTypes.DNA + '" for DNA and "'
                                   + ptreegen.SeqTypes.RNA + '" for RNA.'
    )
    arg_parser.add_argument('-d'
                            , '--dist-measure'
                            , action='store'
                            , default=ptreegen.DistMeasures.JUKES_CANTOR
                            , help='Distance function to be used to compute distance '
                                   + 'between sequences (Jukes-Cantor by default): "'
                                   + ptreegen.DistMeasures.P_DISTANCE + '" for p-distance, "'
                                   + ptreegen.DistMeasures.POISSON_CORRECTED + '" for Poisson correction and "'
                                   + ptreegen.DistMeasures.JUKES_CANTOR + '" for Jukes-Cantor.'
    )

    # visualization arguments
    arg_parser.add_argument('-f'
                            , '--out-form'
                            , action='store'
                            , default=ptreegen.OutputForm.PRINT + "," + ptreegen.OutputForm.NEWICK
                            , help='The output formats for the resulting tree as a comma separated list. '
                                   + 'Possible options: "'
                                   + ptreegen.OutputForm.PRINT + '" prints the tree to command line, "'
                                   + ptreegen.OutputForm.NEWICK + '" saves the tree in newick format to a file in the input directory, "'
                                   + ptreegen.OutputForm.IMAGE_PNG + '" saves the tree as a PNG image in the input directory, "'
                                   + ptreegen.OutputForm.IMAGE_SVG + '" saves the tree as a PNG image in the input directory, "'
                                   + ptreegen.OutputForm.GUI + '" shows the tree in a graphical viewer.'
    )
    arg_parser.add_argument('-t'
                            , '--tree-type'
                            , action='store'
                            , default=ptreegen.TreeType.CIRC
                            , help='The type of the tree to be rendered. '
                                   + 'The default is circular. Can be one of: "'
                                   + ptreegen.TreeType.CIRC + '" for circular or "'
                                   + ptreegen.TreeType.RECT + '" for rectangular.'
    )


    arguments = arg_parser.parse_args()
    results = ptreegen.Computation(vars(arguments))
    results.showResults()

    return 0

if __name__ == "__main__":
    exit(main())


