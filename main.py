import sys
from argparse import ArgumentParser

import ptreegen


def main(args):
    description="This is a simple tool for generating phylogenetic trees from multiple sequence alignments. " \
                "It was created as a school project for a bioinformatics algorithms course."
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
    arg_parser.add_argument('-t'
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
                            , help='Type of the sequences in the alignment: "'
                                   + ptreegen.SeqTypes.AA + '" for proteins, "'
                                   + ptreegen.SeqTypes.DNA + '" for DNA and "'
                                   + ptreegen.SeqTypes.RNA + '" for RNA.'
    )

    arguments = arg_parser.parse_args()
    results = ptreegen.Computation(vars(arguments))
    results.tree.show()

if __name__ == "__main__":
    main(sys.argv)


