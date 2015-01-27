import sys

from enums import *
from Computation import Computation

def main(args):
    options = {
        "algorithm" : TBAlgorithms.NJ
        , "seq_type" : SeqTypes.AA
        , "ali_path" : "testfiles/keratins_ali.fasta"
        , "include_gaps" : True
        , "remove_poor" : True
        , "gap_cutoff" : 0.8
        , "pair_cutoff" : 0.3
        , "gap_penalty" : 0.5
    }
    computation = Computation(options)
    computation.tree.show()


if __name__ == "__main__":
    main(sys.argv)


