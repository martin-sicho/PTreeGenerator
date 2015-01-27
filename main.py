import sys
from enums import *
from Bio import AlignIO
from Bio.Alphabet import generic_protein, generic_dna

class Computation:

    def __init__(self, options=None):
        self.algorithm = options["algorithm"]
        self.includeGaps = options["include_gaps"]
        self.removePoor = options["remove_poor"]


        self.seqType = None
        self.seqType = options["seq_type"]
        if self.seqType == SeqTypes.AA:
            self.aligmnment = AlignIO.read(options["ali_path"], "fasta", alphabet=generic_protein)
        if self.seqType == SeqTypes.DNA:
            self.aligmnment = AlignIO.read(options["ali_path"], "fasta", alphabet=generic_dna)
        assert self.seqType

        self.parseAlignment(self.aligmnment)

    def parseAlignment(self, alignment):
        pass

def main(args):
    options = {
        "algorithm" : TBAlgorithms.NJ
        , "seq_type" : SeqTypes.AA
        , "ali_path" : "testfiles/keratins_ali.fasta"
        , "include_gaps" : True
        , "remove_poor" : True

    }
    computation = Computation(options)


if __name__ == "__main__":
    main(sys.argv)


