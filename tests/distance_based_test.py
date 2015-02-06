import sys
from unittest import TestCase

import ptreegen
from ptreegen.NeigborJoining import NeigborJoining


class TestNeigborJoining(TestCase):

    def setUp(self):
        self.test_matrix = [
            [0, 5, 4, 7, 6, 8]
            , [5, 0, 7, 10, 9, 11]
            , [4, 7, 0, 7, 6, 8]
            , [7, 10, 7, 0, 5, 9]
            , [6, 9, 6, 5, 0, 8]
            , [8, 11, 8, 9, 8, 0]
        ]

    def test___call__(self):
        functor = NeigborJoining(self.test_matrix, ["A", "B", "C", "D", "E", "F"])
        functor().show()

class TestComputation(TestCase):

    def setUp(self):
        self.options = {
            "algorithm" : ptreegen.TBAlgorithms.NJ
            , "seq_type" : ptreegen.SeqTypes.AA
            , "ali_path" : "test_files/keratins_ali.fasta"
            , "include_gaps" : True
            , "remove_poor" : True
            , "gap_cutoff" : 0.8
            , "pair_cutoff" : 0.3
            , "gap_penalty" : 0.5
        }

    def test___init__(self):
        computation = ptreegen.Computation(self.options)
        computation.tree.show()

        self.options["include_gaps"] = False
        computation = ptreegen.Computation(self.options)
        computation.tree.show()

def main(args):
    pass


if __name__ == "__main__":
    main(sys.argv)