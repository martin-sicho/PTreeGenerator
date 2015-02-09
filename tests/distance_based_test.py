import unittest
from unittest import TestCase

import ptreegen
from ptreegen.neigbor_joining import NeigborJoining


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
        nj = NeigborJoining(self.test_matrix, names=["A", "B", "C", "D", "E", "F"])
        nj.tree.show()

class TestComputation(TestCase):

    def setUp(self):
        self.options = {
            "method" : ptreegen.TreeBuildAlgorithms.NJ
            , "sequence_type" : ptreegen.SeqTypes.AA
            , "alignment_file" : "test_files/keratins_ali.fasta"
            , "no_gaps" : False
            , "no_cleaning" : False
            , "gap_cutoff" : 0.8
            , "pair_cutoff" : 0.3
            , "gap_penalty" : 0.5
            , "dist_measure" : ptreegen.DistMeasures.P_DISTANCE
        }

    def test___init__(self):
        computation = ptreegen.Computation(self.options)
        computation.tree.show()

        self.options["no_gaps"] = True
        computation = ptreegen.Computation(self.options)
        computation.tree.show()


if __name__ == '__main__':
    unittest.main()