import sys
from unittest import TestCase
from NeigborJoining import NeigborJoining

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

def main(args):
    pass


if __name__ == "__main__":
    main(sys.argv)