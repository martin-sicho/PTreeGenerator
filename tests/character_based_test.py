import unittest

from ete2 import Tree, TreeStyle
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from ptreegen.parsimony import SmallParsimony, LargeParsimony


class TestSmallParsimony(unittest.TestCase):

    def setUp(self):
        tree = Tree()
        root = tree.get_tree_root()
        root.dist = 0
        root.name = "root"
        node = root.add_child(name="Left")
        node.add_child(name="Alpha")
        node.add_child(name="Beta")
        node = root.add_child(name="Right")
        node.add_child(name="Gamma")
        node.add_child(name="Delta")
        for desc in tree.iter_descendants():
            desc.dist = 0

        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = False
        ts.mode = "c"
        ts.arc_start = 0
        ts.arc_span = 360

        self.circular_style = ts
        self.exampleTree = tree
        self.alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAG", generic_dna), id="Alpha"),
                SeqRecord(Seq("AGA", generic_dna), id="Beta"),
                SeqRecord(Seq("AAA", generic_dna), id="Gamma"),
                SeqRecord(Seq("GGA", generic_dna), id="Delta"),
            ]
        )

    def test_solve(self):
        parsimony = SmallParsimony(self.exampleTree, self.alignment)
        self.assertEqual(parsimony.cost, 4)
        # self.exampleTree.show(tree_style=self.circular_style)

class testLargeParsimony(unittest.TestCase):

    def setUp(self):
        self.alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAG", generic_dna), id="Alpha"),
                SeqRecord(Seq("AGA", generic_dna), id="Beta"),
                SeqRecord(Seq("AAA", generic_dna), id="Gamma"),
                SeqRecord(Seq("GGA", generic_dna), id="Delta"),
                SeqRecord(Seq("TTA", generic_dna), id="Epsilon"),
            ]
        )

    def test_solve(self):
        parsimony = LargeParsimony(self.alignment)


if __name__ == '__main__':
    unittest.main()
