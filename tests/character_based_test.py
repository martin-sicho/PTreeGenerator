import unittest

from ete2 import Tree, TreeStyle
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from ptreegen.parsimony import SmallParsimony


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
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = False
        # ts.mode = "c"
        ts.arc_start = 0
        ts.arc_span = 360
        self.circular_style = ts
        tree.set_style(ts)
        # tree.show()
        # tree.show(tree_style=ts)
        self.tree = tree
        self.alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAG", generic_dna), id="Alpha"),
                SeqRecord(Seq("AGA", generic_dna), id="Beta"),
                SeqRecord(Seq("AAA", generic_dna), id="Gamma"),
                SeqRecord(Seq("GGA", generic_dna), id="Delta"),
            ]
        )

    def test_solve(self):
        # self.tree.show()
        parsimony = SmallParsimony(self.tree, self.alignment)
        self.assertEqual(parsimony.cost, 4)


if __name__ == '__main__':
    unittest.main()
