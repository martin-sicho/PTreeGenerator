## @package parsimony
# Contains two classes (ptreegen::parsimony::SmallParsimony
# and ptreegen::parsimony::LargeParsimony) that implement
# the basic steps of the parsimony approach to tree building.
#

from random import shuffle

from ete2 import Tree
from Bio.Align import MultipleSeqAlignment

from utilities import *

##
# Represents a solution to the small parsimony problem
# (tree is known and we are intersted
# in the parsimony score of the tree).
#
# It implements the Fitch's algorithm to score the tree,
# therefore the input tree should be a rooted binary tree,
# even though this implementation is extended to work with any
# type of tree.
#
class SmallParsimony:

    ##
    # Constructor takes tree and a corresponding alignment as input
    # and saves references to them as SmallParsimony::_tree
    # and SmallParsimony::_alignment.
    #
    # @param tree the tree to be scored
    # @param alignment alignment corresponding to the input tree
    def __init__(self, tree, alignment):
        self._tree = tree
        self._alignment = alignment
        self._seqencesDict = {x.id : x.seq for x in list(alignment)}
        leaf_names = set(tree.get_leaf_names())
        if len(self._seqencesDict) != len(leaf_names):
            raise RuntimeError("Number of sequnces in alignment:\n" + repr(alignment)
                               + "\ndoes not match the number of leaves in tree:\n" + repr(tree))
        for seq_id in self._seqencesDict:
            if seq_id not in leaf_names:
                raise RuntimeError("Sequence ID (" + seq_id
                                   + ") does not match any leaf in tree:\n" + repr(tree))

        self._treeCharacterDict = dict()
        self._cost = float('inf')

    ##
    # A getter for the cost value. Calls the
    # SmallParsimony::_solve() method to compute it.
    #
    @property
    def cost(self):
        self._cost = 0
        self._solve()
        return self._cost

    ##
    # Iteraterates over each column of the alignment
    # and computes the parsimony score for each.
    # It calls the SmallParsimony::_assign() method
    # to score each column and add the score to
    # the SmallParsimony::_cost member.
    #
    def _solve(self):
        for col_idx in range(self._alignment.get_alignment_length()):
            self._treeCharacterDict = dict()
            self._assign(self._tree, col_idx)

    ##
    # Recursive method implementing the Fitch's algorithm.
    # It traverses the tree from an arbitrary node (preferably root)
    # and computes the overall parsimony score.
    #
    # @param node the node to start traversing the tree from (preferably root)
    # @param site_idx index of the column in the alignment that is being processed
    def _assign(self, node, site_idx):
        if node.is_leaf():
            self._treeCharacterDict[node] = set(self._seqencesDict[node.name][site_idx])
        else:
            for child in node.children:
                self._assign(child, site_idx)

            character_set = set(self._treeCharacterDict[node.children[0]])
            for child in node.children:
                character_set.intersection_update(self._treeCharacterDict[child])

            if character_set:
                self._treeCharacterDict[node] = character_set
            else:
                for child in node.children:
                    character_set.update(self._treeCharacterDict[child])
                self._treeCharacterDict[node] = character_set
                self._cost += 1
    ##
    # @var _tree
    # Refernce to the tree being scored.
    # @var _alignment
    # Refernce to the coresponding alignment.
    # @var _treeCharacterDict
    # A dictionary that stores the character sets asigned to each node
    # everytime the SmallParsimony::_assign() method is executed
    # for each column of the alignment.
    # @var _cost
    # Variable to store the intermediate results.
    #

class LargeParsimony:

    def __init__(self, alignment, steps=1000):
        self._qpSteps = int(steps)
        self._alignment = alignment
        self._sequencesDict = {x.id : x.seq for x in list(self._alignment)}
        quartet_combinations = []
        for combination in combination_utils.uniqueCombinationsGenerator(self._sequencesDict.keys(), 4):
            quartet_combinations.append(combination)
        self._optimalQuartets = self.getOptimalQuartets(quartet_combinations)
        self.tree = self.quartetPuzzling(self._qpSteps)



    def quartetPuzzling(self, steps):
        seq_ids = self._sequencesDict.keys()
        trees = []
        for step in range(steps):
            shuffle(seq_ids)
            first_quartet = self._optimalQuartets[self.getQuartetID(seq_ids[0:4])]["topology"]
            rooted_tree = self.treeFromQuartet(first_quartet)
            tree = rooted_tree.children[0]
            tree.add_child(rooted_tree.children[1])
            # tree.show()

            for i in range(4,len(seq_ids)):
                tree_utils.initEdgeLengths(tree, 0)

                quartets = []
                for triplet in combination_utils.combinationsGenerator(seq_ids[0:i], 3):
                    triplet.append(seq_ids[i])
                    quartets.append(tuple(triplet))

                qt_topos_found = set()
                for quartet in quartets:
                    optimal_qt_topo_id = self._optimalQuartets[self.getQuartetID(quartet)]["topology_id"]
                    qt_topo_id = self.getTopologyID(quartet)
                    if qt_topo_id == optimal_qt_topo_id and qt_topo_id not in qt_topos_found:
                        qt_topos_found.add(qt_topo_id)
                        tree_utils.increaseCostOnPath(tree, quartet[0], quartet[1])

                # choose edge with minimum cost, delete it and add new leaf seq_ids[i]
                shortest_edge = tree_utils.findShortestEdge(tree)
                new_node = Tree(name=shortest_edge[0].name + "_" + shortest_edge[1].name)
                new_node.add_child(name=seq_ids[i])
                detached = shortest_edge[1].detach()
                shortest_edge[0].add_child(new_node)
                new_node.add_child(detached)
                # tree.show()

            tree_utils.initEdgeLengths(tree, 1)
            trees.append(tree)

        # find consensus tree
        return tree_utils.findConsensusTree(trees)

    def getOptimalQuartets(self, quartets):
        optimal_quartets = dict()
        for quartet in quartets:
            quartet_id = self.getQuartetID(quartet)
            assert quartet_id not in optimal_quartets

            trees = {tuple(quartet) : self.treeFromQuartet(quartet)}
            for i in range(0,2):
                temp = quartet[i]
                quartet[i] = quartet[2]
                quartet[2] = temp
                trees[tuple(quartet)] = self.treeFromQuartet(quartet)

            min_cost = float("inf")
            for quartet_key, tree in trees.iteritems():
                alignment = MultipleSeqAlignment([])
                for record in self._alignment:
                    if record.id in quartet:
                        alignment.append(record)
                small_parsimony = SmallParsimony(tree, alignment)
                if small_parsimony.cost < min_cost:
                    min_cost = small_parsimony.cost
                    optimal_quartets[quartet_id] = {
                        "topology" : quartet_key
                        , "topology_id" : self.getTopologyID(quartet_key)
                    }
        return optimal_quartets

    @staticmethod
    def getQuartetID(quartet):
        return "".join(sorted(quartet))

    @staticmethod
    def getTopologyID(quartet):
        return LargeParsimony.getQuartetID(quartet[0:2]) + LargeParsimony.getQuartetID(quartet[2:4])

    @staticmethod
    def treeFromQuartet(quartet):
        root = Tree()
        root.name = "root"
        left = root.add_child(name="Left")
        left.add_child(name=quartet[0])
        left.add_child(name=quartet[1])
        right = root.add_child(name="Right")
        right.add_child(name=quartet[2])
        right.add_child(name=quartet[3])
        for desc in root.iter_descendants():
            desc.dist = 0
        return root