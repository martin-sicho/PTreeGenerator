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
# \sa LargeParsimony
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
        self._sequencesDict = {x.id : x.seq for x in list(alignment)}
        leaf_names = set(tree.get_leaf_names())
        if len(self._sequencesDict) != len(leaf_names):
            raise RuntimeError("Number of sequnces in alignment:\n" + repr(alignment)
                               + "\ndoes not match the number of leaves in tree:\n" + repr(tree))
        for seq_id in self._sequencesDict:
            if seq_id not in leaf_names:
                raise RuntimeError("Sequence ID (" + seq_id
                                   + ") does not match any leaf in tree:\n" + repr(tree))

        self._treeCharacterDict = dict()
        self._cost = float('inf')

    ##
    # A getter for the cost value. Calls the
    # SmallParsimony::_solve() method to compute it.
    #
    # @return parsimony score as a single value
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
            self._treeCharacterDict[node] = set(self._sequencesDict[node.name][site_idx])
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
    # @var _sequencesDict
    # Associates sequence id string with a sequence instance.
    # @var _treeCharacterDict
    # A dictionary that stores the character sets asigned to each node
    # everytime the SmallParsimony::_assign() method is executed
    # for each column of the alignment.
    # @var _cost
    # Variable to store the intermediate results.
    #

##
# Represents a solution to the large parsimony problem
# (the tree toplogy is unknown as oposed to the
# small parsimony problem).
#
# It implements the quartet puzzling heuristic to find
# a set of possible tree topologies and it then builds
# a consensus tree employing the 50% majority rule
# (implemented in ptreegen::utilities::tree_utils::findConsensusTree()).
#
# For more information on the quartet puzzling algorithm see: <br />
# Strimmer, Korbinian, and Arndt Von Haeseler.
# <em>Quartet puzzling: a quartet maximum-likelihood method
# for reconstructing tree topologies.</em>
# Molecular Biology and Evolution 13.7 (1996): 964-969.<br />
# or the appropriate chapter in: <br />
# Bockenhauer, H.; Bongartz, D. <em>Algorithmic Aspects of
# Bioinformatics</em>, 1st ed.; Springer: Berlin, 1998.
#
# \sa SmallParsimony
# \sa LargeParsimony::quartetPuzzling()
#
class LargeParsimony:

    ##
    # Initializes the members used during the computation.
    # It also generates the optimal quartets by first
    # computing the parsimony cost of all possible unique quartet
    # combinations and then selecting the quartet with minimum
    # parsimony cost for every four taxa.
    #
    # @param alignment multiple sequence alignment
    # @param steps number of possible tree topologies that will be generated
    def __init__(self, alignment, steps=1000):
        self._qpSteps = int(steps)
        self._alignment = alignment
        self._sequencesDict = {x.id : x.seq for x in list(self._alignment)}
        quartet_combinations = []
        for combination in combination_utils.uniqueCombinationsGenerator(self._sequencesDict.keys(), 4):
            quartet_combinations.append(combination)
        self._optimalQuartets = self.getOptimalQuartets(quartet_combinations)
        self._tree = None

    ##
    # A getter for the obtained consensus tree.
    #
    # @return reference to the tree instance
    @property
    def tree(self):
        if self._tree:
            return self._tree
        else:
            self._tree = self.quartetPuzzling(self._qpSteps)
            self._tree.get_tree_root().unroot()
        return self._tree

    ##
    # A getter for the parsimony score of the
    # computed tree.
    #
    # @return parsimony cost of the tree as a single value
    @property
    def cost(self):
        return SmallParsimony(self.tree, self._alignment).cost

    ##
    # Iteratively applies the quartet puzzling
    # heuristic to obtain intermediate tree topologies.
    #
    # Implemented as described in: <br />
    # Bockenhauer, H.; Bongartz, D. <em>Algorithmic Aspects of
    # Bioinformatics</em>, 1st ed.; Springer: Berlin, 1998.
    #
    # @param steps number of puzzling steps (equals the number of intermediate
    # tree topologies)
    # @return consensus tree according to the 50% majority rule
    def quartetPuzzling(self, steps):
        seq_ids = self._sequencesDict.keys()
        if len(seq_ids) < 4:
            tree = Tree()
            for seq_id in seq_ids:
                tree.add_child(name=seq_id)
            return tree

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
                        self.increaseCostOnPath(tree, quartet[0], quartet[1])

                # choose edge with minimum cost, delete it and add new leaf seq_ids[i]
                shortest_edge = tree_utils.findShortestEdge(tree)
                # new_node = Tree(name=shortest_edge[0].name + "_" + shortest_edge[1].name)
                new_node = Tree()
                new_node.add_child(name=seq_ids[i])
                detached = shortest_edge[1].detach()
                shortest_edge[0].add_child(new_node)
                new_node.add_child(detached)
                # tree.show()

            tree_utils.initEdgeLengths(tree, 1)
            trees.append(tree)

        # find consensus tree
        return tree_utils.findConsensusTree(trees)

    ##
    # Finds the optimal topologies for every four taxa.
    #
    # @param quartets all possible unique combinations of taxa of size four
    # @return a dictionary that associates every for taxa with the optimal topology
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

    ##
    # Provides a unique identification
    # of the quartet by sorting the
    # names of the taxa in the quartet
    # alphabetically and then concatenating
    # their identification strings, hence
    # allowing quartets to be identified with their
    # respective taxa.
    #
    # Note that the taxa identification
    # strings should be unique in order
    # to eliminate collisions.
    #
    # @param quartet quartet to be identified
    @staticmethod
    def getQuartetID(quartet):
        return "".join(sorted(quartet))

    ##
    # Provides a string that uniquely
    # identifies a quartet topology.
    #
    # Quartet topology does not depend
    # on the order of the first two and
    # the last two elements in the quartet.
    # This method is needed to distinguish
    # two quartets of the same four taxa
    # that have different topologies.
    #
    # It also allows to have the first two
    # and last two elements of the quartet
    # in an arbitrary order.
    #
    # @param quartet quartet to be identified
    @staticmethod
    def getTopologyID(quartet):
        return LargeParsimony.getQuartetID(quartet[0:2]) + LargeParsimony.getQuartetID(quartet[2:4])

    ##
    # Constructs a tree from the provided quartet.
    # The tree is rooted so that it can be used
    # in the Fitch's algorithm for solving the small
    # parsimony problem.
    #
    # @param quartet quartet to be transformed
    # @return a reference to the root of the tree
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

    ##
    # Finds a path between two vertices
    # and then increases the weight of all edges
    # on that path by one.
    #
    # @param tree reference to any vertex of the processed tree
    # @param start reference to the start node
    # @param dest reference to the end node
    @staticmethod
    def increaseCostOnPath(tree, start, dest):
        start_node = tree.search_nodes(name=start)[0]
        end_node = tree.search_nodes(name=dest)[0]
        node_from = dict()
        to_visit = set()
        visited = set()

        to_visit.add(start_node)
        node_from[start_node] = None
        path_found = False
        while to_visit and not path_found:
            current = to_visit.pop()
            neighbors = set(current.children)
            if current.up:
                neighbors.add(current.up)
            for neighbor in neighbors:
                if neighbor not in visited:
                    node_from[neighbor] = current
                    if neighbor == end_node:
                        path_found = True
                        break
                    else:
                        to_visit.add(neighbor)
            visited.add(current)

        node = end_node
        while node:
            if node.name != "Left":
                node.dist += 1
            node = node_from[node]

    ##
    # @var _qpSteps
    # The number of quartet puzzling steps to take.
    # @var _alignment
    # Reference to the multiple sequence alignment instance
    # @var _sequencesDict
    # Associates sequence id string with a sequence instance.
    # @var _optimalQuartets
    # Associates every four taxa with their optimal quartet.
    # @var _tree
    # Reference to the tree that was computed in LargeParsimony::quartetPuzzling().
