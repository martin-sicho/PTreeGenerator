from random import shuffle

from ete2 import Tree
from Bio.Align import MultipleSeqAlignment

class SmallParsimony:

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
        self._cost = 0

        self.solve()

    @property
    def cost(self):
        return self._cost

    def solve(self):
        for col_idx in range(self._alignment.get_alignment_length()):
            self._treeCharacterDict = dict()
            self.assign(self._tree, col_idx)

    def assign(self, node, site_idx):
        if node.is_leaf():
            self._treeCharacterDict[node] = set(self._seqencesDict[node.name][site_idx])
        else:
            for child in node.children:
                self.assign(child, site_idx)

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

class LargeParsimony:

    def __init__(self, alignment):
        self._alignment = alignment
        self._sequencesDict = {x.id : x.seq for x in list(self._alignment)}
        quartet_combinations = []
        for combination in LargeParsimony.uniqueCombinationsGenerator(self._sequencesDict.keys(), 4):
            quartet_combinations.append(combination)

        self._optimalQuartets = self.getOptimalQuartets(quartet_combinations)
        self.tree = self.quartetPuzzling() # TODO: repeat this a few times and choose consensus tree

    def quartetPuzzling(self):
        seq_ids = self._sequencesDict.keys()
        shuffle(seq_ids)
        first_quartet = self._optimalQuartets[self.getQuartetID(seq_ids[0:4])]["topology"]
        rooted_tree = self.treeFromQuartet(first_quartet)
        tree = rooted_tree.children[0]
        tree.add_child(rooted_tree.children[1])
        # tree.show()

        for i in range(4,len(seq_ids)):
            self.initEdgeLengths(tree, 0)

            quartets = []
            for triplet in self.combinationsGenerator(seq_ids[0:i], 3):
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
            shortest_edge = self.findShortestEdge(tree)
            new_node = Tree(name=shortest_edge[0].name + "_" + shortest_edge[1].name)
            new_node.add_child(name=seq_ids[i])
            detached = shortest_edge[1].detach()
            shortest_edge[0].add_child(new_node)
            new_node.add_child(detached)
            # tree.show()

        self.initEdgeLengths(tree, 1)
        return tree

    @staticmethod
    def findShortestEdge(tree):
        shortest_edge = [None, None]
        min_dist = float("inf")
        for node in tree.iter_descendants():
            if node.dist < min_dist:
                shortest_edge[0] = node.up
                shortest_edge[1] = node
        return tuple(shortest_edge)

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

        # tree.show()

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
    def initEdgeLengths(tree, value):
        tree.dist = value
        for desc in tree.iter_descendants():
            desc.dist = value

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

    @staticmethod
    def uniqueCombinationsGenerator(items, n):
            if not n:
                yield []
            else:
                for i in xrange(len(items)):
                    for comb in LargeParsimony.uniqueCombinationsGenerator(items[i+1:], n-1):
                        yield [items[i]] + comb
    @staticmethod
    def combinationsGenerator(items, n):
        if not n:
            yield []
        else:
            for i in xrange(len(items)):
                for comb in LargeParsimony.combinationsGenerator(items[:i]+items[i+1:],n-1):
                    yield [items[i]] + comb