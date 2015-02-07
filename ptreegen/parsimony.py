import copy

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

        self._optimalQuartets = dict()
        self.saveOptimalQuartets(quartet_combinations)
        print quartet_combinations
        print self._optimalQuartets.keys()

    def saveOptimalQuartets(self, quartets):
        for quartet in quartets:
            quartet_id = "".join(quartet)
            quartet_copy = copy.deepcopy(quartet)
            assert quartet_id not in self._optimalQuartets.keys()
            trees = [self.treeFromQuartet(quartet_copy)]
            for i in range(0,2):
                temp = quartet_copy[i]
                quartet_copy[i] = quartet_copy[2]
                quartet_copy[2] = temp
                trees.append(self.treeFromQuartet(quartet_copy))
            min_cost = float("inf")
            for tree in trees:
                alignment = MultipleSeqAlignment([])
                for record in self._alignment:
                    if record.id in quartet_copy:
                        alignment.append(record)
                small_parsimony = SmallParsimony(tree, alignment)
                if small_parsimony.cost < min_cost:
                    min_cost = small_parsimony.cost
                    self._optimalQuartets[quartet_id] = tree # TODO: try better representation

    def treeFromQuartet(self, quartet):
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