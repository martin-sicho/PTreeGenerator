
class SmallParsimony:

    def __init__(self, tree, alignment):
        self._tree = tree
        self._alignment = alignment
        self._seqencesDict = {x.id : x.seq for x in list(alignment)}
        leaf_names = set(tree.get_leaf_names())
        if len(self._seqencesDict) != len(leaf_names):
            raise RuntimeError("Number of sequnces in alignment:\n" + alignment
                               + "\ndoes not match the number of leaves in tree:\n" + tree)
        for seq_id in self._seqencesDict:
            if seq_id not in leaf_names:
                raise RuntimeError("Sequence ID (" + seq_id
                                   + ") does not match any leaf in tree:\n" + tree)

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

