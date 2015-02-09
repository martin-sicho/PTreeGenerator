from ete2 import Tree

from DistanceMatrix import DistanceMatrix

class NeigborJoining:

    def __init__(self, distMatrtix, alignment=None, names=None):
        self._alignment = alignment
        self._names = names
        if self._alignment:
            self._distMatrix = DistanceMatrix(distMatrtix, [x.id for x in self._alignment])
        elif self._names:
            self._distMatrix = DistanceMatrix(distMatrtix, self._names)
        else:
            raise RuntimeError("You must pass either an alignment or a list of names.")

    @property
    def tree(self):
        L = self._distMatrix.columnNames
        tree = Tree()
        tree.name = "root"
        tree.dist = 0
        tree.unroot()
        for seq in L:
            tree.add_child(name=seq, dist=0)

        iter_count = 1
        while len(L) > 2:
            nearest_nbs = self._distMatrix.getNearestNeigbors()
            node_i = tree.search_nodes(name=nearest_nbs[0])[0]
            node_j = tree.search_nodes(name=nearest_nbs[1])[0]
            L.remove(nearest_nbs[0])
            L.remove(nearest_nbs[1])

            node_k = Tree()
            node_k.dist = 0
            node_k.name = "X" + str(iter_count)
            d_ij = self._distMatrix.getDistance(node_i.name, node_j.name)
            assert d_ij > 0
            d_ik = 0.5 * d_ij + 0.5 * (self._distMatrix.getSeparation(node_i.name) - self._distMatrix.getSeparation(node_j.name))
            d_jk = 0.5 * d_ij + 0.5 * (self._distMatrix.getSeparation(node_j.name) - self._distMatrix.getSeparation(node_i.name))
            assert d_jk > 0
            assert d_ik > 0

            tree.remove_child(node_i)
            tree.remove_child(node_j)
            node_k.add_child(node_i, dist=d_ik)
            node_k.add_child(node_j, dist=d_jk)
            tree.add_child(node_k)

            d_km = []
            for node_m in L:
                d_km.append(0.5 * (self._distMatrix.getDistance(node_i.name, node_m) + self._distMatrix.getDistance(node_j.name, node_m) - d_ij) )
                assert d_km > 0

            self._distMatrix.removeData((node_i.name, node_j.name))
            self._distMatrix.appendData(d_km, node_k.name)

            iter_count+=1
            L = self._distMatrix.columnNames

        last_nodes = tree.get_children()
        d_ij = self._distMatrix.getDistance(last_nodes[0].name, last_nodes[1].name)
        last_nodes[0].add_child(last_nodes[1], dist=d_ij)
        tree = last_nodes[0]

        return tree