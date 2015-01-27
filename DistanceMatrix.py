import numpy as np

class DistanceMatrix:

    def __init__(self, matrix, names=None):
        self._distMatrix = np.array(matrix, float)
        assert len(self._distMatrix.shape) == 2
        assert self._distMatrix.shape[0] == self._distMatrix.shape[1]

        if not names:
            generated_names = []
            for i in range(self.size):
                generated_names.append("AUTOGEN_" + str(i+1))
            self._columnNames = list(generated_names)
        else:
            self._columnNames = list(names)

        assert len(self._columnNames) == len(self._distMatrix)

    @property
    def size(self):
        return len(self._distMatrix)

    @property
    def distMatrix(self):
        return np.array(self._distMatrix, float)

    @property
    def columnNames(self):
        return list(self._columnNames)

    def getSeparation(self, name=None):
        dist_sum = 0
        if name:
            idx = self._columnNames.index(name)
            dist_sum = self._distMatrix[idx].sum()
        else:
            dist_sum = self._distMatrix.sum(axis=0)
        return dist_sum / (self.size - 2)

    def getNearestNeigbors(self):
        min_obj_value = None
        nearest_nbrs = tuple()
        separation = self.getSeparation()
        for i in range(self.size):
            for j in range(self.size):
                if j > i:
                    obj_value = self._distMatrix[i, j] - separation[i] - separation[j]
                    if not min_obj_value or obj_value < min_obj_value:
                        min_obj_value = obj_value
                        nearest_nbrs = (self._columnNames[i], self._columnNames[j])
        return nearest_nbrs

    def getDistance(self, name_i, name_j):
        return self._distMatrix[self._columnNames.index(name_i), self._columnNames.index(name_j)]

    def getIdx(self, name):
        return self._columnNames.index(name)

    def getName(self, idx):
        return self._columnNames[idx]

    # def getNames(self):
    #     return list(self._columnNames)

    def removeData(self, names):
        indices = (self._columnNames.index(x) for x in names)
        for idx in indices:
            self._distMatrix = np.delete(self._distMatrix, idx, axis=0)
            self._distMatrix = np.delete(self._distMatrix, idx, axis=1)
            self._columnNames.pop(idx)

    def appendData(self, data, name):
        arr = np.array([data], float)
        self._distMatrix = np.append(self._distMatrix, arr, axis=0)
        xs = []
        for x in data:
            xs.append([x])
        xs.append([0])
        arr = np.array(xs, float)
        self._distMatrix = np.append(self._distMatrix, arr, axis=1)
        self._columnNames.append(name)


