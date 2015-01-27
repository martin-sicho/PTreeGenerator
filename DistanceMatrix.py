import numpy as np

class DistanceMatrix:

    def __init__(self, matrix, names=None):
        self.distMatrix = np.array(matrix, float)
        assert len(self.distMatrix.shape) == 2
        assert self.distMatrix.shape[0] == self.distMatrix.shape[1]

        if not names:
            generated_names = []
            for i in range(self.size):
                generated_names.append("AUTOGEN_" + str(i+1))
            self.columnNames = list(generated_names)
        else:
            self.columnNames = list(names)

        assert len(self.columnNames) == len(self.distMatrix)

    @property
    def size(self):
        return len(self.distMatrix)

    def getSeparation(self, name=None):
        dist_sum = 0
        if name:
            idx = self.columnNames.index(name)
            dist_sum = self.distMatrix[idx].sum()
        else:
            dist_sum = self.distMatrix.sum(axis=0)
        return dist_sum / (self.size - 2)

    def getNearestNeigbors(self):
        min_obj_value = None
        nearest_nbrs = tuple()
        separation = self.getSeparation()
        for i in range(self.size):
            for j in range(self.size):
                if j > i:
                    obj_value = self.distMatrix[i, j] - separation[i] - separation[j]
                    if not min_obj_value or obj_value < min_obj_value:
                        min_obj_value = obj_value
                        nearest_nbrs = (self.columnNames[i], self.columnNames[j])
        return nearest_nbrs

    def getDistance(self, name_i, name_j):
        return self.distMatrix[self.columnNames.index(name_i), self.columnNames.index(name_j)]

    def getIdx(self, name):
        return self.columnNames.index(name)

    def getName(self, idx):
        return self.columnNames[idx]

    def getNames(self):
        return list(self.columnNames)

    def removeData(self, names):
        indices = (self.columnNames.index(x) for x in names)
        for idx in indices:
            self.distMatrix = np.delete(self.distMatrix, idx, axis=0)
            self.distMatrix = np.delete(self.distMatrix, idx, axis=1)
            self.columnNames.pop(idx)

    def appendData(self, data, name):
        arr = np.array([data], float)
        self.distMatrix = np.append(self.distMatrix, arr, axis=0)
        xs = []
        for x in data:
            xs.append([x])
        xs.append([0])
        arr = np.array(xs, float)
        self.distMatrix = np.append(self.distMatrix, arr, axis=1)
        self.columnNames.append(name)


