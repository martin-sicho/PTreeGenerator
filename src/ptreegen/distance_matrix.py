## @package distance_matrix
# Constains just the ptreegen::distance_matrix::DistanceMatrix class.
#

import numpy as np

##
# Basically a wrapper around a numpy array obejct
# representing the alignment distance matrix.
#
# Performs some other additional operations usefull for tree building.
#
class DistanceMatrix:

    ##
    # Takes any "matrix-like" object and tries to convert it to a numpy array.
    #
    # @param matrix a "matrix-like" object
    # @param names optional parameter with column and row names (the taxa names)
    #
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

    ##
    # A getter for the matrix size (number of columns/taxa).
    #
    @property
    def size(self):
        return len(self._distMatrix)

    ##
    # A getter for a copy of the whole distance matrix.
    #
    @property
    def distMatrix(self):
        return np.array(self._distMatrix, float)

    ##
    # A getter for a list of column/taxa names.
    #
    @property
    def columnNames(self):
        return list(self._columnNames)

    ##
    # Returns a separation of value used in the Neigbor-joining algorithm.
    #
    # It can be computed for one sequence only (parameter name)
    # or for all sequences (no parameter).
    #
    # The separation value is computed as follows:
    # <em>sum(d_ik) / (L - 2)</em>, where <em>sum(d_ik)</em> is the
    # sum of distances from one sequence to
    # all the other sequences and L is the total number of sequences.
    #
    # @param name identification of one sequence
    # @return returns separation values for all sequences as a list
    # or one value for one sequence with the specified name
    def getSeparation(self, name=None):
        dist_sum = None
        if name:
            idx = self._columnNames.index(name)
            dist_sum = self._distMatrix[idx].sum()
        else:
            dist_sum = self._distMatrix.sum(axis=0)
        return dist_sum / (self.size - 2)

    ##
    # Finds the pair of nearest sequences.
    #
    # Finds the pair of closest sequences according to the rule
    # from the Neigbor-Joining algorithm.
    #
    # @return tuple of size two that contains the names of two nearest sequences
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

    ##
    # Returns the distance from one sequence to another.
    #
    # Based on the value from the distance matrix.
    #
    # @return one single number
    def getDistance(self, name_i, name_j):
        return self._distMatrix[self._columnNames.index(name_i), self._columnNames.index(name_j)]

    ##
    # Finds the position of a sequence in the distance matrix.
    #
    # @param name the identification of the sequence
    # @return index of a column in the matrix as a single number
    def getIdx(self, name):
        return self._columnNames.index(name)

    ##
    # Finds the name of a sequence based on its position in the matrix.
    #
    # @param idx the position in the matrix
    # @return index the identification of the sequence as string
    def getName(self, idx):
        return self._columnNames[idx]

    ##
    # Removes rows and columns for the specified sequences.
    #
    # @param names the identifications of the sequences as an iterable
    def removeData(self, names):
        indices = (self._columnNames.index(x) for x in names)
        for idx in indices:
            self._distMatrix = np.delete(self._distMatrix, idx, axis=0)
            self._distMatrix = np.delete(self._distMatrix, idx, axis=1)
            self._columnNames.pop(idx)

    ##
    # Adds a row and a column for the specified sequence.
    #
    # @param name the identification of the sequence
    # @param data data to be appended as an iterable
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

    ##
    # @var _distMatrix
    # Distance matrix as a numpy array object.
    # @var _columnNames
    # List of column names (the identification strings of the sequences).
    #
