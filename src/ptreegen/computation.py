## @package computation
# Contains just the ptreegen::computation::Computation class.
#

from Bio import AlignIO
from Bio.Alphabet import generic_protein, generic_dna, generic_rna
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from ptreegen.enums import *
import distance_functions as dfuncs
from neighbor_joining import NeighborJoining
from ptreegen.parsimony import LargeParsimony




##
# Parses user specified options and delegetes
# appropriate actions to other modules. Also serves as
# a data storage of computed results.
class Computation:

    ##
    # Constructor initializes the object variables
    # and calls other modules to do so.
    #
    # @param options computation options in the form of a dictionary-like object
    def __init__(self, options):
        self.algorithm = None
        self.gapPenalty = None
        self.includeGaps = None
        self.removePoor = None
        self.gapCutoff = None
        self.pairCutoff = None
        self.seqType = None
        self.distFunction = None
        self.options = options
        self.parseOptions(self.options)
        self.alignment = self.cleanAlignment(self.alignment)
        self.distanceMatrix = None
        self.tree = self.computeTree()

    ##
    # Parses the options passed to the constructor
    #
    # @param options computation options in the form of a dictionary-like object
    def parseOptions(self, options):
        self.algorithm = options["method"]
        if self.algorithm not in (TreeBuildAlgorithms.NJ, TreeBuildAlgorithms.PARSIMONY):
            raise RuntimeError("Unknown method: " + self.algorithm)
        self.gapPenalty = options["gap_penalty"]
        if self.gapPenalty < 0 or self.gapPenalty > 1:
            raise RuntimeError("Bad gap penalty value. Must be between 0 and 1. Got: " + self.gapPenalty)
        self.includeGaps = not options["no_gaps"]
        self.removePoor = not options["no_cleaning"]
        self.gapCutoff = options["gap_cutoff"]
        if self.gapCutoff < 0 or self.gapCutoff > 1:
            raise RuntimeError("Bad gap cutoff value. Must be between 0 and 1. Got: " + self.gapCutoff)
        self.pairCutoff = options["pair_cutoff"]
        self.seqType = options["sequence_type"]
        if self.seqType == SeqTypes.AA:
            self.alignment = AlignIO.read(options["alignment_file"], "fasta", alphabet=generic_protein)
        elif self.seqType == SeqTypes.DNA:
            self.alignment = AlignIO.read(options["alignment_file"], "fasta", alphabet=generic_dna)
        elif self.seqType == SeqTypes.RNA:
            self.alignment = AlignIO.read(options["alignment_file"], "fasta", alphabet=generic_rna)
        else:
            raise RuntimeError("Unknown sequence type: " + self.seqType)
        if options["dist_measure"] == DistMeasures.P_DISTANCE:
            self.distFunction = dfuncs.p_distance
        elif options["dist_measure"] == DistMeasures.POISSON_CORRECTED:
            self.distFunction = dfuncs.poisson_corrected
        elif options["dist_measure"] == DistMeasures.JUKES_CANTOR:
            self.distFunction = dfuncs.jukes_cantor
        else:
            raise RuntimeError("Unknown distance measure: " + options["dist_measure"])
    ##
    # Method that delegates tree computation
    # to the appropriate module.
    #
    # @return reference to the computed tree object
    def computeTree(self):
        if self.algorithm == TreeBuildAlgorithms.NJ:
            self.distanceMatrix = self.computeDistanceMatrix(self.alignment, self.distFunction)
            return NeighborJoining(self.distanceMatrix, self.alignment).tree
        elif self.algorithm == TreeBuildAlgorithms.PARSIMONY:
            return LargeParsimony(self.alignment).tree
        else:
            raise RuntimeError(self.algorithm + " not implemented.")

    ##
    # Method that can be used to update the results,
    # if the Computation class changes.
    #
    def update(self):
        self.alignment = self.cleanAlignment(self.alignment)
        self.distanceMatrix = self.computeDistanceMatrix(self.alignment, self.distFunction)
        self.tree = self.computeTree()

    ##
    # Computes the distance matrix from the alignment.
    #
    # @param alignment the mutliple sequence alignment instance
    # @param distFunction the distance measure used, can be one
    # of the functions in ptreegen::distance_functions.
    # @return distance matrix as a tuple
    def computeDistanceMatrix(self, alignment, distFunction):
        dist_matrix = []
        for i,record_i in enumerate(alignment):
            seq_i = record_i.seq
            distances = []
            for j,record_j in enumerate(alignment):
                if j > i:
                    seq_j = record_j.seq
                    distances.append(distFunction(seq_i, seq_j, **self.options))
                elif i == j:
                    distances.append(0)
                else:
                    distances.append(dist_matrix[j][i])
            dist_matrix.append(tuple(distances))
        return tuple(dist_matrix)

    ##
    # Method responsible for the correct alignment cleaning.
    #
    # It removes badly conserved regions and/or regions with too many gaps,
    # if requested by the user.
    #
    # @param alignment the mutliple sequence alignment instance
    def cleanAlignment(self, alignment):
        to_remove = set()
        for col_idx in range(alignment.get_alignment_length()):
            column = alignment[:,col_idx]

            # removal of columns with gaps
            if not self.includeGaps:
                if column.find("-") != -1:
                    to_remove.add(col_idx)

            # poorly conserved regions removal
            if self.removePoor and (col_idx not in to_remove):
                # remove column if it contains too many gaps
                nongap_ratio = float(len(column) - column.count("-")) / len(column)
                if nongap_ratio < self.gapCutoff:
                    to_remove.add(col_idx)
                    continue

                # remove column if there are not enough identical residues
                all_pairs = (len(column) - column.count("-")) * ((len(column) - column.count("-")) - 1) / 2.0
                identical_pairs = 0
                for res in set(column):
                    if res != "-":
                        res_count = column.count(res)
                        if res_count > 1:
                            identical_pairs += res_count * (res_count - 1) / 2.0
                identical_pairs_ratio = identical_pairs / all_pairs
                if identical_pairs_ratio < self.pairCutoff:
                    to_remove.add(col_idx)
                    continue

        cleaned_alignment = MultipleSeqAlignment([])
        for record in alignment:
            seq = record.seq.tomutable()
            counter = 0
            for idx in to_remove:
                seq.pop(idx - counter)
                counter+=1
            cleaned_alignment.append(SeqRecord(seq.toseq(), id=record.id))

        return cleaned_alignment

    ##
    # @var algorithm
    # The methodology used to build the tree as one of those in ptreegen::enums.
    # @var gapPenalty
    # Cost of gaps when they are included in the distance computation.
    # @var includeGaps
    # Specifies if columns with gaps should be left in the alignment or deleted.
    # @var removePoor
    # Specifies if poorly conserved regions should be removed from the alignment.
    # @var gapCutoff
    # If at least x% in column are not gaps, then the column is left in the alignment.
    # @var pairCutoff
    # If at least x% of pairs in column are identical, then the column is left in the alignment.
    # @var seqType
    # The type of the input sequence as one of ptreegen::enums.
    # @var distFunction
    # Pointer to function used to compute distances between a two sequences.
    # @var options
    # Reference to the dictionary like object.
    # @var alignment
    # Reference to an object representing the multiple sequence alignment.
    # @var distanceMatrix
    # The distance matrix for the alignment.
    # @var tree
    # The generated tree.
