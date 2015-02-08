from Bio import AlignIO
from Bio.Alphabet import generic_protein, generic_dna, generic_rna
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from ptreegen.enums import *
from neigbor_joining import NeigborJoining
import distance_functions as dfuncs


class Computation:

    def __init__(self, options):
        self.algorithm = None
        self.gapPenalty = None
        self.gapPenalty = None
        self.includeGaps = None
        self.removePoor = None
        self.gapCutoff = None
        self.pairCutoff = None
        self.seqType = None
        self.parseOptions(options)
        self.alignment = self.cleanAlignment(self.alignment)
        self.distanceMatrix = self.computeDistanceMatrix(self.alignment, dfuncs.p_distance)
        self.tree = self.computeTree()

    def parseOptions(self, options):
        self.algorithm = options["method"]
        if self.algorithm not in (TreeBuildAlgorithms.NJ, TreeBuildAlgorithms.PARSIMONY):
            raise RuntimeError("Unknown method: " + self.algorithm)
        self.gapPenalty = options["gap_penalty"] # cost of gaps when they are included to the distance computation
        if self.gapPenalty < 0 or self.gapPenalty > 1:
            raise RuntimeError("Bad gap penalty value. Must be between 0 and 1. Got: " + self.gapPenalty)
        self.includeGaps = not options["no_gaps"] # includes gaps into the tree computation
        self.removePoor = not options["no_cleaning"] # removes poorly conserved regions
        self.gapCutoff = options["gap_cutoff"] # at least x% in column are not gaps
        if self.gapCutoff < 0 or self.gapCutoff > 1:
            raise RuntimeError("Bad gap cutoff value. Must be between 0 and 1. Got: " + self.gapCutoff)
        self.pairCutoff = options["pair_cutoff"] # at least x% of pairs in column are identical
        self.seqType = options["sequence_type"]
        if self.seqType == SeqTypes.AA:
            self.alignment = AlignIO.read(options["alignment_file"], "fasta", alphabet=generic_protein)
        elif self.seqType == SeqTypes.DNA:
            self.alignment = AlignIO.read(options["alignment_file"], "fasta", alphabet=generic_dna)
        elif self.seqType == SeqTypes.RNA:
            self.alignment = AlignIO.read(options["alignment_file"], "fasta", alphabet=generic_rna)
        else:
            raise RuntimeError("Unknown sequence type: " + self.seqType)

    def computeTree(self):
        if self.algorithm == TreeBuildAlgorithms.NJ:
            return NeigborJoining(self.distanceMatrix, [x.id for x in self.alignment])()
        else:
            raise RuntimeError(self.algorithm + " not implemented.")

    def update(self):
        self.alignment = self.cleanAlignment(self.alignment)
        self.distanceMatrix = self.computeDistanceMatrix(self.alignment, dfuncs.p_distance)
        self.tree = self.computeTree()

    def computeDistanceMatrix(self, alignment, distFunction):
        dist_matrix = []
        for i,record_i in enumerate(alignment):
            seq_i = record_i.seq
            distances = []
            for j,record_j in enumerate(alignment):
                if j > i:
                    seq_j = record_j.seq
                    if self.gapPenalty:
                        distances.append(distFunction(seq_i, seq_j, gapPenalty=self.gapPenalty))
                    else:
                        distances.append(distFunction(seq_i, seq_j))
                elif i == j:
                    distances.append(0)
                else:
                    distances.append(dist_matrix[j][i])
            dist_matrix.append(tuple(distances))
        return tuple(dist_matrix)

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