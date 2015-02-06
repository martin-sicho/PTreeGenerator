import math

from Bio import AlignIO
from Bio.Alphabet import generic_protein, generic_dna
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from ptreegen.enums import *
from NeigborJoining import NeigborJoining
import distance_functions as dfuncs

class Computation:

    def __init__(self, options):
        self.algorithm = options["algorithm"]
        self.gapPenalty = None
        self.gapPenalty = options["gap_penalty"] # cost of gaps when they are included to the distance computation
        self.includeGaps = options["include_gaps"] # includes gaps into the tree computation
        self.removePoor = options["remove_poor"] # removes poorly conserved regions
        self.gapCutoff = options["gap_cutoff"] # at least x% in column are not gaps
        self.pairCutoff = options["pair_cutoff"] # at least x% of pairs in column are identical

        self.seqType = None
        self.seqType = options["seq_type"]
        if self.seqType == SeqTypes.AA:
            self.alignment = AlignIO.read(options["ali_path"], "fasta", alphabet=generic_protein)
        elif self.seqType == SeqTypes.DNA:
            self.alignment = AlignIO.read(options["ali_path"], "fasta", alphabet=generic_dna)
        else:
            raise RuntimeError("Unknown sequence type: " + self.seqType)

        self.alignment = self.cleanAlignment(self.alignment)
        self.distanceMatrix = self.computeDistanceMatrix(self.alignment, dfuncs.p_distance)

        self.tree = self.computeTree()

    def computeTree(self):
        if self.algorithm == TBAlgorithms.NJ:
            return NeigborJoining(self.distanceMatrix, [x.id for x in self.alignment])()
        else:
            raise RuntimeError(self.algorithm + "not implemented.")

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