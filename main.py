import sys, math
from enums import *
from Bio import AlignIO
from Bio.Alphabet import generic_protein, generic_dna
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

class Computation:

    def __init__(self, options=None):
        self.algorithm = options["algorithm"]
        self.includeGaps = options["include_gaps"] # includes gaps into the tree computation
        self.removePoor = options["remove_poor"] # removes poorly conserved regions
        self.gapCutoff = options["gap_cutoff"] # at least x% in column are not gaps
        self.pairCutoff = options["pair_cutoff"] # at least x% of pairs in column are identical

        self.seqType = None
        self.seqType = options["seq_type"]
        if self.seqType == SeqTypes.AA:
            self.aligmnment = AlignIO.read(options["ali_path"], "fasta", alphabet=generic_protein)
        if self.seqType == SeqTypes.DNA:
            self.aligmnment = AlignIO.read(options["ali_path"], "fasta", alphabet=generic_dna)
        assert self.seqType

        self.aligmnment = self.cleanAlignment(self.aligmnment)

    def cleanAlignment(self, alignment):
        to_remove = set()
        for col_idx in range(alignment.get_alignment_length()):
            column = alignment[:,col_idx]

            if not self.includeGaps:
                if column.find("-") != -1:
                    to_remove.add(col_idx)
                    continue

            if self.removePoor:
                # remove column if it contains too many gaps
                nongap_ratio = float(len(column) - column.count("-")) / len(column)
                if nongap_ratio < self.gapCutoff:
                    to_remove.add(col_idx)
                    continue

                # remove column if there are not enough identical residues
                all_pairs = math.factorial(len(column) - column.count("-")) / (2.0 * math.factorial((len(column) - column.count("-")) - 2))
                identical_pairs = 0
                for res in set(column):
                    if res != "-":
                        res_count = column.count(res)
                        if res_count > 1:
                            identical_pairs += math.factorial(res_count) / (2.0 * math.factorial(res_count - 2))
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

def main(args):
    options = {
        "algorithm" : TBAlgorithms.NJ
        , "seq_type" : SeqTypes.AA
        , "ali_path" : "testfiles/keratins_ali.fasta"
        , "include_gaps" : True
        , "remove_poor" : True
        , "gap_cutoff" : 0.8
        , "pair_cutoff" : 0.3
    }
    computation = Computation(options)


if __name__ == "__main__":
    main(sys.argv)


