from BasePairUtils import BasePairUtils
from BasePairReadKnapsack import BasePairReadKnapsack
from collections import OrderedDict


class PileupReadKnapsack():

    def __init__(self, base_pair_read_knapsacks):
        self._base_pair_read_knapsacks = base_pair_read_knapsacks

    @property
    def base_pair_read_knapsacks(self):
        return self._base_pair_read_knapsacks

    @property
    def pileupcolumn_names(self):
        return self._base_pair_read_knapsacks.keys()

    @classmethod
    def create(cls, chrom, start, end, sample_bam_file, ref_seq_file):
        base_pair_read_knapsacks = OrderedDict()

        for pileupcolumn in sample_bam_file.pileup(chrom, start, end, truncate=False):
            pileupcolumn_name = BasePairUtils.retrieve_position_name(chrom, start, end)
            if pileupcolumn.pos != start:
                ref_allele = ref_seq_file.fetch(chrom, start=pileupcolumn.pos,
                                                end=pileupcolumn.pos+1)  # [start,end) region is called
                base_pair_read_knapsack = BasePairReadKnapsack.create(chrom=chrom, start=pileupcolumn.pos,
                                                                      end=pileupcolumn.pos+1, pileupcolumn=pileupcolumn,
                                                                      ref_allele=ref_allele)
                base_pair_read_knapsacks[pileupcolumn_name] = base_pair_read_knapsack

        return PileupReadKnapsack(base_pair_read_knapsacks=base_pair_read_knapsacks)