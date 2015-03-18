from collections import OrderedDict
from BasePairReadKnapsack import BasePairReadKnapsack
from BasePairUtils import BasePairUtils


class BasePairMask(object):

    def __init__(self, mask, base_pair_read_knapsacks):
        self._base_pair_mask = mask
        self._base_pair_read_knapsacks = base_pair_read_knapsacks

    @property
    def mask(self):
        return self._base_pair_mask

    @property
    def base_pair_read_knapsacks(self):
        return self._base_pair_read_knapsacks

    def is_position_masked(self, position_name):
        return self._base_pair_mask[position_name]

    @staticmethod
    def create(chrom, start, end, bam_file, ref_seq_file):
        base_pair_read_knapsacks = OrderedDict()
        mask = OrderedDict()
        for pileupcolumn in bam_file.pileup(chrom, start, end, truncate=False):
            if pileupcolumn.pos != start:
                base_pair_read_knapsack = BasePairReadKnapsack.create(chrom=chrom, start=pileupcolumn.pos,
                                                                      end=pileupcolumn.pos+1, pileupcolumn=pileupcolumn,
                                                                      refseqfile=ref_seq_file)
                base_pair_read_knapsacks[base_pair_read_knapsack.name] = base_pair_read_knapsack
                # aggregate_counts = base_pair_read_knapsack.retrieve_base_pair_aggregate_counts()

                # TODO: Using aggregate counts, determine whether the base pair should be masked or not
                # For now, set all of them to True

                mask[base_pair_read_knapsack.name] = False

        return BasePairMask(mask, base_pair_read_knapsacks)