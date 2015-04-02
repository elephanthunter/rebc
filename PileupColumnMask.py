from collections import OrderedDict
from PileupReadKnapsack import PileupReadKnapsack
from BasePairUtils import BasePairUtils


class PileupColumnMask(object):

    def __init__(self, mask):
        self._mask = mask

    @property
    def mask(self):
        return self._mask

    def is_pileupcolumn_masked(self, pileupcolumn_name):
        return ~(pileupcolumn_name in self._mask) | self._mask[pileupcolumn_name]

    @staticmethod
    def create(pileupcolumn_knapsack):
        mask = OrderedDict()
        prev_pileupcolumn_name = None
        for pileupcolumn_name in pileupcolumn_knapsack.pileupread_knapsacks:
            # determine what is alt at this position

            # aggregate_counts = pileupread_knapsack.base_pair_aggregate_counts
            # TODO: using aggregate counts, determine whether the base pair should be masked or not
            # TODO: ensure that mask for insertions is at the base before and after (used later)
            # TODO: ensure that the mask for deletions is covering all the bases
            # TODO: mask will only either hide the entire deletion or that it will hide part of it

            # TODO: correct counts for the adjacent reads
        # For now, set all of them to True
            mask[pileupcolumn_name] = False

        return PileupColumnMask(mask)