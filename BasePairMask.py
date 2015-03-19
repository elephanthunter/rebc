from collections import OrderedDict
from PileupReadKnapsack import PileupReadKnapsack
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
    def create(pileupread_knapsacks):
        base_pair_read_knapsacks = OrderedDict()
        mask = OrderedDict()
        for pileupcolumn_name in pileupread_knapsacks.pileupread_knapsacks:
            base_pair_read_knapsack = pileupread_knapsacks.pileupread_knapsacks[pileupcolumn_name]
                # aggregate_counts = pileupread_knapsack.base_pair_aggregate_counts
                # TODO: Using aggregate counts, determine whether the base pair should be masked or not
                # For now, set all of them to True
            mask[pileupcolumn_name] = False

        return BasePairMask(mask, base_pair_read_knapsacks)