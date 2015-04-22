from collections import OrderedDict
from PileupColumnUtils import PileupColumnUtils
from DataType import DataType
from scipy.stats import binom


class PileupColumnMask(object):

    def __init__(self, mask):
        self._mask = mask

    @property
    def mask(self):
        return self._mask

    @staticmethod
    def create(pileupcolumn_knapsack, data_type=DataType.germline, nearby_somatic_mutations=None):
        nearby_somatic_mutations = [] if nearby_somatic_mutations is None else nearby_somatic_mutations
        mask = OrderedDict()

        if data_type == DataType.germline:  # germline data
            for pileupcolumn_name in pileupcolumn_knapsack.pileupread_knapsacks:
                mask[pileupcolumn_name] = False
                pileupread_knapsack = pileupcolumn_knapsack.pileupread_knapsacks[pileupcolumn_name]
                if pileupread_knapsack.non_ref_allele_count != 0:
                    if binom.cdf(pileupread_knapsack.non_ref_allele_count,
                                 pileupread_knapsack.total_allele_count, 0.45) > 0.05:
                        mask[pileupcolumn_name] = True
        elif data_type == DataType.somatic:  # somatic data
            # mask nearby somatic mutations (may include DNPs, TNPs, etc.)
            mask = OrderedDict([(pileupcolumn_name, False)
                                for pileupcolumn_name in pileupcolumn_knapsack.pileupcolumn_names])
            for nearby_somatic_mutation in nearby_somatic_mutations:
                pileupcolumn_name = PileupColumnUtils.retrieve_pileupcolumn_name(nearby_somatic_mutation.chrom,
                                                                                 nearby_somatic_mutation.start,
                                                                                 nearby_somatic_mutation.end)
                mask[pileupcolumn_name] = True
        else:
            mask = OrderedDict([(pileupcolumn_name, False)
                                for pileupcolumn_name in pileupcolumn_knapsack.pileupcolumn_names])

            # TODO: using aggregate counts, determine whether the base pair should be masked or not
            # TODO: ensure that mask for insertions is at the base before and after (used later)
            # TODO: ensure that the mask for deletions is covering all the bases
            # TODO: mask will only either hide the entire deletion or that it will hide part of it
            # TODO: correct counts for the adjacent reads

        return PileupColumnMask(mask)