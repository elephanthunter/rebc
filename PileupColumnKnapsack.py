from BasePairUtils import BasePairUtils
from PileupReadKnapsack import PileupReadKnapsack
from collections import OrderedDict


class PileupColumnKnapsack():

    def __init__(self, pileupread_knapsacks):
        self._pileupread_knapsacks = pileupread_knapsacks

    @property
    def pileupread_knapsacks(self):
        return self._pileupread_knapsacks

    @property
    def pileupcolumn_names(self):
        return self._pileupread_knapsacks.keys()

    def retrieve_pileupread(self, pileupcolumn_name, alignment_query_name):
        pileupread = None
        if pileupcolumn_name in self._pileupread_knapsacks:
            pileupread_knapsack = self._pileupread_knapsacks[pileupcolumn_name]
            if alignment_query_name in pileupread_knapsack.pileupread_alignment_query_names:
                pileupreads = pileupread_knapsack.pileupreads
                pileupread = pileupreads[alignment_query_name]
        return pileupread

    @classmethod
    def create(cls, chrom, start, end, sample_bam_file, ref_seq_file):
        pileupread_knapsacks = OrderedDict()

        for pileupcolumn in sample_bam_file.pileup(chrom, start, end, truncate=False):
            pileupcolumn_chrom = chrom
            pileupcolumn_start = pileupcolumn.pos
            pileupcolumn_end = pileupcolumn.pos+1

            pileupcolumn_name = BasePairUtils.retrieve_pileupcolumn_name(pileupcolumn_chrom, pileupcolumn_start,
                                                                         pileupcolumn_end)
            if pileupcolumn.pos != start:
                ref_allele = ref_seq_file.fetch(chrom, pileupcolumn_start, pileupcolumn_end)  # [start,end) region is called
                pileupread_knapsacks[pileupcolumn_name] = \
                    PileupReadKnapsack.create(pileupcolumn_chrom, pileupcolumn_start, pileupcolumn_end, ref_allele,
                                              pileupcolumn)

        return PileupColumnKnapsack(pileupread_knapsacks=pileupread_knapsacks)