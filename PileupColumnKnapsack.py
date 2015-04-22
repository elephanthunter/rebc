from PileupColumnUtils import PileupColumnUtils
from PileupReadKnapsack import PileupReadKnapsack
from collections import OrderedDict


class PileupColumnKnapsack():

    def __init__(self, somatic_mutation, pileupread_knapsacks):
        self._somatic_mutation = somatic_mutation
        self._pileupread_knapsacks = pileupread_knapsacks

    @property
    def pileupread_knapsacks(self):
        return self._pileupread_knapsacks

    @property
    def pileupcolumn_names(self):
        return self._pileupread_knapsacks.keys()

    def retrieve_pileupread(self, pileupcolumn_name, alignment_query_name):
        pileupread = None
        try:
            pileupread_knapsack = self._pileupread_knapsacks[pileupcolumn_name]
            return pileupread_knapsack.pileupreads.get(alignment_query_name, None)
        except KeyError:
            return pileupread

    def retrieve_ref_allele(self, pileupcolumn_name):
        ref_allele = None
        if pileupcolumn_name in self._pileupread_knapsacks:
            return self._pileupread_knapsacks[pileupcolumn_name].ref_allele
        return ref_allele

    @classmethod
    #@profile
    def create(cls, somatic_mutation, sample_bam_file, ref_seq_file):
        pileupread_knapsacks = OrderedDict()
        prev_right_soft_clipped = []
        prev_baseknapsack = None
        for pileupcolumn in sample_bam_file.pileup(somatic_mutation.chrom, somatic_mutation.start, somatic_mutation.end,
                                                   truncate=False):
            pileupcolumn_start = pileupcolumn.pos
            pileupcolumn_end = pileupcolumn.pos+1
            pileupcolumn_name = PileupColumnUtils.retrieve_pileupcolumn_name(somatic_mutation.chrom, pileupcolumn_start,
                                                                             pileupcolumn_end)
            if pileupcolumn.pos != somatic_mutation.start:
                ref_allele = ref_seq_file.fetch(somatic_mutation.chrom, pileupcolumn_start, pileupcolumn_end)  # [start,end) region
                pileupread_knapsacks[pileupcolumn_name] = \
                    PileupReadKnapsack.create(somatic_mutation.chrom, pileupcolumn_start, pileupcolumn_end, ref_allele,
                                              pileupcolumn)

                # Re-adjust left and right soft clipped
                baseknapsack = pileupread_knapsacks[pileupcolumn_name].baseknapsack
                if not prev_baseknapsack:
                    baseknapsack["left_soft_clipped"] = []  # dropped
                    prev_right_soft_clipped = baseknapsack["right_soft_clipped"]
                    baseknapsack["right_soft_clipped"] = []  # cleaned
                else:
                    prev_baseknapsack["right_soft_clipped"] = baseknapsack["left_soft_clipped"]
                    baseknapsack["left_soft_clipped"] = prev_right_soft_clipped
                    prev_right_soft_clipped = baseknapsack["right_soft_clipped"]
                    baseknapsack["right_soft_clipped"] = []  # cleaned
                prev_baseknapsack = baseknapsack

                # TODO: adjust for the presence of collapsed soft clipped reads

        return PileupColumnKnapsack(somatic_mutation, pileupread_knapsacks)