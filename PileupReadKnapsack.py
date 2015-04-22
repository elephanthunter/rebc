import math
from collections import OrderedDict


class PileupReadKnapsack(object):

    # TODO: What to do about overlapping reads? Because I am using a dictionary, it should not matter.
    def __init__(self, chrom, start, end, ref_allele, pileupreads):
        self._chrom = chrom
        self._start = start
        self._end = end
        self._ref_allele = ref_allele
        self._pileupreads = pileupreads
        self._pileupread_alignment_query_names = set(pileupreads.keys())
        self._descriptor = PileupReadKnapsack.PileupColumnDescriptor.create(position=start, ref_allele=ref_allele,
                                                                            pileupreads=pileupreads)

    @property
    def baseknapsack(self):
        return self._descriptor.baseknapsack

    @property
    def pileupread_alignment_query_names(self):
        return self._pileupread_alignment_query_names

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def ref_allele(self):
        return self._ref_allele

    @property
    def pileupreads(self):
        return self._pileupreads

    @property
    def total_allele_count(self):
        # Does not include soft clipped bases
        return self._descriptor.num_indels + self._descriptor.num_ref_snps + self._descriptor.num_non_ref_snps

    @property
    def non_ref_allele_count(self):
        # Does not include soft clipped bases
        return self._descriptor.num_indels + self._descriptor.num_non_ref_snps

    @classmethod
    def create(cls, chrom, start, end, ref_allele, pileupcolumn):
        pileupreads = OrderedDict()  # ordered list of pileup reads corresponding to the
        for pileupread in pileupcolumn.pileups:  # won't there be an issue with the soft clipped region?
            pileupreads[pileupread.alignment.query_name] = pileupread
        return PileupReadKnapsack(chrom=chrom, start=start, end=end, ref_allele=ref_allele, pileupreads=pileupreads)

    class PileupColumnDescriptor(object):

        def __init__(self, ref_allele, baseknapsack):
            self._ref_allele = ref_allele
            self._baseknapsack = baseknapsack

        @property
        def num_non_ref_snps(self):
            return self.num_ade + self.num_thy + self.num_cyt + self.num_gua - self.num_ref_snps

        @property
        def num_indels(self):
            return self.num_del + self.num_ins

        @property
        def num_ref_snps(self):
            return len(self._baseknapsack[self._ref_allele])

        @property
        def num_cyt(self):
            return len(self._baseknapsack["C"])

        @property
        def num_ade(self):
            return len(self._baseknapsack["A"])

        @property
        def num_gua(self):
            return len(self._baseknapsack["G"])

        @property
        def num_thy(self):
            return len(self._baseknapsack["T"])

        @property
        def num_ins(self):
            return len(self._baseknapsack["ins"])

        @property
        def num_del(self):
            return len(self._baseknapsack["del"])

        @property
        def num_soft_clipped(self):
            return len(self._baseknapsack["left_soft_clipped"]) + len(self._baseknapsack["right_soft_clipped"])

        @property
        def baseknapsack(self):
            return self._baseknapsack

        @classmethod
        def create(cls, position, ref_allele, pileupreads):
            baseknapsack = OrderedDict()
            baseknapsack["A"] = []
            baseknapsack["T"] = []
            baseknapsack["C"] = []
            baseknapsack["G"] = []
            baseknapsack["ins"] = []
            baseknapsack["del"] = []
            baseknapsack["left_soft_clipped"] = []
            baseknapsack["right_soft_clipped"] = []

            for pileupread_alignment_query_name in pileupreads:
                pileupread = pileupreads[pileupread_alignment_query_name]  # dictionary (overlapping reads are added
                                                                           # only once)
                if not pileupread.is_del and pileupread.indel == 0 and "S" in pileupread.alignment.cigarstring:
                    position_index = pileupread.alignment.positions.index(position)
                    if position_index == 0 and pileupread.alignment.cigarstring.index("S") < \
                            int(math.ceil(len(pileupread.alignment.cigarstring)/2.0)):
                        baseknapsack["left_soft_clipped"] += [pileupread_alignment_query_name]
                        continue
                    elif position_index == len(pileupread.alignment.positions)-1 and \
                            pileupread.alignment.cigarstring.index("S") >= len(pileupread.alignment.cigarstring)-1:
                        baseknapsack["right_soft_clipped"] += [pileupread_alignment_query_name]
                        continue

                if pileupread.indel > 0:  # insertion
                    baseknapsack["ins"] += [pileupread_alignment_query_name]
                elif pileupread.is_del:  # deletion
                    baseknapsack["del"] += [pileupread_alignment_query_name]
                else:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    if base in ("A", "T", "C", "G",):
                        baseknapsack[base] += [pileupread_alignment_query_name]

            return PileupReadKnapsack.PileupColumnDescriptor(ref_allele, baseknapsack)