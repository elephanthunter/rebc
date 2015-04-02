from collections import OrderedDict


class PileupReadKnapsack(object):

    # TODO: What to do about overlapping reads? Because I am using a dictionary, it should not matter.
    def __init__(self, chrom, start, end, ref_allele, pileupreads):
        self._chrom = chrom
        self._start = start
        self._end = end
        self._ref_allele = ref_allele
        self._pileupreads = pileupreads
        self._baseknapsack = PileupReadKnapsack.PileupColumnDescriptor.create(position=start, pileupreads=pileupreads)

    @property
    def pileupread_alignment_query_names(self):
        return self._pileupreads.keys()

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

    @classmethod
    def create(cls, chrom, start, end, ref_allele, pileupcolumn):
        pileupreads = OrderedDict()  # ordered list of pileup reads corresponding to the
        for pileupread in pileupcolumn.pileups:  # won't there be an issue with the soft clipped region?
            pileupreads[pileupread.alignment.query_name] = pileupread
        return PileupReadKnapsack(chrom=chrom, start=start, end=end, ref_allele=ref_allele, pileupreads=pileupreads)

    class PileupColumnDescriptor(object):

        def __init__(self, baseknapsack):
            self._baseknapsack = baseknapsack

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
            return len(self._baseknapsack["soft_clipped"])

        @property
        def baseknapsack(self):
            return self._baseknapsack

        @staticmethod
        def _is_edge_soft_clipped(position, pileupread):
            index = pileupread.alignment.positions.index(position)
            return index == 0 or index == len(pileupread.alignment.positions)-1

        @classmethod
        def create(cls, position, pileupreads):
            baseknapsack = OrderedDict()
            baseknapsack["A"] = []
            baseknapsack["T"] = []
            baseknapsack["C"] = []
            baseknapsack["G"] = []
            baseknapsack["ins"] = []
            baseknapsack["del"] = []
            baseknapsack["soft_clipped"] = []

            for pileupread_alignment_query_name in pileupreads:
                pileupread = pileupreads[pileupread_alignment_query_name]  # dictionary (overlapping reads are added
                                                                           # only once)
                if not pileupread.is_del and pileupread.indel == 0 and "S" in pileupread.alignment.cigarstring:
                    if PileupReadKnapsack.PileupColumnDescriptor._is_edge_soft_clipped(position, pileupread):
                        baseknapsack["soft_clipped"] += [pileupread_alignment_query_name]
                        continue

                if pileupread.indel > 0:  # insertion
                    baseknapsack["ins"] += [pileupread_alignment_query_name]
                elif pileupread.is_del:  # deletion
                    baseknapsack["del"] += [pileupread_alignment_query_name]
                else:
                    readbase = pileupread.alignment.query_sequence[pileupread.query_position]
                    if readbase in ("A", "T", "C", "G",):
                        baseknapsack[readbase] += [pileupread_alignment_query_name]

            return PileupReadKnapsack.PileupColumnDescriptor(baseknapsack)