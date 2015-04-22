import pandas


class SomaticMutation(object):
    def __init__(self, chrom, start, end, ref_allele, alt_allele):
        self._chrom = chrom
        self._start = start
        self._end = end
        self._ref_allele = ref_allele
        self._alt_allele = alt_allele

    def __eq__(self, otr):
        return self._chrom == otr.chrom and self._start == otr.start and self._end == otr.end \
            and self._ref_allele == otr.ref_allele and self._alt_allele == otr.alt_allele

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
    def alt_allele(self):
        return self._alt_allele

    @classmethod
    def create(cls, row):
        return SomaticMutation(chrom=str(row["Chromosome"]), start=int(row["Start_position"])-1,  # 0-based indexing
                               end=int(row["End_position"]), ref_allele=str(row["Reference_Allele"]),
                               alt_allele=str(row["Tumor_Seq_Allele2"]))

    def retrieve_as_series(self):
        return pandas.Series(data=[self._chrom, self._start, self._end, self._ref_allele, self._alt_allele],
                             index=["Chromosome", "Start_position", "End_position", "Reference_Allele",
                                    "Tumor_Seq_Allele2"])