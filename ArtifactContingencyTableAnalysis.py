import math
import pandas
from ArtifactContingencyTableAnalysisUtils import ArtifactContingencyTableAnalysisUtils
from collections import OrderedDict


class ArtifactContingencyTableAnalysis(object):

    def __init__(self, alt_non_ref_pileupread_bp_count, alt_ref_pileupread_bp_count, ref_non_ref_pileupread_bp_count,
                 ref_ref_pileupread_bp_count, alt_soft_clipped_pileupread_bp_count,
                 ref_soft_clipped_pileupread_bp_count, alt_overlapping_aligned_segment_count,
                 ref_overlapping_aligned_segment_count):

        self._alt_non_ref_pileupread_bp_count = alt_non_ref_pileupread_bp_count
        self._alt_ref_pileupread_bp_count = alt_ref_pileupread_bp_count
        self._ref_non_ref_pileupread_bp_count = ref_non_ref_pileupread_bp_count
        self._ref_ref_pileupread_bp_count = ref_ref_pileupread_bp_count
        self._alt_soft_clipped_pileupread_bp_count = alt_soft_clipped_pileupread_bp_count
        self._ref_soft_clipped_pileupread_bp_count = ref_soft_clipped_pileupread_bp_count
        self._alt_overlapping_aligned_segment_count = alt_overlapping_aligned_segment_count
        self._ref_overlapping_aligned_segment_count = ref_overlapping_aligned_segment_count

        # Pre-computed values for contingency table
        self._log_total_alt_col_bp_count = math.log(alt_non_ref_pileupread_bp_count + alt_ref_pileupread_bp_count +
                                                    ArtifactContingencyTableAnalysisUtils.EPS)  # col 1 margin
        self._log_total_ref_col_bp_count = math.log(ref_non_ref_pileupread_bp_count + ref_ref_pileupread_bp_count +
                                                    ArtifactContingencyTableAnalysisUtils.EPS)  # col 2 margin
        self._log_total_non_ref_row_bp_count = math.log(alt_non_ref_pileupread_bp_count +
                                                        ref_non_ref_pileupread_bp_count +
                                                        ArtifactContingencyTableAnalysisUtils.EPS)  # row 1 margin
        self._log_total_ref_row_bp_count = math.log(alt_ref_pileupread_bp_count + ref_ref_pileupread_bp_count +
                                                    ArtifactContingencyTableAnalysisUtils.EPS)  # row 2 margin
        self._log_total_bp_count = math.log(alt_non_ref_pileupread_bp_count + alt_ref_pileupread_bp_count +
                                            ref_non_ref_pileupread_bp_count + ref_ref_pileupread_bp_count +
                                            ArtifactContingencyTableAnalysisUtils.EPS)  # total count

        alt_non_ref_soft_clipped_pileupread_bp_count = alt_non_ref_pileupread_bp_count + \
            alt_soft_clipped_pileupread_bp_count
        ref_non_ref_soft_clipped_pileupread_bp_count = ref_non_ref_pileupread_bp_count + \
            ref_soft_clipped_pileupread_bp_count
        self._alt_non_ref_soft_clipped_pileupread_bp_count = alt_non_ref_soft_clipped_pileupread_bp_count
        self._alt_ref_soft_clipped_pileupread_bp_count = alt_ref_pileupread_bp_count
        self._ref_non_ref_soft_clipped_pileupread_bp_count = ref_non_ref_soft_clipped_pileupread_bp_count
        self._ref_ref_soft_clipped_pileupread_bp_count = ref_ref_pileupread_bp_count
        self._log_total_alt_col_soft_clipped_bp_count = \
            math.log(alt_non_ref_soft_clipped_pileupread_bp_count + alt_ref_pileupread_bp_count +
                     ArtifactContingencyTableAnalysisUtils.EPS)  # col 1 margin
        self._log_total_ref_col_soft_clipped_bp_count = \
            math.log(ref_non_ref_soft_clipped_pileupread_bp_count + ref_ref_pileupread_bp_count +
                     ArtifactContingencyTableAnalysisUtils.EPS)  # col 2 margin
        self._log_total_non_ref_row_soft_clipped_bp_count = \
            math.log(alt_non_ref_soft_clipped_pileupread_bp_count + ref_non_ref_soft_clipped_pileupread_bp_count +
                     ArtifactContingencyTableAnalysisUtils.EPS)  # row 1 margin
        self._log_total_ref_row_soft_clipped_bp_count = \
            math.log(alt_ref_pileupread_bp_count + ref_ref_pileupread_bp_count +
                     ArtifactContingencyTableAnalysisUtils.EPS)  # row 2 margin
        self._log_total_soft_clipped_bp_count = \
            math.log(alt_non_ref_soft_clipped_pileupread_bp_count + alt_ref_pileupread_bp_count +
                     ref_non_ref_soft_clipped_pileupread_bp_count + ref_ref_pileupread_bp_count +
                     ArtifactContingencyTableAnalysisUtils.EPS)  # total count

    @property
    def expected_alt_non_ref_pileupread_bp_count(self):
        return math.exp(self._log_total_alt_col_bp_count + self._log_total_non_ref_row_bp_count -
                        self._log_total_bp_count)

    @property
    def expected_alt_ref_pileupread_bp_count(self):
        return math.exp(self._log_total_alt_col_bp_count + self._log_total_ref_row_bp_count - self._log_total_bp_count)

    @property
    def expected_ref_non_ref_pileupread_bp_count(self):
        return math.exp(self._log_total_ref_col_bp_count + self._log_total_non_ref_row_bp_count -
                        self._log_total_bp_count)

    @property
    def expected_ref_ref_pileupread_bp_count(self):
        return math.exp(self._log_total_ref_col_bp_count + self._log_total_ref_row_bp_count -
                        self._log_total_bp_count)

    @property
    def expected_alt_non_ref_soft_clipped_pileupread_bp_count(self):
        return math.exp(self._log_total_alt_col_soft_clipped_bp_count +
                        self._log_total_non_ref_row_soft_clipped_bp_count -
                        self._log_total_soft_clipped_bp_count)

    @property
    def expected_alt_ref_soft_clipped_pileupread_bp_count(self):
        return math.exp(self._log_total_alt_col_soft_clipped_bp_count +
                        self._log_total_ref_row_soft_clipped_bp_count -
                        self._log_total_soft_clipped_bp_count)

    @property
    def expected_ref_non_ref_soft_clipped_pileupread_bp_count(self):
        return math.exp(self._log_total_ref_col_soft_clipped_bp_count +
                        self._log_total_non_ref_row_soft_clipped_bp_count -
                        self._log_total_soft_clipped_bp_count)

    @property
    def expected_ref_ref_soft_clipped_pileupread_bp_count(self):
        return math.exp(self._log_total_ref_col_soft_clipped_bp_count +
                        self._log_total_ref_row_soft_clipped_bp_count -
                        self._log_total_soft_clipped_bp_count)

    @property
    def alt_overlapping_aligned_segment_count(self):
        return self._alt_overlapping_aligned_segment_count

    @property
    def ref_overlapping_aligned_segment_count(self):
        return self._ref_overlapping_aligned_segment_count

    @property
    def alt_non_ref_pileupread_bp_count(self):
        return self._alt_non_ref_pileupread_bp_count

    @property
    def alt_ref_pileupread_bp_count(self):
        return self._alt_ref_pileupread_bp_count

    @property
    def ref_non_ref_pileupread_bp_count(self):
        return self._ref_non_ref_pileupread_bp_count

    @property
    def ref_ref_pileupread_bp_count(self):
        return self._ref_ref_pileupread_bp_count

    @property
    def alt_soft_clipped_pileupread_bp_count(self):
        return self._alt_soft_clipped_pileupread_bp_count

    @property
    def ref_soft_clipped_pileupread_bp_count(self):
        return self._ref_soft_clipped_pileupread_bp_count

    @property
    def alt_non_ref_soft_clipped_pileupread_bp_count(self):
        return self._alt_non_ref_soft_clipped_pileupread_bp_count

    @property
    def alt_ref_soft_clipped_pileupread_bp_count(self):
        return self._alt_ref_soft_clipped_pileupread_bp_count

    @property
    def ref_non_ref_soft_clipped_pileupread_bp_count(self):
        return self._ref_non_ref_soft_clipped_pileupread_bp_count

    @property
    def ref_ref_soft_clipped_pileupread_bp_count(self):
        return self._ref_ref_soft_clipped_pileupread_bp_count

    @staticmethod
    def retrieve_fields_as_series(dataTable, prefix):
        data = pandas.Series()

        data["alt_non_ref_pileupread_bp_count"] = dataTable.alt_non_ref_pileupread_bp_count
        data["alt_ref_pileupread_bp_count"] = dataTable.alt_ref_pileupread_bp_count
        data["expected_alt_non_ref_pileupread_bp_count"] = dataTable.expected_alt_non_ref_pileupread_bp_count
        data["expected_alt_ref_pileupread_bp_count"] = dataTable.expected_alt_ref_pileupread_bp_count
        data["alt_non_ref_soft_clipped_pileupread_bp_count"] = dataTable.alt_non_ref_soft_clipped_pileupread_bp_count
        data["alt_ref_soft_clipped_pileupread_bp_count"] = dataTable.alt_ref_soft_clipped_pileupread_bp_count
        data["expected_alt_non_ref_soft_clipped_pileupread_bp_count"] = \
            dataTable.expected_alt_non_ref_soft_clipped_pileupread_bp_count
        data["expected_alt_ref_soft_clipped_pileupread_bp_count"] = \
            dataTable.expected_alt_ref_soft_clipped_pileupread_bp_count
        data["alt_soft_clipped_pileupread_bp_count"] = dataTable.alt_soft_clipped_pileupread_bp_count
        data["alt_overlapping_aligned_segment_count"] = dataTable.alt_overlapping_aligned_segment_count
        data["ref_non_ref_pileupread_bp_count"] = dataTable.ref_non_ref_pileupread_bp_count
        data["ref_ref_pileupread_bp_count"] = dataTable.ref_ref_pileupread_bp_count
        data["expected_ref_non_ref_pileupread_bp_count"] = dataTable.expected_ref_non_ref_pileupread_bp_count
        data["expected_ref_ref_pileupread_bp_count"] = dataTable.expected_ref_ref_pileupread_bp_count
        data["ref_non_ref_soft_clipped_pileupread_bp_count"] = dataTable.ref_non_ref_soft_clipped_pileupread_bp_count
        data["ref_ref_soft_clipped_pileupread_bp_count"] = dataTable.ref_ref_soft_clipped_pileupread_bp_count
        data["expected_ref_non_ref_soft_clipped_pileupread_bp_count"] = \
            dataTable.expected_ref_non_ref_soft_clipped_pileupread_bp_count
        data["expected_ref_ref_soft_clipped_pileupread_bp_count"] = \
            dataTable.expected_ref_ref_soft_clipped_pileupread_bp_count
        data["ref_soft_clipped_pileupread_bp_count"] = dataTable.ref_soft_clipped_pileupread_bp_count
        data["ref_overlapping_aligned_segment_count"] = dataTable.ref_overlapping_aligned_segment_count

        index = dict([(name, prefix + name) for name in data.index.tolist()])  # add prefix

        return data.rename(index=index)

    @staticmethod
    def render_contingency_table(dataTable):
        return [[dataTable.alt_non_ref_pileupread_bp_count, dataTable.ref_non_ref_pileupread_bp_count],
                [dataTable.alt_ref_pileupread_bp_count, dataTable.ref_ref_pileupread_bp_count]]

    @staticmethod
    def render_soft_clipped_contingency_table(dataTable):
        return [[dataTable.alt_non_ref_soft_clipped_pileupread_bp_count, dataTable.ref_non_ref_soft_clipped_pileupread_bp_count],
                [dataTable.alt_ref_soft_clipped_pileupread_bp_count, dataTable.ref_ref_soft_clipped_pileupread_bp_count]]

    @classmethod
    def create(cls, ref_allele, alt_allele, pileupcolumn, pileupcolumn_knapsack, pileupcolumn_mask=None):

        alt_non_ref_pileupread_bp_count = 0
        alt_ref_pileupread_bp_count = 0

        ref_non_ref_pileupread_bp_count = 0
        ref_ref_pileupread_bp_count = 0

        alt_soft_clipped_pileupread_bp_count = 0
        ref_soft_clipped_pileupread_bp_count = 0

        ref_overlapping_aligned_segment_count = 0
        alt_overlapping_aligned_segment_count = 0

        pileupread_alignment_query_names = []
        ref_supporting_pileupread_alignment_query_names = []
        alt_supporting_pileupread_alignment_query_names = []

        ref_alleles = ArtifactContingencyTableAnalysisUtils.retrieve_ref_alleles(pileupcolumn_knapsack)
        num_reads = 0
        num_ref_supporting_reads = 0
        num_alt_supporting_reads = 0

        for pileupread in pileupcolumn.pileups:
            num_reads += 1
            pileupread_alignment_query_name = pileupread.alignment.query_name

            is_indel = ArtifactContingencyTableAnalysisUtils.pileupread_has_indel(pileupread)  # not used
            is_ref = ArtifactContingencyTableAnalysisUtils.pileupread_has_ref_allele(pileupread, ref_allele)
            is_alt = ArtifactContingencyTableAnalysisUtils.pileupread_has_alt_allele(pileupread, alt_allele)

            if not is_indel and is_ref and not is_alt and \
                    pileupread_alignment_query_name in ref_supporting_pileupread_alignment_query_names:
                ref_overlapping_aligned_segment_count += 1  # ref supporting mate pairs
                pass  # do nothing
            elif not is_indel and not is_ref and is_alt and \
                    pileupread_alignment_query_name in alt_supporting_pileupread_alignment_query_names:
                alt_overlapping_aligned_segment_count += 1  # alt supporting mate pairs
                pass  # do nothing
            elif (not is_indel and is_ref and not is_alt and
                    pileupread_alignment_query_name in alt_supporting_pileupread_alignment_query_names) or \
                    (not is_indel and not is_ref and is_alt and
                        pileupread_alignment_query_name in ref_supporting_pileupread_alignment_query_names):
                pass  # alt/ref supporting mate pairs
            elif not is_indel and is_ref and not is_alt:  # ref supporting read
                pileupread_alignment_query_names += [pileupread_alignment_query_name]
                ref_supporting_pileupread_alignment_query_names += [pileupread_alignment_query_name]
                num_ref_supporting_reads += 1
            elif not is_indel and not is_ref and is_alt:  # alt supporting read
                pileupread_alignment_query_names += [pileupread_alignment_query_name]
                alt_supporting_pileupread_alignment_query_names += [pileupread_alignment_query_name]
                num_alt_supporting_reads += 1
            else:  # indel based pivots and non-biallelic leftovers
                pass  # do nothing

        # TODO: what happens when the read is ref supporting but it's mate is alt supporting?
        # ref_supporting_pileupread_alignment_query_names = set(ref_supporting_pileupread_alignment_query_names)
        # alt_supporting_pileupread_alignment_query_names = set(alt_supporting_pileupread_alignment_query_names)

        indexed_pileupreads = \
            ArtifactContingencyTableAnalysisUtils.retrieve_indexed_pileupreads(pileupcolumn_knapsack,
                                                                               pileupread_alignment_query_names,
                                                                               pileupcolumn_mask)

        # ref supporting counts
        # ref_ref_pileupread_bp_count = \
        #     ArtifactContingencyTableAnalysisUtils.retrieve_ref_bp_count(indexed_pileupreads, ref_alleles,
        #         ref_supporting_pileupread_alignment_query_names)

        # for pileupcolumn_name in indexed_pileupreads:  # iterate over columns
        #     pileupreads = indexed_pileupreads[pileupcolumn_name]
        #     for pileupread in pileupreads:  # iterate over rows
        #         if pileupread.alignment.query_name in ref_supporting_pileupread_alignment_query_names:
        #             # if pileupread.alignment.query_name.startswith("C0CPNACXX120126:8:1101"):
        #             #     print pileupcolumn_name
        #             print "%s:%s:%s" % (pileupread.alignment.query_name, pileupread.alignment.cigarstring,
        #                                 pileupread.alignment.get_tag("NM"))

        ref_non_ref_pileupread_bp_count = \
            ArtifactContingencyTableAnalysisUtils.retrieve_non_ref_bp_count(indexed_pileupreads, ref_alleles,
                ref_supporting_pileupread_alignment_query_names)

        # for pileupcolumn_name in indexed_pileupreads:  # iterate over columns
        #     pileupreads = indexed_pileupreads[pileupcolumn_name]
        #     for pileupread in pileupreads:  # iterate over rows
        #         if pileupread.alignment.query_name in ref_supporting_pileupread_alignment_query_names:
        #             # if pileupread.alignment.query_name.startswith("C0CPNACXX120126:8:1101"):
        #             #     print pileupcolumn_name
        #             print "%s:%s:%s" % (pileupread.alignment.query_name, pileupread.alignment.cigarstring,
        #                                 pileupread.alignment.get_tag("NM"))


        # # alt supporting counts
        # alt_ref_pileupread_bp_count = \
        #     ArtifactContingencyTableAnalysisUtils.retrieve_ref_bp_count(indexed_pileupreads, ref_alleles,
        #         alt_supporting_pileupread_alignment_query_names)
        # alt_non_ref_pileupread_bp_count = \
        #     ArtifactContingencyTableAnalysisUtils.retrieve_non_ref_bp_count(indexed_pileupreads, ref_alleles,
        #         alt_supporting_pileupread_alignment_query_names)

        import sys
        sys.exit()
        # return ArtifactContingencyTableAnalysis(alt_non_ref_pileupread_bp_count=alt_non_ref_pileupread_bp_count,
        #                                         alt_ref_pileupread_bp_count=alt_ref_pileupread_bp_count,
        #                                         ref_non_ref_pileupread_bp_count=ref_non_ref_pileupread_bp_count,
        #                                         ref_ref_pileupread_bp_count=ref_ref_pileupread_bp_count,
        #                                         alt_soft_clipped_pileupread_bp_count=alt_soft_clipped_pileupread_bp_count,
        #                                         ref_soft_clipped_pileupread_bp_count=ref_soft_clipped_pileupread_bp_count,
        #                                         alt_overlapping_aligned_segment_count=alt_overlapping_aligned_segment_count,
        #                                         ref_overlapping_aligned_segment_count=ref_overlapping_aligned_segment_count)
