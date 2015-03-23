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

        pileupreads = OrderedDict()
        ref_supporting_pileupread_alignment_query_name = []
        alt_supporting_pileupread_alignment_query_name = []

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

            print "%s:%s:%s:%s:%s:%s:%s" % \
                  (pileupread.alignment.query_name, num_reads, ref_allele, alt_allele, is_indel, is_ref,
                   is_alt)
            if pileupread_alignment_query_name in pileupreads:
                # is the read in the ref bag, alt bag, or neither?
                if not is_indel and is_ref and not is_alt:
                    ref_overlapping_aligned_segment_count += 1
                elif not is_indel and not is_ref and is_alt:
                    alt_overlapping_aligned_segment_count += 1
                else:  # indel based pivots and non-biallelic leftovers
                    pass
            else:
                pileupreads[pileupread_alignment_query_name] = pileupread
                # is the read in the ref bag, alt bag, or neither?
                if not is_indel and is_ref and not is_alt:  # ref supporting
                    ref_supporting_pileupread_alignment_query_name += [pileupread_alignment_query_name]
                    num_ref_supporting_reads += 1
                elif not is_indel and not is_ref and is_alt:  # alt supporting
                    alt_supporting_pileupread_alignment_query_name += [pileupread_alignment_query_name]
                    num_alt_supporting_reads += 1
                else:  # indel based pivots and non-biallelic leftovers
                    pass

        # indexed_pileupreads = \
        #     ArtifactContingencyTableAnalysisUtils.retrieve_indexed_pileupreads(pileupcolumn_knapsack,
        #                                                                        pileupreads.keys(), pileupcolumn_mask)
        #
        # ref_supporting_indexed_pileupreads = OrderedDict()
        # # for
        # # ref supporting counts
        # ref_ref_pileupread_bp_count = \
        #         ArtifactContingencyTableAnalysisUtils.retrieve_ref_bp_count()
        # ref_non_ref_pileupread_bp_count = \
        #     ArtifactContingencyTableAnalysisUtils.retrieve_non_ref_bp_count()
        #
        # # alt supporting counts
        # alt_ref_pileupread_bp_count = \
        #         ArtifactContingencyTableAnalysisUtils.retrieve_ref_bp_count()
        # ref_non_ref_pileupread_bp_count = \
        #     ArtifactContingencyTableAnalysisUtils.retrieve_non_ref_bp_count()

                    # print "%s:%s:%s:%s:%s:%s:%s:%s" % \
                    #       (pileupread.alignment.query_name, chrom, start, ref_allele, alt_allele, is_indel, is_ref,
                    #        is_alt)
                    # continue

                    # if pileupread.alignment.query_name in aligned_segment_names:  # remove overlapping reads
                    #     # is the read in the ref or alt bag or neither?
                    #     if not is_indel and is_ref and not is_alt:
                    #         ref_overlapping_aligned_segment_count += 1
                    #     elif not is_indel and not is_ref and is_alt:
                    #         alt_overlapping_aligned_segment_count += 1
                    #     else:  # indel based pivots and non-biallelic leftovers
                    #         pass
                    # else:
                    #     aligned_segment_names.add(pileupread.alignment.query_name)
                    #     if is_ref and not is_indel and not is_alt:
                    # #         ref_non_ref_pileupread_bp_count += \
                    # #             ArtifactContingencyTableAnalysisUtils.retrieve_non_ref_pileupread_bp_count(pileupread,
                    # #                 chrom, binarized_position_names_mask, ref_alleles, binarize_indel_lengths=True)
                    #         ref_ref_pileupread_bp_count += \
                    #             ArtifactContingencyTableAnalysisUtils.retrieve_ref_pileupread_bp_count(pileupread,
                    #                 chrom, binary_position_names_mask, ref_alleles)
                    #         # TODO: subtract 1 in cases where it is ref supporting

                    #         ref_soft_clipped_pileupread_bp_count += \
                    #             ArtifactContingencyTableAnalysisUtils.retrieve_soft_clipped_pileupread_bp_count(pileupread,
                    #                 binarized_position_names_mask, binarize_soft_clipped_lengths=True)
                    #     elif is_alt and not is_indel and not is_ref:
                    #         alt_non_ref_pileupread_bp_count += \
                    #             ArtifactContingencyTableAnalysisUtils.retrieve_non_ref_pileupread_bp_count(pileupread,
                    #                 chrom, binarized_position_names_mask, ref_alleles, binarize_indel_lengths=True)
                            # TODO: subtract 1 in cases where it is alt supporting
                    #         alt_ref_pileupread_bp_count += \
                    #             ArtifactContingencyTableAnalysisUtils.retrieve_ref_pileupread_bp_count(pileupread,
                    #                 chrom, binarized_position_names_mask, ref_alleles)
                    #         alt_soft_clipped_pileupread_bp_count += \
                    #             ArtifactContingencyTableAnalysisUtils.retrieve_soft_clipped_pileupread_bp_count(pileupread,
                    #                 binarized_position_names_mask, binarize_soft_clipped_lengths=True)
                    #     else:  # multi-allelic sites
                    #         pass
            # break  # only one pass permitted

        return ArtifactContingencyTableAnalysis(alt_non_ref_pileupread_bp_count=alt_non_ref_pileupread_bp_count,
                                                alt_ref_pileupread_bp_count=alt_ref_pileupread_bp_count,
                                                ref_non_ref_pileupread_bp_count=ref_non_ref_pileupread_bp_count,
                                                ref_ref_pileupread_bp_count=ref_ref_pileupread_bp_count,
                                                alt_soft_clipped_pileupread_bp_count=alt_soft_clipped_pileupread_bp_count,
                                                ref_soft_clipped_pileupread_bp_count=ref_soft_clipped_pileupread_bp_count,
                                                alt_overlapping_aligned_segment_count=alt_overlapping_aligned_segment_count,
                                                ref_overlapping_aligned_segment_count=ref_overlapping_aligned_segment_count)
