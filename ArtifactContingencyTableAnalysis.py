import math
import pandas
from ArtifactContingencyTableAnalysisUtils import ArtifactContingencyTableAnalysisUtils


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
    def create(cls, chrom, start, end, ref_allele, alt_allele, bam_file, binarized_position_names_mask=None,
               ref_seq_file=None):

        alt_non_ref_pileupread_bp_count = 0
        alt_ref_pileupread_bp_count = 0

        ref_non_ref_pileupread_bp_count = 0
        ref_ref_pileupread_bp_count = 0

        alt_soft_clipped_pileupread_bp_count = 0
        ref_soft_clipped_pileupread_bp_count = 0

        ref_overlapping_aligned_segment_count = 0
        alt_overlapping_aligned_segment_count = 0

        ref_alleles = ArtifactContingencyTableAnalysisUtils.retrieve_ref_alleles(chrom=chrom, start=start, end=end,
                                                                                 refseqfile=ref_seq_file)
        for pileupcolumn in bam_file.pileup(chrom, start, end, truncate=True):
            if pileupcolumn.pos == start:
                aligned_segment_names = set()

                for pileupread in pileupcolumn.pileups:
                    is_indel = ArtifactContingencyTableAnalysisUtils.pileupread_has_indel(pileupread)  # not used
                    is_ref = ArtifactContingencyTableAnalysisUtils.pileupread_has_ref_allele(pileupread, ref_allele)
                    is_alt = ArtifactContingencyTableAnalysisUtils.pileupread_has_alt_allele(pileupread, alt_allele)

                    # print "%s:%s:%s:%s:%s:%s:%s:%s" % \
                    #       (pileupread.alignment.query_name, chrom, start, ref_allele, alt_allele, is_indel, is_ref,
                    #        is_alt)
                    # continue

                    if pileupread.alignment.query_name in aligned_segment_names:  # remove overlapping reads
                        # is the read in the ref or alt bag or neither?
                        if not is_indel and is_ref and not is_alt:
                            ref_overlapping_aligned_segment_count += 1
                        elif not is_indel and not is_ref and is_alt:
                            alt_overlapping_aligned_segment_count += 1
                        else:  # indel based pivots and non-biallelic leftovers
                            pass
                    else:
                        aligned_segment_names.add(pileupread.alignment.query_name)
                        if is_ref and not is_indel and not is_alt:
                    #         ref_non_ref_pileupread_bp_count += \
                    #             ArtifactContingencyTableAnalysisUtils.retrieve_non_ref_pileupread_bp_count(pileupread,
                    #                 chrom, binarized_position_names_mask, ref_alleles, binarize_indel_lengths=True)
                            ref_ref_pileupread_bp_count += \
                                ArtifactContingencyTableAnalysisUtils.retrieve_ref_pileupread_bp_count(pileupread,
                                    chrom, binarized_position_names_mask, ref_alleles)
                            # TODO: subtract 1 in cases where it is ref supporting

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
            break  # only one pass permitted

        return ArtifactContingencyTableAnalysis(alt_non_ref_pileupread_bp_count=alt_non_ref_pileupread_bp_count,
                                                alt_ref_pileupread_bp_count=alt_ref_pileupread_bp_count,
                                                ref_non_ref_pileupread_bp_count=ref_non_ref_pileupread_bp_count,
                                                ref_ref_pileupread_bp_count=ref_ref_pileupread_bp_count,
                                                alt_soft_clipped_pileupread_bp_count=alt_soft_clipped_pileupread_bp_count,
                                                ref_soft_clipped_pileupread_bp_count=ref_soft_clipped_pileupread_bp_count,
                                                alt_overlapping_aligned_segment_count=alt_overlapping_aligned_segment_count,
                                                ref_overlapping_aligned_segment_count=ref_overlapping_aligned_segment_count)
