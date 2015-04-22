import math
from ArtifactAnalysisTableUtils import ArtifactAnalysisTableUtils


class ArtifactAnalysisTable(object):

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
                                                    ArtifactAnalysisTableUtils.EPS)  # col 1 margin
        self._log_total_ref_col_bp_count = math.log(ref_non_ref_pileupread_bp_count + ref_ref_pileupread_bp_count +
                                                    ArtifactAnalysisTableUtils.EPS)  # col 2 margin
        self._log_total_non_ref_row_bp_count = math.log(alt_non_ref_pileupread_bp_count +
                                                        ref_non_ref_pileupread_bp_count +
                                                        ArtifactAnalysisTableUtils.EPS)  # row 1 margin
        self._log_total_ref_row_bp_count = math.log(alt_ref_pileupread_bp_count + ref_ref_pileupread_bp_count +
                                                    ArtifactAnalysisTableUtils.EPS)  # row 2 margin
        self._log_total_bp_count = math.log(alt_non_ref_pileupread_bp_count + alt_ref_pileupread_bp_count +
                                            ref_non_ref_pileupread_bp_count + ref_ref_pileupread_bp_count +
                                            ArtifactAnalysisTableUtils.EPS)  # total count

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
                     ArtifactAnalysisTableUtils.EPS)  # col 1 margin
        self._log_total_ref_col_soft_clipped_bp_count = \
            math.log(ref_non_ref_soft_clipped_pileupread_bp_count + ref_ref_pileupread_bp_count +
                     ArtifactAnalysisTableUtils.EPS)  # col 2 margin
        self._log_total_non_ref_row_soft_clipped_bp_count = \
            math.log(alt_non_ref_soft_clipped_pileupread_bp_count + ref_non_ref_soft_clipped_pileupread_bp_count +
                     ArtifactAnalysisTableUtils.EPS)  # row 1 margin
        self._log_total_ref_row_soft_clipped_bp_count = \
            math.log(alt_ref_pileupread_bp_count + ref_ref_pileupread_bp_count +
                     ArtifactAnalysisTableUtils.EPS)  # row 2 margin
        self._log_total_soft_clipped_bp_count = \
            math.log(alt_non_ref_soft_clipped_pileupread_bp_count + alt_ref_pileupread_bp_count +
                     ref_non_ref_soft_clipped_pileupread_bp_count + ref_ref_pileupread_bp_count +
                     ArtifactAnalysisTableUtils.EPS)  # total count

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

    @classmethod
    #@profile
    def create(cls, ref_allele, alt_allele, pileupcolumn, pileupcolumn_knapsack, pileupcolumn_mask=None):
        ref_overlapping_aligned_segment_count = 0
        alt_overlapping_aligned_segment_count = 0

        pileupread_alignment_query_names = []
        ref_supporting_pileupread_alignment_query_names = []
        alt_supporting_pileupread_alignment_query_names = []

        ref_alleles = ArtifactAnalysisTableUtils.retrieve_ref_alleles(pileupcolumn_knapsack)
        for pileupread in pileupcolumn.pileups:

            # Quality check: use arbitrary threshold of 90% to decide whether to discard the read or not
            if not ArtifactAnalysisTableUtils.pileupread_passes_quality_control(pileupread):
                continue

            pileupread_alignment_query_name = pileupread.alignment.query_name

            is_indel = ArtifactAnalysisTableUtils.pileupread_has_indel(pileupread)  # not used
            is_ref = ArtifactAnalysisTableUtils.pileupread_has_ref_allele(pileupread, ref_allele)
            is_alt = ArtifactAnalysisTableUtils.pileupread_has_alt_allele(pileupread, alt_allele)

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
            elif not is_indel and not is_ref and is_alt:  # alt supporting read
                pileupread_alignment_query_names += [pileupread_alignment_query_name]
                alt_supporting_pileupread_alignment_query_names += [pileupread_alignment_query_name]
            else:  # indel based pivots and non-biallelic leftovers
                pass  # do nothing

        # TODO: what happens when the read is ref supporting but it's mate is alt supporting?
        indexed_pileupreads = ArtifactAnalysisTableUtils.retrieve_indexed_pileupreads(pileupcolumn_knapsack,
                                                                                      pileupread_alignment_query_names,
                                                                                      pileupcolumn_mask)

        # ref supporting counts
        ref_supporting_pileupread_alignment_query_names = set(ref_supporting_pileupread_alignment_query_names)
        ref_ref_pileupread_bp_count = \
            ArtifactAnalysisTableUtils.retrieve_ref_bp_count(indexed_pileupreads, ref_alleles,
                ref_supporting_pileupread_alignment_query_names)
        ref_non_ref_pileupread_bp_count = \
            ArtifactAnalysisTableUtils.retrieve_non_ref_bp_count(indexed_pileupreads, ref_alleles,
                ref_supporting_pileupread_alignment_query_names)
        ref_soft_clipped_pileupread_bp_count = \
            ArtifactAnalysisTableUtils.retrieve_soft_clipped_bp_count(indexed_pileupreads,
                ref_supporting_pileupread_alignment_query_names)

        # alt supporting counts
        alt_supporting_pileupread_alignment_query_names = set(alt_supporting_pileupread_alignment_query_names)
        alt_ref_pileupread_bp_count = \
            ArtifactAnalysisTableUtils.retrieve_ref_bp_count(indexed_pileupreads, ref_alleles,
                alt_supporting_pileupread_alignment_query_names)
        alt_non_ref_pileupread_bp_count = \
            ArtifactAnalysisTableUtils.retrieve_non_ref_bp_count(indexed_pileupreads, ref_alleles,
                alt_supporting_pileupread_alignment_query_names)
        alt_soft_clipped_pileupread_bp_count = \
            ArtifactAnalysisTableUtils.retrieve_soft_clipped_bp_count(indexed_pileupreads,
                alt_supporting_pileupread_alignment_query_names)

        return ArtifactAnalysisTable(alt_non_ref_pileupread_bp_count=alt_non_ref_pileupread_bp_count,
                                     alt_ref_pileupread_bp_count=alt_ref_pileupread_bp_count,
                                     ref_non_ref_pileupread_bp_count=ref_non_ref_pileupread_bp_count,
                                     ref_ref_pileupread_bp_count=ref_ref_pileupread_bp_count,
                                     alt_soft_clipped_pileupread_bp_count=alt_soft_clipped_pileupread_bp_count,
                                     ref_soft_clipped_pileupread_bp_count=ref_soft_clipped_pileupread_bp_count,
                                     alt_overlapping_aligned_segment_count=alt_overlapping_aligned_segment_count,
                                     ref_overlapping_aligned_segment_count=ref_overlapping_aligned_segment_count)
