from BasePairUtils import BasePairUtils
from collections import OrderedDict
import math
import pandas


class ArtifactAnalysisTableUtils():

    EPS = 2.2204460492503131e-16
    BAM_CSOFT_CLIP = 4

    @staticmethod
    def pileupread_has_indel(pileupread):
        return bool(pileupread.is_del | (pileupread.indel > 0))  # del or ins

    @staticmethod
    def pileupread_has_alt_allele(pileupread, alt_allele):
        return bool(~ArtifactAnalysisTableUtils.pileupread_has_indel(pileupread) &
                    (pileupread.alignment.query_sequence[pileupread.query_position] == alt_allele))

    @staticmethod
    def pileupread_has_ref_allele(pileupread, ref_allele):
        # soft clip regions are discarded
        return bool(~ArtifactAnalysisTableUtils.pileupread_has_indel(pileupread) &
                    (pileupread.alignment.query_sequence[pileupread.query_position] == ref_allele))

    @staticmethod
    def _retrieve_non_ref_length(pileupread, ref_allele, binarize_length=True):
        length = 0
        if not pileupread.is_del and pileupread.indel > 0:
            if binarize_length:
                length += 1  # for the subsequent insertion
            else:
                length += int(pileupread.indel)
            if pileupread.alignment.query_sequence[pileupread.query_position] != ref_allele:
                length += 1  # for a SNP (adds 1 for the SNP in front of an insertion)
        elif not pileupread.is_del and pileupread.indel < 0:
            if binarize_length:
                length += 1  # for the subsequent deletion
            else:
                length += int(math.fabs(pileupread.indel))
            if pileupread.alignment.query_sequence[pileupread.query_position] != ref_allele:
                length += 1  # for a SNP (adds 1 for the SNP in front of a deletion)
        elif not pileupread.is_del:
            if pileupread.alignment.query_sequence[pileupread.query_position] != ref_allele:
                length += 1  # for a SNP
        return length

    @staticmethod
    def retrieve_indexed_pileupreads(pileupcolumn_knapsack, pileupread_alignment_query_names, pileupcolumn_mask):
        indexed_pileupreads = OrderedDict()
        for pileupcolumn_name in pileupcolumn_knapsack.pileupcolumn_names:
            if pileupcolumn_name in pileupcolumn_mask and \
                    not pileupcolumn_mask[pileupcolumn_name]:  # pileupcolumn is not masked
                indexed_pileupreads[pileupcolumn_name] = []
                for pileupread_alignment_query_name in pileupread_alignment_query_names:
                    pileupread = pileupcolumn_knapsack.retrieve_pileupread(pileupcolumn_name,
                                                                           pileupread_alignment_query_name)
                    if pileupread:
                        indexed_pileupreads[pileupcolumn_name] += [pileupread]
        return indexed_pileupreads

    @staticmethod
    def retrieve_ref_alleles(pileupcolumn_knapsack):
        ref_alleles = OrderedDict()
        for pileupcolumn_name in pileupcolumn_knapsack.pileupcolumn_names:
            ref_allele = pileupcolumn_knapsack.retrieve_ref_allele(pileupcolumn_name)
            if ref_allele:
                ref_alleles[pileupcolumn_name] = ref_allele
        return ref_alleles

    @staticmethod
    def retrieve_ref_bp_count(indexed_pileupreads, ref_alleles, pileupread_alignment_query_names=None):
        pileupread_alignment_query_names = [] if not pileupread_alignment_query_names else \
            pileupread_alignment_query_names
        count = 0
        for pileupcolumn_name in indexed_pileupreads:  # iterate over columns
            ref_allele = ref_alleles[pileupcolumn_name]
            pileupreads = indexed_pileupreads[pileupcolumn_name]
            for pileupread in pileupreads:  # iterate over rows
                if pileupread.alignment.query_name in pileupread_alignment_query_names and \
                        ArtifactAnalysisTableUtils.pileupread_has_ref_allele(pileupread, ref_allele):
                    count += 1
        return count

    @staticmethod  # does not include soft clipped counts
    def retrieve_non_ref_bp_count(indexed_pileupreads, ref_alleles, pileupread_alignment_query_names=None,
                                  binarize_length=True):
        pileupread_alignment_query_names = [] if not pileupread_alignment_query_names else \
            pileupread_alignment_query_names
        count = 0
        for pileupcolumn_name in indexed_pileupreads:  # iterate over columns
            ref_allele = ref_alleles[pileupcolumn_name]
            pileupreads = indexed_pileupreads[pileupcolumn_name]
            for pileupread in pileupreads:  # iterate over rows
                if pileupread.alignment.query_name in pileupread_alignment_query_names:
                    count += ArtifactAnalysisTableUtils._retrieve_non_ref_length(pileupread, ref_allele,
                                                                                 binarize_length)
        return count

    @staticmethod
    def _retrieve_soft_clipped_length(pileupread, binarize_length=True):
        length = 0
        for cigartuple in pileupread.alignment.cigartuples:
            if cigartuple[0] == ArtifactAnalysisTableUtils.BAM_CSOFT_CLIP:
                if not binarize_length:
                    length += cigartuple[1]
                else:
                    length += 1
        return length

    @staticmethod
    def retrieve_soft_clipped_bp_count(indexed_pileupreads, pileupread_alignment_query_names=None,
                                       binarize_length=True):
        pileupread_alignment_query_names = [] if not pileupread_alignment_query_names else \
            pileupread_alignment_query_names
        count = 0
        for pileupcolumn_name in indexed_pileupreads:  # iterate over columns
            pileupreads = indexed_pileupreads[pileupcolumn_name]
            for pileupread in pileupreads:  # iterate over rows
                if pileupread.alignment.query_name in pileupread_alignment_query_names:
                    count += ArtifactAnalysisTableUtils._retrieve_soft_clipped_length(pileupread, binarize_length)
        return count

    @staticmethod
    def render_contingency_table(dataTable):
        return [[dataTable.alt_non_ref_pileupread_bp_count, dataTable.ref_non_ref_pileupread_bp_count],
                [dataTable.alt_ref_pileupread_bp_count, dataTable.ref_ref_pileupread_bp_count]]

    @staticmethod
    def render_soft_clipped_contingency_table(dataTable):
        return [[dataTable.alt_non_ref_soft_clipped_pileupread_bp_count,
                 dataTable.ref_non_ref_soft_clipped_pileupread_bp_count],
                [dataTable.alt_ref_soft_clipped_pileupread_bp_count,
                 dataTable.ref_ref_soft_clipped_pileupread_bp_count]]

    @staticmethod
    def retrieve_table_as_series(dataTable, prefix=None, suffix=None):
        prefix = "" if prefix is None else prefix
        suffix = "" if suffix is None else suffix

        data = OrderedDict()

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

        return pandas.Series(data=OrderedDict([("%s%s%s" % (prefix, key, suffix), val)
                                               for key, val in data.iteritems()]))