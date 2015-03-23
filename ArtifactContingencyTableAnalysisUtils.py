from BasePairUtils import BasePairUtils
from collections import OrderedDict
import math


class ArtifactContingencyTableAnalysisUtils():

    EPS = 2.2204460492503131e-16

    @staticmethod
    def pileupread_has_indel(pileupread):
        return bool(pileupread.is_del | (pileupread.indel > 0))  # del or ins

    @staticmethod
    def pileupread_has_alt_allele(pileupread, alt_allele):
        return bool(~ArtifactContingencyTableAnalysisUtils.pileupread_has_indel(pileupread) &
                    (pileupread.alignment.query_sequence[pileupread.query_position] == alt_allele))

    @staticmethod
    def pileupread_has_ref_allele(pileupread, ref_allele):
        # soft clip regions are discarded
        return bool(~ArtifactContingencyTableAnalysisUtils.pileupread_has_indel(pileupread) &
                    (pileupread.alignment.query_sequence[pileupread.query_position] == ref_allele))

    @staticmethod
    def _retrieve_non_ref_length(pileupread, ref_allele, binarize_length=True):
        # TODO: investigate cases where the SNP lies at beginning of an insertion (what happens then?)
        length = 0
        if not pileupread.is_del and pileupread.indel > 0:
            if binarize_length:
                length += 1  # for the subsequent insertion
            else:
                length += int(pileupread.indel)
            if pileupread.alignment.query_sequence[pileupread.query_position] != ref_allele:  # for a SNP
                length += 1  # for a SNP
        elif not pileupread.is_del and pileupread.indel < 0:
            if binarize_length:
                length += 1  # for the subsequent deletion
            else:
                length += int(math.fabs(pileupread.indel))
            if pileupread.alignment.query_sequence[pileupread.query_position] != ref_allele:
                length += 1  # for a SNP
        elif not pileupread.is_del:
            if pileupread.alignment.query_sequence[pileupread.query_position] != ref_allele:
                length += 1  # for a SNP
        return length

    @staticmethod
    def retrieve_soft_clipped_pileupread_bp_count(indexed_pileupreads):
        # TODO: not yet implemented
        count = 0
        for indexed_pileupread in indexed_pileupreads:
            count += 1
        return count

    @staticmethod
    def retrieve_indexed_pileupreads(pileupcolumn_knapsack, pileupread_alignment_query_names,
                                     pileupcolumn_mask):
        indexed_pileupreads = OrderedDict()
        for pileupcolumn_name in pileupcolumn_knapsack.pileupcolumn_names:
            if not pileupcolumn_mask[pileupcolumn_name]:  # pileupcolumn is not masked
                indexed_pileupreads[pileupcolumn_name] = []
                for pileupread_alignment_query_name in pileupread_alignment_query_names:
                    pileupread = pileupcolumn_knapsack.retrieve_pileupread(pileupcolumn_name,
                                                                           pileupread_alignment_query_name)
                    if not pileupread:
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
    def retrieve_ref_bp_count(indexed_pileupreads, ref_alleles):
        count = 0
        for index in indexed_pileupreads:  # iterate over columns
            ref_allele = ref_alleles[index]
            pileupreads = indexed_pileupreads[index]
            for pileupread in pileupreads:  # iterate over rows
                if ArtifactContingencyTableAnalysisUtils.pileupread_has_ref_allele(pileupread, ref_allele):
                    count += 1
        return count

    @staticmethod
    def retrieve_non_ref_bp_count(indexed_pileupreads, ref_alleles):  # does not include soft clip counts
        count = 0
        for index in indexed_pileupreads:  # iterate over columns
            ref_allele = ref_alleles[index]
            pileupreads = indexed_pileupreads[index]
            for pileupread in pileupreads:  # iterate over rows
                count += ArtifactContingencyTableAnalysisUtils._retrieve_non_ref_length(pileupread, ref_allele,
                                                                                        binarize_length=True)
                # if the next is deletion and not indexed then don't call length (next iteration)
        return count