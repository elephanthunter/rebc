from BasePairUtils import BasePairUtils
from collections import OrderedDict


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
        return bool(~ArtifactContingencyTableAnalysisUtils.pileupread_has_indel(pileupread) &
                    (pileupread.alignment.query_sequence[pileupread.query_position] == ref_allele))

    @staticmethod
    def retrieve_soft_clipped_pileupread_bp_count(pileupread, binarized_position_names_mask=None,
                                                  binarize_soft_clipped_lengths=True):
        count = 0
        if not binarized_position_names_mask:
            for cigartuple in pileupread.alignment.cigartuples:
                if cigartuple[0] == BasePairUtils.BAM_CSOFT_CLIP:  # add soft clipped values
                    if not binarize_soft_clipped_lengths:
                        count += cigartuple[1]
                    else:
                        count += 1
        else:
            aligned_pair_positions = pileupread.alignment.get_reference_positions(full_length=True)
        return count


    # NOTE: edit distance is 0 in cases where insertion or deletion related collapse

    @staticmethod
    def retrieve_non_ref_pileupread_bp_count(pileupread, chrom, binarized_position_names_mask=None,
                                             ref_alleles=None, binarize_indel_lengths=True):
        count = 0
        if not binarized_position_names_mask or not ref_alleles:  # analysis without using a mask
            try:  # number of non-ref base pairs can be computed directly from edit distance
                count = pileupread.alignment.get_tag("NM")  # edit distance (excludes soft clip length but includes
                # ins/del lengths)
            except ValueError:
                pass
            # remove the lengths of the ins/del; count them as 1 only
            if binarize_indel_lengths:
                for cigartuple in pileupread.alignment.cigartuples:
                    if cigartuple[0] == BasePairUtils.BAM_CINS:  # remove ins and add binary 1
                        count -= cigartuple[1] - 1
                    elif cigartuple[0] == BasePairUtils.BAM_CDEL:  # remove del and add binary 1
                        count -= cigartuple[1] - 1

        elif len(binarized_position_names_mask) > 0 and len(ref_alleles) > 0:  # analysis with a mask
            is_prev_aligned_pair_position_masked = False
            is_prev_del_aligned_pair_position_masked = False

            # iterate over aligned base pairs and determine how many are non-ref, ins and dels
            aligned_pair_index = 0
            for aligned_pair in pileupread.aligned_pairs:  # soft clips are removed, ins positions appear as "None"
                aligned_pair_position = aligned_pair[1]  # retrieve position ("None" for ins)
                if aligned_pair[0] is not None and aligned_pair_position is not None:  # position either an alt or a ref
                    aligned_pair_position_name = BasePairUtils.retrieve_pileupcolumn_name(chrom, aligned_pair_position,
                                                                                      aligned_pair_position+1)
                    is_prev_del_aligned_pair_position_masked = False
                    # as these are SNPs, we do not need to check for binarization
                    if aligned_pair_position_name in binarized_position_names_mask and \
                            aligned_pair_position_name in ref_alleles:
                        if not binarized_position_names_mask[aligned_pair_position_name]:  # position was not masked
                            is_prev_aligned_pair_position_masked = False
                            ref_allele = ref_alleles[aligned_pair_position_name]
                            if ref_allele != pileupread.alignment.query_alignment_sequence[aligned_pair_index]:
                                count += 1
                        else:  # position was masked
                            is_prev_aligned_pair_position_masked = True
                elif aligned_pair_position is None:  # position is an ins (not aligned to ref)
                    is_prev_del_aligned_pair_position_masked = False
                    if not binarize_indel_lengths:
                        if not is_prev_aligned_pair_position_masked:
                            count += 1
                    else:
                        if not is_prev_aligned_pair_position_masked:
                            is_prev_aligned_pair_position_masked = True
                            count += 1
                elif aligned_pair[0] is None:  # position is a del (aligned to ref but not within read)
                    aligned_pair_position_name = BasePairUtils.retrieve_pileupcolumn_name(chrom, aligned_pair_position,
                                                                                      aligned_pair_position+1)
                    if aligned_pair_position_name in binarized_position_names_mask and \
                            aligned_pair_position_name in ref_alleles:
                        if not binarized_position_names_mask[aligned_pair_position_name] \
                                and not is_prev_del_aligned_pair_position_masked:  # position was not masked
                            is_prev_del_aligned_pair_position_masked = False
                            if not binarize_indel_lengths:
                                count += 1
                                is_prev_del_aligned_pair_position_masked = True
                            else:
                                count += 1
                        else:  # position was masked
                            is_prev_del_aligned_pair_position_masked = True
                aligned_pair_index += 1
        return count

    @staticmethod
    def retrieve_ref_pileupread_bp_count(pileupread, chrom, position_names_mask=None, ref_alleles=None):
        count = 0
        if not position_names_mask or not ref_alleles:
            # theoretically, the lower bound of the code below can be negative; however, in practice, we will not see
            # this as the read will fail to align
            try:  # edit distance (excludes soft clip length, includes ins/del lengths)
                count = -pileupread.alignment.get_tag("NM")
            except ValueError:
                pass
            for cigartuple in pileupread.alignment.cigartuples:
                if cigartuple[0] == BasePairUtils.BAM_CMATCH:  # add match counts
                    count += cigartuple[1]
                elif cigartuple[0] == BasePairUtils.BAM_CINS:  # negate insertion counts
                    count += cigartuple[1]
                elif cigartuple[0] == BasePairUtils.BAM_CDEL:  # negate deletion counts
                    count += cigartuple[1]
        elif len(position_names_mask) > 0 and len(ref_alleles) > 0:  # remove counts for masked positions
            l = 0
            # iterate over aligned base pairs and determine how many are non-ref, ins and dels
            for aligned_pair_index in xrange(len(pileupread.alignment.aligned_pairs)):  # soft clips are removed, ins positions appear as "None"
                aligned_pair = pileupread.alignment.aligned_pairs[aligned_pair_index]
                aligned_pair_position = aligned_pair[1]  # retrieve position ("None" for ins)
                if aligned_pair[0] is not None and aligned_pair_position is not None:  # position is either an alt or a ref
                    aligned_pair_position_name = BasePairUtils.retrieve_pileupcolumn_name(chrom, aligned_pair_position,
                                                                                      aligned_pair_position+1)
                    # as these are SNPs, we do not need to check for binarization
                    if aligned_pair_position_name in position_names_mask and \
                            aligned_pair_position_name in ref_alleles:
                        if not position_names_mask[aligned_pair_position_name]:  # position was not masked
                            ref_allele = ref_alleles[aligned_pair_position_name]
                            if ref_allele == pileupread.alignment.query_alignment_sequence[aligned_pair_index]:
                                count += 1  # remove counts for masked positions
                    else:
                        print "%s:%s:%s" % (aligned_pair_position_name,
                                            aligned_pair_position_name in position_names_mask,
                                            aligned_pair_position_name in ref_alleles)
                        l += 0

            if count != 76-l:
                print "----------------------%s-------------------------" % count
        return count

    @staticmethod
    def squeeze_pileupread_knapsack():
        return None


    # walk across the positions
    # for a given position, retrieve pileupreads at that location
    # iterate over those pileupreads and determine whether they are part of the the pivot pileup
    # if the read belongs to that set, get the number of alts for it
    @staticmethod
    def retrieve_ref_bp_count(pileupread_knapsack, binary_pileupcolumn_names_mask):

        for pileupcolumn_name in binary_pileupcolumn_names_mask:
            if binary_pileupcolumn_names_mask[pileupcolumn_name]:  # position is masked
                continue
            pileupread_knapsack = pileupread_knapsack[pileupcolumn_name]

            pass

        return -1


    @staticmethod
    def retrieve_ref_alleles(chrom, start, end, refseqfile):
        ref_alleles = OrderedDict()
        start -= 100
        end += 101
        ref_seq = refseqfile.fetch(reference=chrom, start=start, end=end)
        for offset in range(0, len(ref_seq)):
            position_name = BasePairUtils.retrieve_pileupcolumn_name(chrom=chrom, start=start+offset,
                                                                 end=start+offset+1)
            ref_alleles[position_name] = ref_seq[offset]
        return ref_alleles