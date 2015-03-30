from collections import OrderedDict


class BasePairUtils(object):


    BAM_CMATCH = 0
    BAM_CINS = 1
    BAM_CDEL = 2
    BAM_CSOFT_CLIP = 4

    @staticmethod
    def retrieve_pileupcolumn_name(chrom, start, end):
        return chrom + ":" + str(start) + ":" + str(end)

    @staticmethod
    def intersect_pileupcolumn_masks(pileupcolumn_mask1, pileupcolumn_mask2):
        mask = OrderedDict()
        for pileupcolumn_name in pileupcolumn_mask1.mask:
            if pileupcolumn_name in pileupcolumn_mask2.mask:
                mask[pileupcolumn_name] = \
                    pileupcolumn_mask1.mask[pileupcolumn_name] | pileupcolumn_mask2.mask[pileupcolumn_name]
        return mask

    @staticmethod
    def determine_soft_clipped_region_length(position, aligned_segment):
        aligned_pair_positions = [aligned_pair[1] for aligned_pair in aligned_segment.aligned_pairs]
        if aligned_pair_positions[0] > position and aligned_segment.cigartuples[0][0] == BasePairUtils.BAM_CSOFT_CLIP:
            return aligned_segment.cigartuples[0][1]
        elif aligned_pair_positions[-1] < position and aligned_segment.cigartuples[0][0] == \
                BasePairUtils.BAM_CSOFT_CLIP:
            return aligned_segment.cigartuples[-1][1]
        return 0

    @staticmethod
    def aligned_segment_has_ins(aligned_pairs, position):

        pass

    @staticmethod
    def aligned_segment_has_del(aligned_segment_positions, position):
        return position not in aligned_segment_positions