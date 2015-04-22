from collections import OrderedDict


class PileupColumnUtils(object):
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