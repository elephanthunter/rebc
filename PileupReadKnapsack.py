from BasePairUtils import BasePairUtils
from collections import OrderedDict


class PileupReadKnapsack(object):

    # TODO: What to do about overlapping reads? Because I am using a dictionary, it should not matter.

    def __init__(self, chrom, start, end, ref_allele, pileupreads):
        self._chrom = chrom
        self._start = start
        self._end = end
        self._ref_allele = ref_allele
        self._pileupreads = pileupreads
        self._alignment_sequence_base_descriptors = OrderedDict()

        # Create and insert base pair descriptors from each pileup read
        for pileupread_alignment_query_name in pileupreads:
            pileupread = pileupreads[pileupread_alignment_query_name]
            self._alignment_sequence_base_descriptors[pileupread_alignment_query_name] = \
                PileupReadKnapsack.AlignmentSequenceBaseDescriptor.create(start, ref_allele, pileupread)

    @property
    def pileupread_alignment_query_names(self):
        return self._alignment_sequence_base_descriptors.keys()

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
    def pileupreads(self):
        return self._pileupreads

    @property
    def base_pair_aggregate_counts(self):
        return PileupReadKnapsack.BasePairAggregateCounts.create(self._pileupreads,
                                                                 self._alignment_sequence_base_descriptors)

    @classmethod
    def create(cls, chrom, start, end, ref_allele, pileupcolumn):
        pileupreads = OrderedDict()
        for pileupread in pileupcolumn.pileups:
            pileupreads[pileupread.alignment.query_name] = pileupread
        return PileupReadKnapsack(chrom=chrom, start=start, end=end, ref_allele=ref_allele, pileupreads=pileupreads)

    class AlignmentSequenceBaseDescriptor(object):
        def __init__(self, ref_allele, alt_allele, is_ref, is_alt, is_ins, is_del, is_soft_clipped):
            self._ref_allele = ref_allele
            self._alt_allele = alt_allele
            self._is_alt = is_alt
            self._is_ref = is_ref
            self._is_ins = is_ins
            self._is_del = is_del
            self._is_soft_clipped = is_soft_clipped

        @staticmethod
        def create(position, ref_allele, pileupread):
            alt_allele = None
            is_ref = is_alt = is_ins = is_del = is_soft_clipped = False

            # soft clips appear either at the beginning or at the end of the aligned read
            if pileupread.alignment.aligned_pairs[0][1] > position or \
                    pileupread.alignment.aligned_pairs[-1][1] < position:
                is_soft_clipped = True
                return PileupReadKnapsack.AlignmentSequenceBaseDescriptor(ref_allele=ref_allele, alt_allele=alt_allele,
                                                                          is_ref=is_ref, is_alt=is_alt, is_ins=is_ins,
                                                                          is_del=is_del,
                                                                          is_soft_clipped=is_soft_clipped)

            if pileupread.is_del:  # deletions' positions are NOT included in the aligned read
                is_del = True
                return PileupReadKnapsack.AlignmentSequenceBaseDescriptor(ref_allele=ref_allele, alt_allele=alt_allele,
                                                                          is_ref=is_ref, is_alt=is_alt, is_ins=is_ins,
                                                                          is_del=is_del,
                                                                          is_soft_clipped=is_soft_clipped)

            # insertions appear as "None"'s w.r.t. position (positions around the ins are inserts)
            aligned_pair_positions = [aligned_pair[1] for aligned_pair in pileupread.alignment.aligned_pairs]
            aligned_pair_position_index = aligned_pair_positions.index(position)
            if not (aligned_pair_position_index+1 < len(aligned_pair_positions) and
                    not aligned_pair_positions[aligned_pair_position_index+1]) or \
                    (aligned_pair_position_index-1 > 0 and not aligned_pair_positions[aligned_pair_position_index-1]):
                is_ins = True
                return PileupReadKnapsack.AlignmentSequenceBaseDescriptor(ref_allele=ref_allele, alt_allele=alt_allele,
                                                                          is_ref=is_ref, is_alt=is_alt, is_ins=is_ins,
                                                                          is_del=is_del,
                                                                          is_soft_clipped=is_soft_clipped)

            if pileupread.alignment.query_sequence[pileupread.query_position] == ref_allele:
                is_ref = True
                return PileupReadKnapsack.AlignmentSequenceBaseDescriptor(ref_allele=ref_allele, alt_allele=alt_allele,
                                                                          is_ref=is_ref, is_alt=is_alt, is_ins=is_ins,
                                                                          is_del=is_del,
                                                                          is_soft_clipped=is_soft_clipped)

            alt_allele = pileupread.alignment.query_sequence[pileupread.query_position]
            is_alt = True
            return PileupReadKnapsack.AlignmentSequenceBaseDescriptor(ref_allele=ref_allele, alt_allele=alt_allele,
                                                                      is_ref=is_ref, is_alt=is_alt, is_ins=is_ins,
                                                                      is_del=is_del, is_soft_clipped=is_soft_clipped)

        @property
        def ref_allele(self):
            return self._ref_allele

        @property
        def alt_allele(self):
            return self._alt_allele

        @property
        def is_soft_clipped(self):
            return self._is_soft_clipped

        @property
        def is_alt(self):
            return self._is_alt

        @property
        def is_ref(self):
            return self._is_ref

        @property
        def is_ins(self):
            return self._is_ins

        @property
        def is_del(self):
            return self._is_del

    class BasePairAggregateCounts():
        def __init__(self, num_ref, num_alt, num_ins, num_del, num_soft_clipped):
            self._num_ref = num_ref
            self._num_alt = num_alt
            self._num_ins = num_ins
            self._num_del = num_del
            self._num_soft_clipped = num_soft_clipped

        @property
        def num_ref(self):
            return self._num_ref

        @property
        def num_alt(self):
            return self._num_alt

        @property
        def num_ins(self):
            return self._num_ins

        @property
        def num_del(self):
            return self._num_del

        @property
        def num_soft_clipped(self):
            return self._num_soft_clipped

        @property
        def alt_allelic_fraction(self):
            return -1

        @property
        def num_alt_adenine(self):
            return None

        @property
        def num_alt_thymine(self):
            return None

        @property
        def num_alt_guanine(self):
            return None

        @property
        def num_alt_cytosine(self):
            return None

        @staticmethod
        def create(start, ref_allele, pileupreads, alignment_sequence_base_descriptors):
            num_alt_adenine = num_alt_thymine = num_alt_guanine = num_alt_cytosine = None
            num_ref = num_alt = num_ins = num_del = num_soft_clipped = 0

            soft_clipped_aligned_segment_names = []
            for aligned_segment_name in alignment_sequence_base_descriptors:
                if aligned_segment_name in pileupreads:
                    alignment_sequence_base_descriptor = alignment_sequence_base_descriptors[aligned_segment_name]
                    alt_allele = alignment_sequence_base_descriptor.alt_allele  # None for del, ins and soft clipped

                    if not alt_allele and alignment_sequence_base_descriptor.is_soft_clipped:
                        soft_clipped_aligned_segment_names += [aligned_segment_name]
                        continue

                    if num_alt_adenine is None and num_alt_thymine is None and num_alt_guanine is None and \
                            num_alt_cytosine is None:  # initialize alt allele based on ref allele
                        if ref_allele == "A":
                            num_alt_adenine = None
                            num_alt_thymine = num_alt_guanine = num_alt_cytosine = 0
                        elif ref_allele == "C":
                            num_alt_cytosine = None
                            num_alt_thymine = num_alt_adenine = num_alt_guanine = 0
                        elif ref_allele == "T":
                            num_alt_thymine = None
                            num_alt_cytosine = num_alt_adenine = num_alt_guanine = 0
                        elif ref_allele == "G":
                            num_alt_guanine = None
                            num_alt_cytosine = num_alt_adenine = num_alt_thymine = 0

                    # alt_allele is either A(denine), G(uanine), C(ytosine) or T(thymine)
                    num_alt_adenine += 1 if alt_allele == "A" else 0
                    num_alt_guanine += 1 if alt_allele == "G" else 0
                    num_alt_cytosine += 1 if alt_allele == "C" else 0
                    num_alt_thymine += 1 if alt_allele == "T" else 0

                    # alt_allele is either an insertion or a deletion
                    num_ins += 1 if not alt_allele and alignment_sequence_base_descriptor.is_ins else 0
                    num_del += 1 if not alt_allele and alignment_sequence_base_descriptor.is_del else 0

            if len(soft_clipped_aligned_segment_names) > 0:  # change in read is next to a soft clip
                for soft_clipped_aligned_segment_name in soft_clipped_aligned_segment_names:
                    pileupread = pileupreads[soft_clipped_aligned_segment_name]
                    soft_clipped_region_length = \
                        BasePairUtils.determine_soft_clipped_region_length(position=start,
                                                                           aligned_segment=pileupread.aligned_segment)
                    if soft_clipped_region_length == 1:  # assigned to the dominant transformation
                        pass
                    else:  # assign to dominant ins/del
                        pass
                    # # alt_allele is either A, G, C or T
                    # num_alt_adenine += 1 if alt_allele == "A" else 0
                    # num_alt_guanine += 1 if alt_allele == "G" else 0
                    # num_alt_cytosine += 1 if alt_allele == "C" else 0
                    # num_alt_thymine += 1 if alt_allele == "T" else 0
                    #
                    # # alt_allele is either an insertion or a deletion
                    # num_ins += 1 if not alt_allele and alignment_sequence_base_descriptor.is_ins else 0
                    # num_del += 1 if not alt_allele and alignment_sequence_base_descriptor.is_del else 0

                #BasePairUtils.determine_soft_clipped_region_length(pos)

            # for aligned_segment_name in alignment_sequence_base_descriptors:
            #     alignment_sequence_base_descriptor = alignment_sequence_base_descriptors[aligned_segment_name]
            #     if alignment_sequence_base_descriptor.is_soft_clipped and aligned_segment_name in pileupreads:
            #         pileupread = pileupreads[aligned_segment_name]
            #
            #
            #         num_soft_clipped += 1 if alignment_sequence_base_descriptor.is_soft_clipped else 0

            return PileupReadKnapsack.BasePairAggregateCounts(num_ref=num_ref, num_alt=num_alt, num_ins=num_ins,
                                                                num_del=num_del, num_soft_clipped=num_soft_clipped)