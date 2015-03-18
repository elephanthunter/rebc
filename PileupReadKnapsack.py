# from BasePairUtils import BasePairUtils
#
#
# class PileupReadKnapsack():
#
#     def __init__(self):
#         pass
#
#     @classmethod
#     def create(cls, chrom, start, end, sample_bam_file, ref_seq_file):
#         #ref_allele = ref_seq_file.fetch(chrom, start=start, end=end)  # [start,end) region is called
#         for pileupcolumn in sample_bam_file.pileup(chrom, start, end, truncate=False):
#             if pileupcolumn.pos != start:
#                 base_pair_read_knapsack = BasePairReadKnapsack.create(chrom=chrom, start=pileupcolumn.pos,
#                                                                       end=pileupcolumn.pos+1, pileupcolumn=pileupcolumn,
#                                                                       refseqfile=ref_seq_file)
#                 base_pair_read_knapsacks[base_pair_read_knapsack.name] = base_pair_read_knapsack
#                 # aggregate_counts = base_pair_read_knapsack.retrieve_base_pair_aggregate_counts()
#
#                 # TODO: Using aggregate counts, determine whether the base pair should be masked or not
#                 # For now, set all of them to True
#
#                 mask[base_pair_read_knapsack.name] = False
#         pass