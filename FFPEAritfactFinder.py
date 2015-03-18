import scipy
import scipy.stats
import argparse
from pysam import AlignmentFile
from pysam import FastaFile
import pandas
import math
from BasePairMask import BasePairMask
from collections import OrderedDict
from BasePairUtils import BasePairUtils
from ArtifactContingencyTableAnalysis import ArtifactContingencyTableAnalysis
from BasePairReadKnapsack import BasePairReadKnapsack
from PileupReadKnapsack import PileupReadKnapsack

# Q. Does reverse strand impact query sequence?
# A. No, it does not. Pysam orders them correctly.

# Q. Two reasons for the clipping? base quality went down and they were alternate

CODING_VARIANT_CLASSIFICATION = ["Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Silent", "Splice_Site",
                                 "In_Frame_Ins", "In_Frame_Del", "Nonsense_Mutation", "Start_Codon_Del"]


def retrieve_features(chrom, start, end, ref_allele, alt_allele, bam_file, prefix="case_", base_pair_mask=None):

    dataTable = None
    for pileupcolumn in bam_file.pileup(chrom, start, end, truncate=True):
        if pileupcolumn.pos == start:
            dataTable = ArtifactContingencyTable.create(ref_allele=ref_allele, alt_allele=alt_allele, pileupcolumn=pileupcolumn,
                                         base_pair_mask=base_pair_mask)
            break

    contingency_table = dataTable.render_contingency_table(dataTable=dataTable)
    _, two_sided_pvalue = scipy.stats.fisher_exact(contingency_table, alternative="two-sided")
    _, greater_pvalue = scipy.stats.fisher_exact(contingency_table, alternative="greater")

    clipped_contingency_table = dataTable.render_clipped_contingency_table(dataTable=dataTable)
    _, two_sided_clipped_pvalue = scipy.stats.fisher_exact(clipped_contingency_table, alternative="two-sided")
    _, greater_clipped_pvalue = scipy.stats.fisher_exact(clipped_contingency_table, alternative="greater")

    # row = pandas.Series({"Chromosome": chrom, "Start_position": start, "End_position": end,
    #                      "Reference_Allele": ref_allele, "Tumor_Seq_Allele2": alt_allele
    #                      "Tumor_Sample_Barcode": row["Tumor_Sample_Barcode"],
    #                      "Matched_Norm_Sample_Barcode": row["Matched_Norm_Sample_Barcode"]})
    row = ArtifactContingencyTable.retrieve_as_series(dataTable, prefix)

    two_sided_pvalue += ArtifactContingencyTable.EPS
    greater_pvalue += ArtifactContingencyTable.EPS
    two_sided_clipped_pvalue += ArtifactContingencyTable.EPS
    greater_clipped_pvalue += ArtifactContingencyTable.EPS

    row[prefix + "log_two_sided_p_value"] = math.log(two_sided_pvalue, 10)
    row[prefix + "log_greater_p_value"] = math.log(greater_pvalue, 10)
    row[prefix + "clipped_log_two_sided_p_value"] = math.log(two_sided_clipped_pvalue, 10)
    row[prefix + "clipped_log_greater_p_value"] = math.log(greater_clipped_pvalue, 10)

    return row


def main():
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("--input_maf_filename", dest="input_maf_filename", action="store", required=True,
                        help="Input MAF filename.")
    parser.add_argument("--ref_seq_filename", dest="ref_seq_filename", action="store", required=False,
                        help="Reference genome sequence fasta",
                        default="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
    parser.add_argument("--case_sample_bam_filename", dest="case_sample_bam_filename", action="store", required=False,
                        help="Case sample bam filename")
    parser.add_argument("--control_sample_bam_filename", dest="control_sample_bam_filename", action="store",
                        required=False, help="Control sample bam filename")
    parser.add_argument("--sample_bam_filename", dest="sample_bam_filename", action="store",
                        required=False, help="List of samples and associated bam filenames")
    parser.add_argument("--output_maf_filename", dest="output_maf_filename", action="store",
                        required=True, help="Name of the output MAF filename.")

    # TODO: add option for both germline and somatic mask, etc.

    args, _ = parser.parse_known_args()

    mutations_dataframe = pandas.read_csv(args.input_maf_filename, sep="\t", header=0, comment="#")
    mutations_dataframe = mutations_dataframe[(mutations_dataframe["Variant_Type"] == "SNP")]  # use SNPs
    mutations_dataframe = \
        mutations_dataframe[mutations_dataframe["Variant_Classification"].isin(CODING_VARIANT_CLASSIFICATION)]  # use coding regions
    mutations_2_write = pandas.DataFrame()

    # mutations = mutations.drop_duplicates()

    if (not args.case_sample_bam_filename or not args.control_sample_bam_filename) \
            and args.sample_bam_filename is not None:
        samples = pandas.read_csv(args.sample_bam_filename, sep="\t")
        samples = samples[["sample_id", "clean_bam_file_capture"]]  # subset
        mutations_dataframe = mutations_dataframe.sort(["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"])  # sort by case and control sample names
        mutations_dataframe = mutations_dataframe[(mutations_dataframe["Tumor_Sample_Barcode"].isin(samples["sample_id"])) &
                                                  (mutations_dataframe["Matched_Norm_Sample_Barcode"].isin(samples["sample_id"]))]

        mutations_dataframe_grouped = mutations_dataframe.groupby(["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"])
        column_names = None

        for _, mutations_dataframe_group in mutations_dataframe_grouped:
            case_sample_bam_filename = \
                samples[samples["sample_id"].isin([mutations_dataframe_group["Tumor_Sample_Barcode"][0]])].iloc[0, 1]
            control_sample_bam_filename = \
                samples[samples["sample_id"].isin([mutations_dataframe_group["Matched_Norm_Sample_Barcode"][0]])].iloc[0, 1]

            mutational_features = retrieve_mutational_features(mutations_dataframe=mutations_dataframe_group,
                                                               case_sample_bam_filename=case_sample_bam_filename,
                                                               control_sample_bam_filename=control_sample_bam_filename,
                                                               ref_seq_filename=args.ref_seq_filename)

        #     mutation_group = pandas.merge(control_mutations, case_mutations, how="inner")
        #     column_names = mutation_group.columns.tolist() if column_names is None else column_names
        #     mutations_2_write = pandas.concat([mutations_2_write, mutation_group])
        # mutations_2_write = mutations_2_write[column_names]
    else:
        mutational_features = retrieve_mutational_features(mutations_dataframe=mutations_dataframe,
                                                           case_sample_bam_filename=args.case_sample_bam_filename,
                                                           control_sample_bam_filename=args.control_sample_bam_filename,
                                                           ref_seq_filename=args.ref_seq_filename)
    mutations_2_write.to_csv(args.output_maf_filename, sep="\t", index=False)


def retrieve_mutational_features(mutations_dataframe, case_sample_bam_filename, control_sample_bam_filename,
                                 ref_seq_filename):
    case_sample_bam_file = AlignmentFile(case_sample_bam_filename, "rb")
    control_sample_bam_file = AlignmentFile(control_sample_bam_filename, "rb")
    ref_seq_file = FastaFile(ref_seq_filename)

    for index, mutation_row in mutations_dataframe.iterrows():
        chrom = str(mutation_row["Chromosome"])
        start = mutation_row["Start_position"] - 1  # subtracted to account for zero based indexing when using pysam
        end = mutation_row["End_position"]

        ref_allele = mutation_row["Reference_Allele"]
        alt_allele = mutation_row["Tumor_Seq_Allele2"]

        # Gather data for cases
        case_pileupread_knapsack = PileupReadKnapsack.create(chrom, start, end, case_sample_bam_file, ref_seq_file)
        case_base_pair_mask = BasePairMask.create(case_pileupread_knapsack)

        # Gather data for controls
        control_pileupread_knapsack = PileupReadKnapsack.create(chrom, start, end, control_sample_bam_file,
                                                                ref_seq_file)
        control_base_pair_mask = BasePairMask.create(control_pileupread_knapsack)

        # Determine what positions to mask
        position_names_mask = BasePairUtils.intersect_base_pair_masks(case_base_pair_mask, control_base_pair_mask)

        # case_contingency_table = \
        #     ArtifactContingencyTableAnalysis.create(chrom=chrom, start=start, end=end, ref_allele=ref_allele,
        #                                             alt_allele=alt_allele, bam_file=case_sample_bam_file,
        #                                             binarized_position_names_mask=binarized_position_names_mask,
        #                                             ref_seq_file=ref_seq_file)
        # control_contingency_table = \
        #     ArtifactContingencyTableAnalysis.create(chrom=chrom, start=start, end=end, ref_allele=ref_allele,
        #                                             alt_allele=alt_allele, bam_file=control_sample_bam_file,
        #                                             binarized_position_names_mask=binarized_position_names_mask,
        #                                             ref_seq_file=ref_seq_file)


    case_sample_bam_file.close()
    control_sample_bam_file.close()
    ref_seq_file.close()


if __name__ == "__main__":
    main()