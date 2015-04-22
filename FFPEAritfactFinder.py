import argparse
from pysam import AlignmentFile
from pysam import FastaFile
import pandas
import sys
import time
import numpy
from PileupColumnMask import PileupColumnMask
from PileupColumnUtils import PileupColumnUtils
from ArtifactAnalysisTable import ArtifactAnalysisTable
from ArtifactAnalysisTable import ArtifactAnalysisTableUtils
from PileupColumnKnapsack import PileupColumnKnapsack
from SomaticMutation import SomaticMutation
from DataType import DataType
import math
import collections


# Q. Does reverse strand impact query sequence?
# A. No, it does not. Pysam orders them correctly.

# Q. Two reasons for the clipping?
# A. 1. Base quality went down; 2. There were alternates.

CODING_VARIANT_CLASSIFICATION = ["Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Silent", "Splice_Site",
                                 "In_Frame_Ins", "In_Frame_Del", "Nonsense_Mutation", "Start_Codon_Del"]
VALID_CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
                     "19", "20", "21", "22", "X", "Y", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                     "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                     "chr20", "chr21", "chr22", "chrX", "chrY"]


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
    somatic_mutations_dataframe = pandas.read_csv(args.input_maf_filename, sep="\t", header=0, comment="#",
                                                  usecols=["Chromosome", "Start_position", "End_position",
                                                           "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type",
                                                           "Variant_Classification", "Tumor_Sample_Barcode",
                                                           "Matched_Norm_Sample_Barcode"],
                                                  dtype={"Chromosome": numpy.str, "Start_position": numpy.int,
                                                         "End_position": numpy.int, "Reference_Allele": numpy.str,
                                                         "Tumor_Seq_Allele2": numpy.str, "Variant_Type": numpy.str,
                                                         "Tumor_Sample_Barcode": numpy.str,
                                                         "Matched_Norm_Sample_Barcode": numpy.str,
                                                         "Variant_Classification": numpy.str}, low_memory=False)

    criteria = (somatic_mutations_dataframe["Variant_Type"].isin(["SNP"])) & \
               (somatic_mutations_dataframe["Variant_Classification"].isin(CODING_VARIANT_CLASSIFICATION)) & \
               (somatic_mutations_dataframe["Chromosome"].isin(VALID_CHROMOSOMES))
    somatic_mutations_dataframe = somatic_mutations_dataframe[criteria]
    somatic_mutations_dataframe = somatic_mutations_dataframe.sort(["Tumor_Sample_Barcode",
                                                                    "Matched_Norm_Sample_Barcode"])

    mutational_features_dataframe = None
    if (not args.case_sample_bam_filename or not args.control_sample_bam_filename) \
            and args.sample_bam_filename is not None:
        samples = pandas.read_csv(args.sample_bam_filename, sep="\t")
        samples = samples[["sample_id", "clean_bam_file_capture"]]  # subset
        somatic_mutations_dataframe = \
            somatic_mutations_dataframe[(somatic_mutations_dataframe["Tumor_Sample_Barcode"].isin(samples["sample_id"]))
                                        & (somatic_mutations_dataframe["Matched_Norm_Sample_Barcode"].isin(samples["sample_id"]))]

        somatic_mutations_dataframe_grouped = somatic_mutations_dataframe.groupby(["Tumor_Sample_Barcode",
                                                                                   "Matched_Norm_Sample_Barcode"])
        column_names = None
        for _, sample_somatic_mutations_dataframe in somatic_mutations_dataframe_grouped:
            case_sample_bam_filename = \
                samples[samples["sample_id"].isin([sample_somatic_mutations_dataframe["Tumor_Sample_Barcode"].iloc[0]])].iloc[0, 1]
            control_sample_bam_filename = \
                samples[samples["sample_id"].isin([sample_somatic_mutations_dataframe["Matched_Norm_Sample_Barcode"].iloc[0]])].iloc[0, 1]
            sample_mutational_features_dataframe = \
                retrieve_mutations_features_as_dataframe(somatic_mutations_dataframe=sample_somatic_mutations_dataframe,
                                                         case_sample_bam_filename=case_sample_bam_filename,
                                                         control_sample_bam_filename=control_sample_bam_filename,
                                                         ref_seq_filename=args.ref_seq_filename)
            if column_names is None:
                column_names = sample_mutational_features_dataframe.columns.tolist()
            if mutational_features_dataframe is None:
                mutational_features_dataframe = pandas.DataFrame(columns=column_names)

            mutational_features_dataframe = pandas.concat([mutational_features_dataframe,
                                                           sample_mutational_features_dataframe])
        mutational_features_dataframe = mutational_features_dataframe[column_names]
    else:
        sample_mutational_features_dataframe = \
            retrieve_mutations_features_as_dataframe(somatic_mutations_dataframe, args.case_sample_bam_filename,
                                                     args.control_sample_bam_filename, args.ref_seq_filename)
        mutational_features_dataframe = sample_mutational_features_dataframe

    if mutational_features_dataframe is not None:
        mutational_features_dataframe.to_csv(args.output_maf_filename, sep="\t", index=False)


#@profile
def retrieve_mutations_features_as_dataframe(somatic_mutations_dataframe, case_sample_bam_filename,
                                             control_sample_bam_filename, ref_seq_filename):
    somatic_mutations = _create_somatic_mutations(somatic_mutations_dataframe)
    case_sample_bam_file = AlignmentFile(case_sample_bam_filename, "rb")
    control_sample_bam_file = AlignmentFile(control_sample_bam_filename, "rb")
    ref_seq_file = FastaFile(ref_seq_filename)

    mutations_features_data = collections.OrderedDict()
    for index in range(0, len(somatic_mutations)):
        t0 = time.time()  # measure start time
        somatic_mutation = somatic_mutations[index]

        # Gather data for cases
        nearby_somatic_mutations = _retrieve_nearby_somatic_mutations(somatic_mutation, somatic_mutations)
        case_pileupcolumn_knapsack = PileupColumnKnapsack.create(somatic_mutation, case_sample_bam_file, ref_seq_file)

        case_pileupcolumn_mask = PileupColumnMask.create(case_pileupcolumn_knapsack, DataType.somatic,
                                                         nearby_somatic_mutations)

        # Gather data for controls
        control_pileupcolumn_knapsack = PileupColumnKnapsack.create(somatic_mutation, control_sample_bam_file,
                                                                    ref_seq_file)
        control_pileupcolumn_mask = PileupColumnMask.create(control_pileupcolumn_knapsack, DataType.germline)

        # Determine what positions to mask
        pileupcolumn_names_mask = PileupColumnUtils.intersect_pileupcolumn_masks(case_pileupcolumn_mask,
                                                                                 control_pileupcolumn_mask)
        if len(pileupcolumn_names_mask) != 0:
            case_data_series = None
            control_data_series = None
            for pileupcolumn in case_sample_bam_file.pileup(somatic_mutation.chrom, somatic_mutation.start,
                                                            somatic_mutation.end, truncate=True):
                case_data_table = ArtifactAnalysisTable.create(somatic_mutation.ref_allele, somatic_mutation.alt_allele,
                                                               pileupcolumn, case_pileupcolumn_knapsack,
                                                               pileupcolumn_names_mask)
                case_data_series = ArtifactAnalysisTableUtils.retrieve_table_as_series(case_data_table, prefix="case_")
                break

            for pileupcolumn in control_sample_bam_file.pileup(somatic_mutation.chrom, somatic_mutation.start,
                                                               somatic_mutation.end, truncate=True):
                control_data_table = ArtifactAnalysisTable.create(somatic_mutation.ref_allele,
                                                                  somatic_mutation.alt_allele, pileupcolumn,
                                                                  control_pileupcolumn_knapsack,
                                                                  pileupcolumn_names_mask)
                control_data_series = ArtifactAnalysisTableUtils.retrieve_table_as_series(control_data_table,
                                                                                          prefix="control_")
                break

            if case_data_series is not None and control_data_series is not None and not case_data_series.empty \
                    and not control_data_series.empty:
                mutation_features = pandas.concat([somatic_mutation.retrieve_as_series(), control_data_series,
                                                   case_data_series])
                if len(mutations_features_data) == 0:
                    for index in mutation_features.index.tolist():
                        mutations_features_data[index] = [mutation_features[index]]
                else:
                    for index in mutation_features.index.tolist():
                        mutations_features_data[index] += [mutation_features[index]]

        t1 = time.time()  # measure end time

        chrom = somatic_mutation.chrom
        start = somatic_mutation.start
        end = somatic_mutation.end
        ref_allele = somatic_mutation.ref_allele
        alt_allele = somatic_mutation.alt_allele
        sys.stdout.write("chrom:%s, start:%s, end:%s, ref:%s, alt:%s, time:%s\n" % (chrom, start, end, ref_allele,
                                                                                    alt_allele, t1-t0))

    case_sample_bam_file.close()
    control_sample_bam_file.close()
    ref_seq_file.close()

    return pandas.DataFrame(data=mutations_features_data)


def _create_somatic_mutations(somatic_mutations_dataframe):
    somatic_mutations_dataframe = somatic_mutations_dataframe.reset_index(drop=True)
    somatic_mutations_dataframe = somatic_mutations_dataframe[["Chromosome", "Start_position", "End_position",
                                                               "Reference_Allele", "Tumor_Seq_Allele2"]]
    return somatic_mutations_dataframe.apply(SomaticMutation.create, axis=1)


def _retrieve_nearby_somatic_mutations(somatic_mutation, somatic_mutations, pileupread_length=75):
    nearby_somatic_mutations = []
    for index in range(0, len(somatic_mutations)):
        if somatic_mutation != somatic_mutations[index] \
                and int(math.fabs(somatic_mutations[index].start - somatic_mutation.start)) < pileupread_length:
            nearby_somatic_mutations += [somatic_mutations[index]]
    return nearby_somatic_mutations


if __name__ == "__main__":
    main()