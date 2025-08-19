import glob
import gzip
import argparse

from os import path
from multiprocessing import Pool
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from functools import partial

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    prog="ASAP-to-kite, version 2",
    usage='Script to reformat raw sequencing \ndata from CellRanger-ATAC demultiplexing to a format \ncompatible with kite (kallisto|bustools)"',
)
parser.add_argument(
    "-f",
    "--fastqs",
    help="Path of folder created by mkfastq or bcl2fastq; can be comma separated that will be collapsed into one output.",
    required=True,
)
parser.add_argument(
    "--sample",
    "-s",
    help="Prefix of the filenames of FASTQs to select; can be comma separated that will be collapsed into one output",
    required=True,
)
parser.add_argument(
    "--id",
    "-o",
    default="asap2kite",
    help="A unique run id, used to name output.",
    required=True,
)
parser.add_argument(
    "--conjugation",
    "-j",
    default="TotalSeqA",
    choices=["TotalSeqA", "TotalSeqB"],
    help="String specifying antibody conjugation; either TotalSeqA (default) or TotalSeqB",
)
parser.add_argument(
    "--cores",
    "-c",
    default=4,
    help="Number of cores for parallel processing.",
    type=int,
)
parser.add_argument(
    "--nreads",
    "-n",
    default=10_000_000,
    help="Maximum number of reads to process in one iteration. Decrease this if in a low memory environment (e.g. laptop).",
    type=int,
)
parser.add_argument(
    "--no-rc-R2",
    "-r",
    action="store_true",
    default=False,
    help="By default, the reverse complement of R2 (barcode) is performed (when sequencing with, for example, the NextSeq). Throw this flag to keep R2 as is-- no reverse complement (rc).",
)
args = parser.parse_args()


# ----------
# Functions - I/O
# ----------


# Function to verify R1/R2/R3 are present for nominated samples
def verify_sample_from_R1(list_of_R1s):
    verified_R1s = []
    for R1file in list_of_R1s:
        R2file = R1file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
        R3file = R1file.replace("_R1_001.fastq.gz", "_R3_001.fastq.gz")
        if path.exists(R2file) and path.exists(R3file):
            verified_R1s.append((R1file, R2file, R3file))
    return verified_R1s


# identify all sequencing data that should be parsed for conversion
def parse_directories(folder_of_fastqs, sample_name):
    list_folders = folder_of_fastqs.split(",")
    list_samples = sample_name.split(",")

    all_R1s = []

    # Look into all supplied folders for specific files:
    for path_to_folder in list_folders:

        # Look at all of the possible sample names
        for sample_name_one in list_samples:
            # matching_R1s = glob.glob(path_to_folder+"/*" + sample_name_one + "*" + "_R1_001.fastq.gz")
            matching_R1s = glob.glob(
                f"{path_to_folder}/*{sample_name_one}*_R1_001.fastq.gz"
            )
            all_R1s.extend(matching_R1s)
    return verify_sample_from_R1(all_R1s)


# Process through iterator
def batch_iterator(iterator, batch_size):
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


# Reformat read for export
def formatRead(title, sequence, quality):
    # return("@%s\n%s\n+\n%s\n" % (title, sequence, quality))
    return f"@{title}\n{sequence}\n+\n{quality}\n"


def asap_to_kite_v1(trio, *, rc_R2, conjugation):
    listRead1, listRead2, listRead3 = trio

    # Parse aspects of existing read
    title1, sequence1, quality1 = listRead1
    title2, sequence2, quality2 = listRead2
    _, sequence3, quality3 = listRead3

    # process R2
    if rc_R2:
        # Update sequence with reverse complement
        bio_r2 = str(Seq(sequence2).reverse_complement())
        sequence2 = bio_r2

        # update the quality
        quality2 = quality2[::-1]

    # Recombine attributes based on conjugation logic
    if conjugation == "TotalSeqA":
        new_sequence1 = sequence2 + sequence1[0:10]
        new_sequence2 = sequence3

        new_quality1 = quality2 + quality1[0:10]
        new_quality2 = quality3

    elif conjugation == "TotalSeqB":
        new_sequence1 = sequence2 + sequence3[0:10] + sequence3[25:34]
        new_sequence2 = sequence3[10:25]

        new_quality1 = quality2 + quality3[0:10] + quality3[25:34]
        new_quality2 = quality3[10:25]

    # Prepare reads for exporting
    out_fq1 = formatRead(title1, new_sequence1, new_quality1)
    out_fq2 = formatRead(title2, new_sequence2, new_quality2)

    return out_fq1, out_fq2


def execute_file_conversion(R1s, rc_R2, conjugation, out, *, n_cpu, n_reads):

    for r, _, _ in R1s:
        print(r.replace("_R1_001.fastq.gz", ""))

    outfq1file = out + "_R1.fastq.gz"
    outfq2file = out + "_R2.fastq.gz"
    with gzip.open(outfq1file, "wt") as out_f1, gzip.open(outfq2file, "wt") as out_f2:
        for R1file, R2file, R3file in R1s:

            # Read in fastq in chunks the size of the maximum user tolerated number
            it1 = batch_iterator(FastqGeneralIterator(gzip.open(R1file, "rt")), n_reads)
            it2 = batch_iterator(FastqGeneralIterator(gzip.open(R2file, "rt")), n_reads)
            it3 = batch_iterator(FastqGeneralIterator(gzip.open(R3file, "rt")), n_reads)

            for batches in zip(it1, it2, it3, strict=True):
                asap2kite_partial = partial(
                    asap_to_kite_v1, rc_R2=rc_R2, conjugation=conjugation
                )
                # create and configure the process pool
                with Pool(processes=n_cpu) as pool:
                    fq_data = pool.map(asap2kite_partial, zip(*batches))

                    # process and write out
                    fq_data = tuple(zip(*fq_data))
                    out_f1.writelines(fq_data[0])
                    out_f2.writelines(fq_data[1])


# Main loop -- process input reads and write out the processed fastq files
if __name__ == "__main__":
    folder_of_fastqs = args.fastqs
    sample_name = args.sample
    out = args.id
    n_cpu = args.cores
    n_reads = args.nreads
    rc_R2 = not (args.no_rc_R2)

    conjugation = args.conjugation

    print("\nASAP-to-kite, version 2\n")
    print("User specified options: ")
    print(args)
    # Import files
    R1s_for_analysis = parse_directories(folder_of_fastqs, sample_name)
    execute_file_conversion(
        R1s_for_analysis, rc_R2, conjugation, out, n_cpu=n_cpu, n_reads=n_reads
    )
    print("\nDone!\n")
