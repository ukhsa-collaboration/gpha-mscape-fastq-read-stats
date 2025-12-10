#!/usr/bin/env python3
import os
import gzip
import sys
import argparse
import functools
import pyfastx
from functional import seq
from typing import NamedTuple
from collections import Counter

__version__ = "0.1.0"


class QualityStats(NamedTuple):
    """
    Parameters
    ----------
        read_id : str
            ID from the header in the fastq
        mean_phred_score : float
            average of the per-base phred scores for a read,
            representing the average probability that a bases
            are incorrectly called, expressed as a ratio
        read_length : int
            length of the read in bases
        gc_content : float
            ratio of GC bases vs total length
        compression_ratio : float
            ratio of length in bytes of a gzip compressed
            version of a read to the uncompressed version
        paired_base_complexity : float
            ratio of bases that differ from the previous
            base vs total length
    """

    read_id: str
    mean_phred_score: float
    read_length: int
    gc_content: float
    compression_ratio: float
    paired_base_complexity: float


def phred_ascii_to_int(inputPhredAsciiChar):
    """
    Takes in a phred score represented in
    ASCII form. Converts to the ordinal representation
    then applies an offset so we can have a regular
    integer representation.

    See https://en.wikipedia.org/wiki/Phred_quality_score
    """
    return ord(inputPhredAsciiChar) - 33


def phred_to_ratio(inputPhredScore: int) -> float:
    """
    Takes in a phred score (0-40)
    returns the probability of an incorrect base call
    as a ratio (0 to 1)
    """
    return 10 ** ((-1 * inputPhredScore) / 10)

@functools.cache
def phred_ascii_to_ratio(inputPhredAsciiChar):
    return phred_to_ratio(phred_ascii_to_int(inputPhredAsciiChar))

def get_gc_content(inputSequence: str) -> int | float:
    """
    Takes in an uppercased sequence
    Counts the occurences of G and C bases
    Divides that by the length to give a ratio of GC
    to other bases
    """
    seq_counter = Counter(inputSequence)

    return (seq_counter.get("G", 0) + seq_counter.get("C", 0)) / seq_counter.total()


def get_compression_ratio(inputSequence: str) -> float:
    """
    Takes in a sequence and compares the length
    in bytes of a gzip-compressed version of that
    sequence to the uncompressed sequence.

    Returns a ratio.
    """
    inputStr = str(inputSequence)
    binaryInput = bytes(inputStr, "ascii")
    compression_ratio = len(gzip.compress(binaryInput)) / len(binaryInput)

    return compression_ratio


def get_paired_base_complexity(inputSequence: str) -> int | float:
    """
    Count of bases that differ from the previous base
    divided by length of sequence minus 1.

    Gives some unintuitive answers for small sequences, but
    tends towards something sane as the length increases.
    """
    seqlen = len(inputSequence)
    
    ## avoids a ZeroDivisionError
    if seqlen == 1:
        return 0

    return seq(zip(inputSequence, inputSequence[1:])).map(
        lambda x: x[0] != x[1]
    ).sum() / (seqlen - 1)


def get_quality_stats(input_seq_record) -> QualityStats:
    """
    Takes in a pyfastx fastq record object.
    Generates a tuple of quality stats.
    """
    read_id, sequence, quality = input_seq_record

    return_rec = QualityStats(
        read_id=read_id,
        mean_phred_score=(
            seq(list(quality)).map(phred_ascii_to_ratio).average()
        ),
        read_length=len(sequence),
        gc_content=get_gc_content(sequence.upper()),
        compression_ratio=get_compression_ratio(sequence),
        paired_base_complexity=get_paired_base_complexity(sequence),
    )

    return return_rec


def parse_seq_and_get_quality_stats(input_filepath: str) -> None:
    """
    Takes in a path to a fastq.gz file
    Generates quality stats and writes them directly
    to the stdout in tab-delimited format.
    """
    recs = pyfastx.Fastq(input_filepath, build_index=False)
    base_filename = os.path.basename(input_filepath)
    stats = seq(recs).map(get_quality_stats)

    for idx, rec in enumerate(stats):
        ## header column if we're at idx 0
        if not idx:
            sys.stdout.write("\t".join(["input_filename", *rec._fields]))
            sys.stdout.write("\n")
        ## write out records
        sys.stdout.write("\t".join([base_filename, *[str(y) for y in rec]]))
        sys.stdout.write("\n")


def init_argparser():
    parser = argparse.ArgumentParser(
        prog="fastq_read_stats",
        description="Generates basic statistics given a fastq.gz as input",
    )

    parser.add_argument("input", help="Path to an input fastq.gz file")

    return parser


def main():
    args = init_argparser().parse_args()

    parse_seq_and_get_quality_stats(args.input)


if __name__ == "__main__":
    main()
