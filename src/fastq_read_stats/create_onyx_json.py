#!/usr/bin/env python3

"""Script and associated functions to take unclassified read QC
results and create an onyx analysis table format json file.
"""

# Imports
import argparse
import logging
import os
import sys
from ast import literal_eval
from pathlib import Path
from typing import List

import pandas as pd
import yaml
from onyx import OnyxConfig, OnyxEnv
from onyx_analysis_helper import onyx_analysis_helper_functions as oa

# Set up onyx config
CONFIG = OnyxConfig(
    domain=os.environ[OnyxEnv.DOMAIN],
    token=os.environ[OnyxEnv.TOKEN],
)


# Functions
# Args and logging
def get_args():
    """Get command line arguments"""
    parser = argparse.ArgumentParser(
        prog="Create read QC analysis table",
        description="""Simple script to generate onyx analysis fields from unclassified read QC results
        """,
    )
    parser.add_argument("--climbid", "-i", type=str, required=True, help="Sample ID")
    parser.add_argument(
        "--server",
        "-s",
        type=str,
        required=True,
        choices=["mscape", "synthscape", "devscape"],
        help="Specify server code is being run on",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Folder to save downloaded files to",
    )
    parser.add_argument(
        "--pipeline_info",
        "-p",
        type=str,
        required=True,
        help="Comma separated str in format: 'pipeline name, version, homepage'",
    )
    parser.add_argument(
        "--readQC_result",
        "-r",
        type=str,
        required=True,
        help="Argument stating whether unclassified read QC ran successfully",
    )
    parser.add_argument(
        "--readQC_file",
        "-f",
        type=str,
        required=True,
        nargs="+",
        help="Filepath(s) for read QC results. Provide one file for single-end or two files for paired-end (R1 R2).",
    )
    parser.add_argument(
        "--context",
        "-c",
        dest="orange_box_version",
        type=str,
        required=False,
        help="Orange box version used to determine if should run.",
    )

    return parser.parse_args()


def set_up_logger(stdout_file):
    """Creates logger for component - all logging messages go to stdout
    log file, error messages also go to stderr log. If component runs
    correctly, stderr is empty.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")

    out_handler = logging.FileHandler(stdout_file, mode="a")
    out_handler.setFormatter(formatter)
    logger.addHandler(out_handler)

    return logger


def get_unclassified_results(unclass_data_files: List[Path]) -> dict:
    """Parses one or two result_data file(s) from unclassified read QC to generate summary statistics.

    Arguments:
        unclass_data_files -- list of Path(s) to file(s) containing unclassified read QC results.
                              Provide one file for single-end, two files for paired-end.
    Returns:
        Dictionary with combined summary statistics and per-read summaries:
          - Combined (top-level) keys are the same structure as before (column -> {stat_name: value})
          - 'unclassified_read_count' top-level key is the combined count across files
          - 'read1' (and 'read2' if provided) contain per-file summary dicts and per-file counts
    """

    if not isinstance(unclass_data_files, list):
        unclass_data_files = [unclass_data_files]

    if len(unclass_data_files) == 0:
        raise ValueError("No read QC files provided")

    if len(unclass_data_files) > 2:
        raise ValueError(
            "This script supports at most two input read QC files for paired-end runs"
        )

    # dfs = []
    per_read_results = {}

    for idx, file_path in enumerate(unclass_data_files):
        p = Path(file_path)
        if not p.exists():
            raise FileNotFoundError(f"Read QC file not found: {p}")

        readqc = pd.read_csv(p, sep="\t")
        # dfs.append(readqc)

        stats = readqc.agg(
            {
                "mean_phred_score": ["min", "max", "mean", "std", "median"],
                "read_length": ["min", "max", "mean", "std", "median"],
                "gc_content": ["min", "max", "mean", "std", "median"],
                "compression_ratio": ["min", "max", "mean", "std", "median"],
                "paired_base_complexity": ["min", "max", "mean", "std", "median"],
            }
        ).to_dict()

        # per-file unclassified count
        stats["unclassified_read_count"] = int(readqc.shape[0])
        stats["unclassified_read_qc_table"] = file_path
        label = f"r{idx + 1}"
        per_read_results[label] = stats

    # Build final dictionary: top-level combined stats plus per-read entries
    result_dict = {}
    for k, v in per_read_results.items():
        result_dict[k] = v

    return result_dict


def make_pipeline_dict(pipeline_str) -> dict:
    """Create dict of pipeline info from input pipeline info string.
    Arguments:
        pipeline_str -- Comma separated string of pipeline info in format
                        'pipeline_name,pipeline_version,pipeline_homepage'
    Returns:
        pipeline_dict -- Dictionary containing pipeline information
    """

    pipeline_list = pipeline_str.split(",")
    pipeline_dict = {
        "name": pipeline_list[0],
        "version": pipeline_list[1],
        "homePage": pipeline_list[2],
    }

    return pipeline_dict


def create_analysis_fields(
    record_id: str,
    unclassified_qc_results: dict,
    server: str,
    pipeline_info: dict,
    orange_box_version: str,
) -> tuple[oa.OnyxAnalysis, int]:
    onyx_analysis = oa.OnyxAnalysis()
    onyx_analysis.add_analysis_details(
        analysis_name="unclassified-read-qc",
        analysis_description="This is an analysis to perform read QC for unclassified reads",
    )
    onyx_analysis.pipeline_name = pipeline_info["name"].strip()
    onyx_analysis.pipeline_version = pipeline_info["version"].strip()
    onyx_analysis.pipeline_url = pipeline_info["homePage"].strip()
    onyx_analysis.outputs = "No s3 outputs"  # ty:ignore[invalid-assignment]

    methods_fail = onyx_analysis.add_methods(methods_dict={})

    _, onyx_versions, query_exitcode = oa.get_data_and_versions_from_onyx(
        sample_id=record_id, server=server, fields=["climb_id"]
    )

    if query_exitcode == 0:
        methods_versions_fail = onyx_analysis.add_versions_to_methods(
            tool_versions={"orange_box_version": orange_box_version},
            onyx_versions=onyx_versions,
        )
    else:
        methods_versions_fail = False

    def get_count(value):
        if value is None:
            return 0
        if isinstance(value, str):
            digits = "".join(ch for ch in value if ch.isdigit())
            return int(digits) if digits else 0
        return int(value)

    read1_count = 0
    read2_count = 0
    total_unclassified_reads = 0

    if "unclassified_read_count" in unclassified_qc_results:
        total_unclassified_reads = get_count(
            unclassified_qc_results["unclassified_read_count"]
        )
        headline_result = f"Total unclassified read count: {total_unclassified_reads}"
    else:
        if "r1" in unclassified_qc_results:
            read1_count = get_count(
                unclassified_qc_results["r1"].get("unclassified_read_count", 0)
            )
        if "r2" in unclassified_qc_results:
            read2_count = get_count(
                unclassified_qc_results["r2"].get("unclassified_read_count", 0)
            )

        total_unclassified_reads = read1_count + read2_count

        if read1_count and read2_count and read1_count == read2_count:
            headline_result = (
                f"Total unclassified read count: (r1 = r2 = {read1_count})"
            )
        elif read1_count or read2_count:
            headline_result = (
                f"Total unclassified read count: {total_unclassified_reads} "
                f"(r1 = {read1_count}, r2 = {read2_count})"
            )
        else:
            headline_result = (
                f"Total unclassified read count: {total_unclassified_reads}"
            )

    results_fail = onyx_analysis.add_results(
        top_result=headline_result,
        results_dict=unclassified_qc_results,
    )

    onyx_analysis.add_server_records(sample_id=record_id, server_name=server)
    required_field_fail, attribute_fail = onyx_analysis.check_analysis_object(
        publish_analysis=False
    )

    if any(
        [
            methods_fail,
            methods_versions_fail,
            results_fail,
            required_field_fail,
            attribute_fail,
        ]
    ):
        exitcode = 1
    else:
        exitcode = 0

    return onyx_analysis, exitcode


# Main script
def main():
    "Main function to process unclassified read QC results ready for submission to onyx"
    args = get_args()
    exitcode = 0

    # Set up log file
    log_file = (
        Path(args.output)
        / f"{args.climbid}.unclassified.readQC.analysis_fields.log.txt"
    )
    set_up_logger(log_file)

    orange_box_version = (
        args.orange_box_version if args.orange_box_version else "unknown"
    )

    if args.readQC_result == "True":
        # Get information from unclassified read QC result files
        result_dict = get_unclassified_results(args.readQC_file)

    else:
        result_dict = {"unclassified_read_count": "No data - empty file"}
    # Convert pipeline info to dict
    pipeline_dict = make_pipeline_dict(args.pipeline_info)

    # Create onyx analysis dict
    onyx_analysis, exitcode = create_analysis_fields(
        record_id=args.climbid,
        unclassified_qc_results=result_dict,
        server=args.server,
        pipeline_info=pipeline_dict,
        orange_box_version=orange_box_version,
    )

    # Exit if analysis object not made correctly
    if exitcode == 1:
        logging.error(
            "Invalid attribute in analysis fields submitted, check logs for details"
        )
        return exitcode

    # Store analysis table info in json file

    onyx_json_file = (
        Path(args.output) / f"{args.climbid}.unclassified_readQC.analysis_fields.json"
    )
    result_file = onyx_analysis.write_analysis_to_json(result_file=onyx_json_file)
    logging.info("Onyx analysis fields written to file %s", result_file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
