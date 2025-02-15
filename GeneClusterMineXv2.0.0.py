#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Leandro de Mattos Pereira
Junior Researcher, CNP Laboratory
Pedro Leao - Team Leader
Date: June 06, 2023 (Updated: October 03, 2024)
Description:
Script to run antiSMASH on multiple .fna or .fasta files within a directory,
generating a result directory for each input file, with support for parallel processing.
Fixes the "could not parse records from GFF3 file" error when using --genefinding-tool glimmerhmm.
"""

import os
import shutil
import argparse
import subprocess
from pathlib import Path
import logging
import sys
import multiprocessing
from functools import partial

# Default output directory and log file
DEFAULT_LOG_FILE = "logs.txt"


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Script to run antiSMASH on multiple .fna or .fasta files in a directory with custom options."
    )

    # Positional arguments
    parser.add_argument(
        "input_dir",
        type=str,
        help="Directory containing .fna or .fasta files for analysis."
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Directory where results will be saved."
    )

    # antiSMASH help option
    parser.add_argument(
        "--antismash-help",
        action="store_true",
        help="Display antiSMASH help and exit."
    )

    # Basic analysis options
    parser.add_argument(
        "-t", "--taxon",
        choices=["bacteria", "fungi"],
        default="bacteria",
        help="Taxonomic classification of the input sequence (default: bacteria)."
    )
    parser.add_argument(
        "-c", "--cpus",
        type=int,
        default=4,  # Default is 4, change as needed
        help="Number of CPUs for parallel use by antiSMASH (default: 4)."
    )
    parser.add_argument(
        "--databases",
        type=str,
        default="/home/mattoslmp/anaconda3/envs/antismash/lib/python3.9/site-packages/antismash/databases",
        help="Root directory of antiSMASH database files."
    )

    # Additional output options
    parser.add_argument(
        "--output-basename",
        type=str,
        help="Base name for output files within the output directory."
    )
    parser.add_argument(
        "--html-title",
        type=str,
        help="Custom title for the HTML output page."
    )
    parser.add_argument(
        "--html-description",
        type=str,
        help="Custom description to add to the HTML output."
    )
    parser.add_argument(
        "--html-start-compact",
        action="store_true",
        help="Use a compact view by default on the overview page."
    )
    parser.add_argument(
        "--html-ncbi-context",
        dest="html_ncbi_context",
        action="store_true",
        help="Show links to the NCBI genomic context for genes."
    )
    parser.add_argument(
        "--no-html-ncbi-context",
        dest="html_ncbi_context",
        action="store_false",
        help="Do not show links to the NCBI genomic context for genes."
    )
    parser.set_defaults(html_ncbi_context=False)

    # Additional analyses
    parser.add_argument(
        "--all",
        action="store_true",
        help="Enable all additional analyses available in antiSMASH."
    )
    parser.add_argument("--fullhmmer", action="store_true", help="Run HMMer analysis on the entire genome using Pfam profiles.")
    parser.add_argument("--cassis", action="store_true", help="Motif-based prediction of SM gene cluster boundaries.")
    parser.add_argument("--clusterhmmer", action="store_true", help="Run HMMer analysis limited to clusters using Pfam profiles.")
    parser.add_argument("--tigrfam", action="store_true", help="Annotate clusters using TIGRFam profiles.")
    parser.add_argument("--asf", action="store_true", help="Run active site finder analysis.")
    parser.add_argument("--cc-mibig", action="store_true", help="Compare identified clusters with the MIBiG database.")
    parser.add_argument("--cb-general", action="store_true", help="Compare identified clusters with a database of predicted clusters from antiSMASH.")
    parser.add_argument("--cb-subclusters", action="store_true", help="Compare identified clusters with known subclusters that synthesize precursors.")
    parser.add_argument("--cb-knownclusters", action="store_true", help="Compare clusters with known clusters from the MIBiG database.")
    parser.add_argument("--pfam2go", action="store_true", help="Map Pfam to Gene Ontology.")
    parser.add_argument("--rre", action="store_true", help="Run RREFinder in precision mode on all RiPP clusters.")
    parser.add_argument("--smcog-trees", action="store_true", help="Generate phylogenetic trees for secondary metabolite orthologous groups.")
    parser.add_argument("--tfbs", action="store_true", help="Run the transcription factor binding site finder (TFBS) on all clusters.")
    parser.add_argument("--tta-threshold", type=float, default=0.65, help="Minimum GC content for annotating TTA codons (default: 0.65).")

    # Gene prediction options
    parser.add_argument(
        "--genefinding-tool",
        choices=["glimmerhmm", "prodigal", "prodigal-m", "none", "error"],
        default="error",
        help="Gene prediction tool to be used (default: error)."
    )
    parser.add_argument(
        "--genefinding-gff3",
        type=str,
        help="Specify a GFF3 file to extract features from."
    )

    # Logging options
    parser.add_argument(
        "--log-file",
        type=str,
        default=DEFAULT_LOG_FILE,
        help=f"Log file to store messages (default: {DEFAULT_LOG_FILE})."
    )

    # Parallel processes option
    parser.add_argument(
        "--parallel-processes",
        type=int,
        default=4,
        help="Number of parallel processes for running antiSMASH on different files (default: 4)."
    )

    args = parser.parse_args()

    # If user requests antiSMASH help, display and exit
    if args.antismash_help:
        help_command = ["antismash", "--help"]
        subprocess.run(help_command)
        sys.exit()

    # If --all is specified, enable all additional analyses
    if args.all:
        args.fullhmmer = True
        args.cassis = True
        args.clusterhmmer = True
        args.tigrfam = True
        args.asf = True
        args.cc_mibig = True
        args.cb_general = True
        args.cb_subclusters = True
        args.cb_knownclusters = True
        args.pfam2go = True
        args.rre = True
        args.smcog_trees = True
        args.tfbs = True

    return args


def setup_logging(log_file):
    """
    Configure logging to record messages with timestamps.
    """
    logging.basicConfig(
        filename=log_file,
        filemode='a',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    # Also add a handler to display messages on the console
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


def process_file(fasta_path, args, lock):
    """
    Process a single FASTA file using antiSMASH.
    """
    try:
        fasta = fasta_path.name
        # Remove file extension to name the result directory
        fasta_stem = fasta_path.stem
        result_dir = args.output_dir / f"Result_{fasta_stem}"
        index_html_path = result_dir / "index.html"

        # Check if the results directory already exists
        if result_dir.exists():
            if index_html_path.exists():
                with lock:
                    logging.info(f"File '{fasta}' was already processed and 'index.html' exists. Skipping.")
                return (fasta, "already processed")
            else:
                with lock:
                    logging.warning(f"File '{fasta}' was partially processed before and 'index.html' is missing. Reprocessing.")
                try:
                    shutil.rmtree(result_dir)
                    with lock:
                        logging.info(f"Directory '{result_dir}' removed for reprocessing.")
                except Exception as e:
                    with lock:
                        logging.error(f"Error removing directory '{result_dir}': {e}")
                    return (fasta, "failure")

        # Create the result directory
        try:
            result_dir.mkdir(parents=True, exist_ok=True)
            with lock:
                logging.info(f"Created result directory: '{result_dir}'.")
        except Exception as e:
            with lock:
                logging.error(f"Error creating directory '{result_dir}': {e}")
            return (fasta, "failure")

        # Adjust taxon if the gene-finding tool is glimmerhmm
        taxon = args.taxon
        if args.genefinding_tool == "glimmerhmm":
            taxon = "fungi"
            with lock:
                logging.info("Gene-finding tool 'glimmerhmm' selected. Setting taxon to 'fungi'.")

        # Build the antiSMASH command
        command = [
            "antismash",
            str(fasta_path),
            "--taxon", taxon,
            "--cpus", str(args.cpus),
            "--databases", args.databases,
            "--output-dir", str(result_dir)
        ]

        # Add output options
        if args.output_basename:
            command += ["--output-basename", args.output_basename]
        if args.html_title:
            command += ["--html-title", args.html_title]
        if args.html_description:
            command += ["--html-description", args.html_description]
        if args.html_start_compact:
            command.append("--html-start-compact")
        if args.html_ncbi_context:
            command.append("--html-ncbi-context")
        else:
            command.append("--no-html-ncbi-context")

        # ----------------------------------------------------------
        # ADJUSTMENT: disable 'cassis' if using glimmerhmm to avoid GFF3 error
        cassis_enabled = args.cassis
        if args.genefinding_tool == "glimmerhmm" and cassis_enabled:
            cassis_enabled = False
            with lock:
                logging.warning("Disabling --cassis due to conflict with 'glimmerhmm' on eukaryotes.")
        # ----------------------------------------------------------

        # Add additional analyses
        if args.fullhmmer:
            command.append("--fullhmmer")
        if cassis_enabled:
            # If gene-finding tool is 'prodigal', there are limitations:
            if args.genefinding_tool == "prodigal":
                with lock:
                    logging.warning("CASSIS disabled because gene-finding tool is 'prodigal'.")
            else:
                command.append("--cassis")
        if args.clusterhmmer:
            command.append("--clusterhmmer")
        if args.tigrfam:
            command.append("--tigrfam")
        if args.asf:
            command.append("--asf")
        if args.cc_mibig:
            command.append("--cc-mibig")
        if args.cb_general:
            command.append("--cb-general")
        if args.cb_subclusters:
            command.append("--cb-subclusters")
        if args.cb_knownclusters:
            command.append("--cb-knownclusters")
        if args.pfam2go:
            command.append("--pfam2go")
        if args.rre:
            command.append("--rre")
        if args.smcog_trees:
            command.append("--smcog-trees")
        if args.tfbs:
            command.append("--tfbs")
        if args.tta_threshold is not None:
            command += ["--tta-threshold", str(args.tta_threshold)]

        # Add gene-finding options
        command += ["--genefinding-tool", args.genefinding_tool]
        if args.genefinding_gff3:
            command += ["--genefinding-gff3", args.genefinding_gff3]

        # Execute the antiSMASH command
        try:
            with lock:
                logging.info(f"Starting processing of file '{fasta}' with antiSMASH.")
            subprocess.run(command, check=True)
            with lock:
                logging.info(f"Successfully processed file '{fasta}'.")
            return (fasta, "success")
        except subprocess.CalledProcessError as e:
            with lock:
                logging.error(f"Error processing file '{fasta}': {e}")
            return (fasta, "failure")
        except Exception as e:
            with lock:
                logging.error(f"Unexpected error processing file '{fasta}': {e}")
            return (fasta, "failure")
    except Exception as e:
        with lock:
            logging.error(f"Unexpected error in processing file '{fasta_path}': {e}")
        return (fasta_path.name, "failure")


def main():
    args = parse_arguments()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    log_file = Path(args.log_file).resolve()
    args.output_dir = output_dir  # Make args.output_dir a Path

    # Set up logging
    setup_logging(log_file)
    logging.info("Starting processing.")

    # Check if input directory exists
    if not input_dir.is_dir():
        logging.error(f"Input directory '{input_dir}' does not exist or is not a directory.")
        sys.exit(1)

    # Create output directory if it doesn't exist
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"Output directory set to '{output_dir}'.")
    except Exception as e:
        logging.error(f"Error creating output directory '{output_dir}': {e}")
        sys.exit(1)

    # Find all .fna and .fasta files in the input directory
    fna_files = list(input_dir.glob("*.fna")) + list(input_dir.glob("*.fasta"))
    if not fna_files:
        logging.warning(f"No .fna or .fasta files found in '{input_dir}'.")
        sys.exit(1)

    logging.info(f"Found {len(fna_files)} files for processing.")

    manager = multiprocessing.Manager()
    lock = manager.Lock()
    results = []

    # Prepare the partial function for multiprocessing
    process_file_partial = partial(process_file, args=args, lock=lock)

    # Set the number of parallel processes
    num_processes = args.parallel_processes

    # Create a pool of processes
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(process_file_partial, fna_files)

    # Summaries
    success_count = sum(1 for _, status in results if status == "success")
    failure_count = sum(1 for _, status in results if status == "failure")
    already_processed_count = sum(1 for _, status in results if status == "already processed")

    # Finalize processing with a summary
    logging.info("Processing complete.")
    logging.info(f"Total genomes processed successfully: {success_count}")
    logging.info(f"Total genomes failed to process: {failure_count}")
    logging.info(f"Total genomes already processed: {already_processed_count}")

    # Optional: print summary to console
    print("\nProcessing Summary:")
    print(f"Total genomes processed successfully: {success_count}")
    print(f"Total genomes failed to process: {failure_count}")
    print(f"Total genomes already processed: {already_processed_count}")


if __name__ == "__main__":
    main()
