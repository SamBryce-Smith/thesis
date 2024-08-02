#!/usr/bin/env python3

import glob
import os
import sys
import argparse

'''
Simple command line script to generate a sample table of minimal QAPA snakemake pipeline run outputs for downstream preprocessing

Assumes you have a directory where subdirectories respond to 'run name' or 'experiment name'

Under each subdirectory, you need a minimum of:
- sample_table CSV file used as input to QAPA_snakemake pipeline
- Output file of qapa quant for all samples in sample table (output by QAPA_snakemake under all_samples.pau_results.txt)
- contrasts CSV file used as input to QAPA_snakemake pipeline

This script will locate those files and produce a csv of:
experiment_name | sample_table | contrasts_table | qapa_quant_tsv
'''

import os
import glob
import csv


# Function to find files in a directory matching a specific pattern
def find_file(directory, pattern):
    files = glob.glob(os.path.join(directory, pattern))
    if files:
        return files[0]  # Assuming only one file matches the pattern
    else:
        return None


def main(input_dir, output_file, sample_tbl_pattern, contrasts_pattern, qapa_pattern):
    # List to hold the results
    results = []

    # Iterate over subdirectories
    for experiment_name in os.listdir(input_dir):
        experiment_dir = os.path.join(input_dir, experiment_name)
        if not os.path.isdir(experiment_dir):
            continue

        # Find files for each type
        sample_table = find_file(experiment_dir, sample_tbl_pattern)
        contrasts_table = find_file(experiment_dir, contrasts_pattern)
        qapa_quant = find_file(experiment_dir, qapa_pattern)

        # Check if all files are found
        for var_name, var, pattern in zip(["sample_table", "contrasts_table", "qapa_quant"],
                                 [sample_table, contrasts_table, qapa_quant],
                                 [sample_tbl_pattern, contrasts_pattern, qapa_pattern]):
            if not var:
                raise Exception(f"Error: Missing expected file {var_name} in experiment '{experiment_name}' using pattern {pattern}")

        # Add to results
        results.append((experiment_name, sample_table, contrasts_table, qapa_quant))

    # Write results to a CSV file
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['experiment_name', 'sample_table', 'contrasts_table', 'qapa_quant'])
        writer.writerows(results)

    print(f"Results written to '{output_file}'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find files in directories.")
    parser.add_argument("input_dir", help="Input base directory containing experiment-wise subdirectories of QAPA_snakemake outputs.")
    parser.add_argument("output_file", help="Output CSV sample table mapping experiment names to file paths.")
    parser.add_argument("--sampletbl-pattern", default='*_sample_tbl.csv', help="Glob pattern for the sample table file. Default: '*_sample_tbl.csv'")
    parser.add_argument("--contrasts-pattern", default='contrasts_*.csv', help="Glob pattern for the contrasts table file. Default: 'contrasts_*.csv'")
    parser.add_argument("--qapa-pattern", default='all_samples.pau_results.txt', help="Filename for the qapa quant file. Default: 'all_samples.pau_results.txt'")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    main(args.input_dir, args.output_file, args.sampletbl_pattern, args.contrasts_pattern, args.qapa_pattern)