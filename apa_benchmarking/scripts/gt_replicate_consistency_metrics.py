#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
from functools import reduce
import os
import sys
import argparse


def jaccard_extend(gr1: pr.PyRanges, gr2: pr.PyRanges, extend_by: int = 0):
    '''Calculate jaccard index between two PyRanges with optional extension

    Parameters
    ----------
    gr1 : pr.PyRanges
        _description_
    gr2 : pr.PyRanges
        _description_
    extend_by : int, optional
        _description_, by default 0
    '''

    try:
        jaccard = gr1.extend(extend_by).statistics.jaccard(gr2.extend(extend_by), strandedness="same")

    except AttributeError:
        # try referring to stats (changed with pyranges version)
        jaccard = gr1.extend(extend_by).stats.jaccard(gr2.extend(extend_by), strandedness="same")

    return jaccard

 
def comparisons_jaccard_index(df: pd.DataFrame, gr_dict: dict, extension_window: int = 0, progress: bool = True):
    '''Calculate jaccard index for a dataframe containing samples to compare

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    gr_dict : dict
        dict of {sample_name: pr.PyRanges}
    extension_window : int, optional
        _description_, by default 0
    progress : bool, optional
        _description_, by default True

    Returns
    -------
    _type_
        _description_
    '''

    assert "comparison_name" in df.columns
    assert "sample_name1" in df.columns
    assert "sample_name2" in df.columns

    jaccard_dict = {}
    i = 1
    num_comparisons = len(df)

    for comparison_name, s1, s2 in zip(df["comparison_name"], df["sample_name1"], df["sample_name2"]):
        if progress:
            print(f"Calculating Jaccard index with extension window {extension_window} for comparison number {i} / {num_comparisons}")
        
        gr1 = gr_dict[s1]
        gr2 = gr_dict[s2]

        jaccard_dict[comparison_name] = jaccard_extend(gr1, gr2, extend_by=extension_window)
        i += 1

    result_df = pd.DataFrame.from_dict(jaccard_dict, orient="index", columns=["jaccardindex_" + str(extension_window)])
    result_df = result_df.reset_index().rename(columns={"index": "comparison_name"})
    
    return result_df


def comparisons_jaccard_index_wrapper(df: pd.DataFrame, gr_dict: dict, extension_sizes=[0, 10, 25, 50, 100], progress=True):
    '''Calculate jaccard index for a dataframe containing samples to compare across a series of extension sizes prior to overlap

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    gr_dict : dict
        _description_
    extension_sizes : list, optional
        _description_, by default [0, 10, 25, 50, 100]
    progress : bool, optional
        _description_, by default True

    Returns
    -------
    _type_
        _description_
    '''
    assert all(isinstance(ext, int) for ext in extension_sizes)

    jaccard_calcs = []

    num_exts = len(extension_sizes)
    for i, ext in enumerate(extension_sizes):
        if progress:
            print(f"Computing Jaccard indices for all comparisons with extension window {ext} - {i + 1} / {num_exts}")
        
        jacc = comparisons_jaccard_index(df, gr_dict, extension_window=ext, progress=False)
        jaccard_calcs.append(jacc)

    merged_df = reduce(lambda left, right: pd.merge(left, right, on="comparison_name", how="outer"), jaccard_calcs)
    
    return merged_df



def main(path_sample2path: str,
         path_pairwise: str,
         output_prefix: str,
         window_sizes: list
         ):
    '''_summary_

    Parameters
    ----------
    path_sample2path : str
        _description_
    path_pairwise : st
        _description_
    window_sizes : list, optional
        _description_, by default [0,10,25,50,100]
    '''

    paths = pd.read_csv(path_sample2path, sep="\t")
    pairwise = pd.read_csv(path_pairwise, sep="\t")

    # Extract sample_names from comparison_name
    pairwise[["sample_name1", "sample_name2"]] = pairwise.comparison_name.str.split("--", expand=True)

    # get list of sample names in groups with at least 2 replicates (for comparison)
    sample_names = set(pairwise["sample_name1"]).union(set(pairwise["sample_name2"]))

    print(f"Number of samples provided - {len(sample_names)}")

    # read in BED files
    # dict of sample_name: pyranges
    print("Reading in BED files...")
    sample2bed = {sample_name: pr.read_bed(file_path) for sample_name, file_path in zip(paths["sample_name"], paths["file_path"]) if sample_name in sample_names}

    # Calculate jaccard indices for each provided comparison across range of input window sizes for overlap
    # df of comparison_name and columns for each window_size (jaccardindex_<window_size>)
    print(f"Number of provided pair-wise comparisons - {pairwise['comparison_name'].nunique()}")
    jaccards = comparisons_jaccard_index_wrapper(pairwise, sample2bed, extension_sizes=window_sizes)

    # add group_id to jaccards_df
    jaccards = pairwise[["group_id", "comparison_name"]].merge(jaccards, on="comparison_name")

    # output jaccards to TSV
    print(f"Writing computed jaccard indices to file - {output_prefix + '.jaccard_index.tsv'}")
    jaccards.to_csv(output_prefix + ".jaccard_index.tsv", sep="\t", index=False, header=True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Compute metrics for assessing ground truth PAS replicate consistency")
    parser.add_argument("sample2path", type=str, help="Path to TSV file mapping sample names to paths (e.g. computed by beds_dir_to_samplesheet.py)")
    parser.add_argument("pairwise", type=str, help="Path to TSV file containing pairwise comparisons (e.g. computed by beds_dir_to_samplesheet.py)")
    parser.add_argument("output_prefix", type=str, help="Prefix for output files - <prefix>.jaccard_index.tsv for jaccard indices")
    parser.add_argument("--window-sizes", nargs='+', type=int, default=[0], help="List of extension window sizes (default: [0])")

    # Parse command-line arguments
    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args.sample2path, args.pairwise, args.output_prefix, args.window_sizes)