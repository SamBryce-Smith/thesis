#!/usr/bin/env python3

import pandas as pd
import os
import sys
import itertools
import argparse

def get_group_id(file_name):
    # Extract the group identifier from the file name
    # The group identifier corresponds to the first two fields separated by underscores
    parts = file_name.split('_')
    if len(parts) >= 2:
        return '_'.join(parts[:2])
    else:
        return None


def group_files_by_identifier(directory):
    # Dictionary to store {group_id: [<file1>, <filen>]}
    grouped_files = {}

    # Iterate over files in the specified directory
    for file_name in os.listdir(directory):
        if file_name.endswith('.bed'):
            file_path = os.path.join(directory, file_name)
            group_id = get_group_id(file_name)
            if group_id:
                # Append the file to the list corresponding to the group_id
                if group_id not in grouped_files:
                    grouped_files[group_id] = []
                grouped_files[group_id].append(file_path)

    return grouped_files


def get_pairwise_combinations(file_list):
    # list itertools.combinations returns a list of tuples of combinations
    # convert to df (1 row per tuple) and name the output columns
    return pd.DataFrame(list(itertools.combinations(file_list, 2)), columns=['file_path1', 'file_path2'])


def extract_sample_name(df, file_path_col, output_col):
    """
    Extract sample names from file paths in a DataFrame column and create a new column.

    Parameters:
    - df (pd.DataFrame): Input DataFrame.
    - file_path_col (str): Name of the column containing file paths.
    - output_col (str): Name of the new column to store the extracted sample names.
    """
    
    df[output_col] = df[file_path_col].apply(lambda x: os.path.basename(x).split('.')[0])

    return df

def main(input_dir, output_prefix):

    # Get the dictionary of grouped file paths {group: [file1,file2...filen]}
    # group inferred from APAeval sample names (first two fields when separate by underscore)
    grouped_files_dict = group_files_by_identifier(input_dir)

    # Convert grouped_files_dict to a DataFrame of group_id | file_path
    grouped_files_df = pd.DataFrame([(group_id, file_path) for group_id, file_list in grouped_files_dict.items() for file_path in file_list],
                                     columns=['group_id', 'file_path'])
    
    # print(grouped_files_df)
    
    # Iterate over each group and compute pairwise combinations
    # Dictionary to store {group_id: df of (file_path1 | file_path2)}
    grouped_combinations_dict = {}

    # Iterate over each group and compute pairwise combinations of file paths
    for group_id, file_list in grouped_files_dict.items():
        combinations_df = get_pairwise_combinations(file_list)
        grouped_combinations_dict[group_id] = combinations_df

    # Convert grouped_combinations_dict to a single DataFrame
    grouped_combinations_df = pd.concat(grouped_combinations_dict.values(), keys=grouped_combinations_dict.keys(),names=["group_id"]).reset_index("group_id").reset_index(drop=True)

    # print(grouped_combinations_df)

    # extract sample names from file name (inferred as everythign before first '.')
    grouped_files_df = extract_sample_name(grouped_files_df, "file_path", "sample_name")
    grouped_files_df = grouped_files_df.loc[:, ["group_id", "sample_name", "file_path"]]
    print(grouped_files_df)

    # create 'combination_name' field to label each combination (sample_name1--sample_name2s)
    grouped_combinations_df = extract_sample_name(grouped_combinations_df, "file_path1", "sample_name1")
    grouped_combinations_df = extract_sample_name(grouped_combinations_df, "file_path2", "sample_name2")
    grouped_combinations_df['comparison_name'] = grouped_combinations_df['sample_name1'] + '--' + grouped_combinations_df['sample_name2']
    grouped_combinations_df.drop(columns=['sample_name1', 'sample_name2'], inplace=True)
    grouped_combinations_df = grouped_combinations_df.loc[:, ["group_id", "comparison_name", "file_path1", "file_path2"]]


    print(grouped_combinations_df)

    print("writing to TSV...")
    grouped_files_df.to_csv(output_prefix + ".sample2path.tsv", sep="\t", index=False, header=True)
    grouped_combinations_df.to_csv(output_prefix + ".pairwise_comparisons.tsv", sep="\t", index=False, header=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process .bed files and extract sample names.')
    parser.add_argument('input_dir', type=str, help='Path to input directory containing .bed files')
    parser.add_argument('output_prefix', type=str, help='Output prefix for saved CSV files')

    args = parser.parse_args()

    main(args.input_dir, args.output_prefix)