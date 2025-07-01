#!/usr/bin/env python3
"""
Script to create a simple scatter plot comparing a chosen numerical value
between two FinnGen data files for common entries, with no color-coding by groups.
"""

import pandas as pd
import matplotlib.pyplot as plt
import gzip
import argparse
import numpy as np
from pathlib import Path

def read_finngen_file(filepath, nrows=None):
    """
    Read FinnGen data file (handles both .gz and regular files)
    
    Args:
        filepath (str): Path to the file
        nrows (int, optional): Number of rows to read. If None, read all rows.
    """
    if filepath.endswith('.gz'):
        with gzip.open(filepath, 'rt') as f:
            df = pd.read_csv(f, sep='\t', nrows=nrows)
    else:
        df = pd.read_csv(filepath, sep='\t', nrows=nrows)
    
    return df

def create_merged_comparison_plot(df1, df2, merge_columns, scatter_value_col, label1_ax, label2_ax):
    """
    Create a single scatter plot comparing a chosen numerical column of common entries between two dataframes.
    
    Args:
        df1, df2: DataFrames to compare.
        merge_columns: List of column names to merge the dataframes.
        scatter_value_col (str): The name of the numerical column to plot on X and Y axes.
        label1_ax (str): Label for the X-axis (file 1's value).
        label2_ax (str): Label for the Y-axis (file 2's value).
    """
    # Select and rename the scatter_value_col to avoid column name conflict after merge
    df1_processed = df1[merge_columns + [scatter_value_col]].copy()
    df2_processed = df2[merge_columns + [scatter_value_col]].copy()

    df1_processed = df1_processed.rename(columns={scatter_value_col: f'{scatter_value_col}_file1'})
    df2_processed = df2_processed.rename(columns={scatter_value_col: f'{scatter_value_col}_file2'})

    # Merge the two dataframes on the specified columns
    # Using 'inner' merge to only compare common entries
    merged_df = pd.merge(df1_processed, df2_processed, on=merge_columns, how='inner')
    print(merged_df)
    if merged_df.empty:
        print(f"Warning: No common entries found between the two files based on columns: {merge_columns}. Cannot create plot.")
        return None

    plt.figure(figsize=(10, 10)) # Square plot for easier comparison of values

    # Scatter plot all points with a single color
    plt.scatter(merged_df[f'{scatter_value_col}_file1'], merged_df[f'{scatter_value_col}_file2'], 
                alpha=0.7, s=50, color='blue', label='Data Points') 

    plt.xlabel(label1_ax)
    plt.ylabel(label2_ax)
    plt.title(f'{scatter_value_col} Comparison (Merged on: {", ".join(merge_columns)})')
    
    # Add a diagonal line for perfect agreement
    val1_max = merged_df[f'{scatter_value_col}_file1'].max()
    val1_min = merged_df[f'{scatter_value_col}_file1'].min()
    val2_max = merged_df[f'{scatter_value_col}_file2'].max()
    val2_min = merged_df[f'{scatter_value_col}_file2'].min()

    max_val = max(val1_max, val2_max)
    min_val = min(val1_min, val2_min)
    
    # Extend the line slightly beyond min/max for better visualization, unless range is tiny
    range_val = max_val - min_val
    if range_val > 0:
        buffer = range_val * 0.05 # 5% buffer
        plot_min = min_val - buffer
        plot_max = max_val + buffer
    else: # Handle cases where all values are the same (e.g., test data)
        plot_min = min_val - 1 # Just a small range
        plot_max = max_val + 1

    plt.plot([plot_min, plot_max], [plot_min, plot_max], 'k--', alpha=0.5, label='y=x (Perfect Agreement)')


    plt.grid(True, alpha=0.3)
    plt.gca().set_aspect('equal', adjustable='box') # Ensure scales are equal

    # Legend for data points and y=x line
    plt.legend(loc='upper left') 

    plt.tight_layout()
    return plt

def main():
    parser = argparse.ArgumentParser(description='Create a scatter plot comparing a chosen numerical value between two FinnGen data files for common entries.')
    parser.add_argument('file1', help='First FinnGen data file (can be .gz)')
    parser.add_argument('file2', help='Second FinnGen data file (can be .gz)')
    parser.add_argument('--output-dir', '-o', default='.', help='Output directory for plots')
    parser.add_argument('--merge-on', '-m', nargs='+', default=["PHENO", "ID"], 
                        help='List of columns to merge the two dataframes on. These columns must exist in both files. Use "None" to merge on all common columns except the scatter value column.')
    parser.add_argument('--scatter-value-col', '-v', default='LOG10P',
                        help='The numerical column to use for the scatter plot values (X and Y axes). Default: LOG10P.')
    parser.add_argument('--label1', default=None, help='Label for the X-axis (first file). Defaults to basename of file1.')
    parser.add_argument('--label2', default=None, help='Label for the Y-axis (second file). Defaults to basename of file2.')
    parser.add_argument('--test', '-t', nargs='?', type=int, const=1000, default=None, metavar='N',
                       help='Test mode: only read first N rows from each file (default: 1000 if no N specified)')
    
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Read files
    print(f"Reading {args.file1}...")
    df1 = read_finngen_file(args.file1, nrows=args.test)
    print(f"Loaded {len(df1)} rows from {args.file1}")
    
    print(f"Reading {args.file2}...")
    df2 = read_finngen_file(args.file2, nrows=args.test)
    print(f"Loaded {len(df2)} rows from {args.file2}")

    # --- Pre-checks for necessary columns ---
    # Check if the chosen scatter value column exists
    if args.scatter_value_col not in df1.columns:
        print(f"Error: '{args.scatter_value_col}' column not found in {args.file1}. Please choose an existing column.")
        return
    if args.scatter_value_col not in df2.columns:
        print(f"Error: '{args.scatter_value_col}' column not found in {args.file2}. Please choose an existing column.")
        return
    
    # Determine the actual merge columns after initial checks
    if 'None' in args.merge_on:
        # Find common columns, exclude the chosen scatter_value_col from merge columns
        common_cols = list(set(df1.columns) & set(df2.columns))
        merge_columns = [col for col in common_cols if col != args.scatter_value_col]
        if not merge_columns:
            print(f"Error: No common columns found between files to merge on, besides '{args.scatter_value_col}'.")
            print("Please specify columns to merge using --merge-on, or ensure files have common identifiers.")
            return
        print(f"Merging on all common columns (excluding '{args.scatter_value_col}'): {merge_columns}")
    else:
        merge_columns = args.merge_on
        # Check if all merge_columns exist in both dataframes
        missing_cols_df1 = [col for col in merge_columns if col not in df1.columns]
        missing_cols_df2 = [col for col in merge_columns if col not in df2.columns]
        
        if missing_cols_df1 or missing_cols_df2:
            print(f"Error: The following merge columns are missing:")
            if missing_cols_df1:
                print(f"  In {args.file1}: {missing_cols_df1}")
            if missing_cols_df2:
                print(f"  In {args.file2}: {missing_cols_df2}")
            return

    # Set default labels if not provided
    label1_ax = args.label1 if args.label1 is not None else Path(args.file1).stem
    label2_ax = args.label2 if args.label2 is not None else Path(args.file2).stem
    
    plt_obj = create_merged_comparison_plot(df1, df2, merge_columns, args.scatter_value_col, label1_ax, label2_ax)

    if plt_obj: # Only save if a plot was created
        # Generate an informative output filename
        merge_suffix = f"merged_on_{'_'.join(merge_columns)}" if merge_columns else "merged_on_all_common"
        output_file = output_dir / f"{args.scatter_value_col}_comparison_{Path(args.file1).stem}_vs_{Path(args.file2).stem}_{merge_suffix}.png"
        
        plt_obj.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Comparison plot saved to {output_file}")
        plt_obj.show()

if __name__ == "__main__":
    main()
