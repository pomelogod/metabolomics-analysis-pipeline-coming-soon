# %%
"""
statistic_analysis.py
language: Python3
author: L. Xie <xiel4@vcu.edu>
"""


import json
import math
import pandas as pd
from pathlib import Path
import numpy as np
import re
from typing import List, Dict, Tuple, Optional
from itertools import combinations
import multiprocessing
import warnings
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats

from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns
import os
import copy


# Import necessary functions from enrichment_analysis.py
from enrichment_analysis import (
    generate_enrichment_table,
    create_bar_chart,
    create_dot_plot,
    create_network_view,
    generate_summary_table,
    pathway_db,
    id_db,
    category_styles,
    pathway_categories
)

def preprocess_name(name):
    """Generate variations of a compound name for improved matching."""
    if pd.isna(name):
        return []
    
    name = str(name).lower()
    variations = [name]

    name = re.sub(r'^d-\(\+\)-', 'l-', name)
    name = re.sub(r'^d-', '', name)

    name_no_non_digit_parentheses = re.sub(r'\(([^1-9]*)\)', '', name)
    name_digits_parentheses_removed = re.sub(r'\((\d+)\)', r'\1', name_no_non_digit_parentheses)
    name = name_digits_parentheses_removed

    replacements = {
        'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'δ': 'delta',
        'ω': 'omega', 'ε': 'epsilon', 'κ': 'kappa', 'τ': 'tau',
        'î': 'i', 'ô': 'o', 'û': 'u', 'ê': 'e', 'â': 'a',
        ' ': '', '-': ' ', ',': '', "'": '',
        'dextro': 'd', 'levo': 'l', 'cis': 'z', 'trans': 'e',
        'dl-': 'l-'
    }
    for old, new in replacements.items():
        variations.append(name.replace(old, new))
    
    variations.extend([
        name,
        re.sub(r'\b(acid|ester|salt|hydrate|anhydrous)\b', '', name),
        re.sub(r'[^a-z0-9]', '', name),
    ])
    
    if 'vitamin' in name:
        variations.append(re.sub(r'vitamin\s*', '', name))
        variations.append(re.sub(r'vitamin\s*', 'vit', name))
    
    endings = ['ate', 'ic', 'ous', 'ide', 'ine']
    for ending in endings:
        if name.endswith(ending):
            variations.append(name[:-len(ending)])
    
    return list(set(var.strip() for var in variations if var.strip()))

def clean_formula(formula):
    """Remove spaces from formula strings for consistent matching."""
    if pd.isna(formula):
        return ""
    return re.sub(r'\s+', '', str(formula))

def parse_formula(formula):
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    return {elem: int(count) if count else 1 for elem, count in elements}

def is_formula_match(query_formula, candidate_formula):
    query_counts = parse_formula(query_formula)
    candidate_counts = parse_formula(candidate_formula)
    
    for element, count in query_counts.items():
        if element == 'H':
            continue
        if candidate_counts.get(element) != count:
            return False
    
    query_h_count = query_counts.get('H', 0)
    candidate_h_count = candidate_counts.get('H', 0)
    return abs(query_h_count - candidate_h_count) <= 1

"""
def create_name_to_hmdb_dict(id_db: pd.DataFrame) -> Dict[str, str]:
    name_to_hmdb = {}
    for _, row in id_db.iterrows():
        hmdb = str(row['HMDB'])
        if pd.isna(hmdb) or hmdb == 'nan':
            continue
        for col in [f'NAME{i}' for i in range(1, 282)]:
            if pd.notna(row[col]):
                variations = preprocess_name(row[col])
                for var in variations:
                    if var:
                        if var not in name_to_hmdb or len(var) > len(name_to_hmdb[var]):
                            name_to_hmdb[var] = hmdb
    return name_to_hmdb
"""

def create_formula_name_to_hmdb_dict(id_db):
    formula_name_to_hmdb = {}
    name_to_hmdb = {}
    for _, row in id_db.iterrows():
        formula = clean_formula(row['FORMULA'])
        hmdb = str(row['HMDB'])
        if pd.isna(hmdb) or hmdb == 'nan':
            continue
        if formula and formula != 'nan':
            if formula not in formula_name_to_hmdb:
                formula_name_to_hmdb[formula] = {}
            for col in [f'NAME{i}' for i in range(1, 282)]:
                if pd.notna(row[col]):
                    variations = preprocess_name(row[col])
                    for var in variations:
                        if var:
                            formula_name_to_hmdb[formula][var] = hmdb
        else:
            for col in [f'NAME{i}' for i in range(1, 282)]:
                if pd.notna(row[col]):
                    variations = preprocess_name(row[col])
                    for var in variations:
                        if var:
                            name_to_hmdb[var] = hmdb
    return formula_name_to_hmdb, name_to_hmdb

def create_hmdb_to_kegg_dict(id_db: pd.DataFrame) -> Dict[str, str]:
    hmdb_to_kegg = {}
    for _, row in id_db.iterrows():
        hmdb = str(row['HMDB'])
        kegg = str(row['KEGG'])
        if pd.notna(hmdb) and pd.notna(kegg) and hmdb != 'nan' and kegg != 'nan':
            hmdb_to_kegg[hmdb] = kegg
    return hmdb_to_kegg

"""
def get_hmdb(name: str, name_to_hmdb: Dict[str, str]) -> Optional[str]:
    variations = preprocess_name(name)
    for var in variations:
        if var in name_to_hmdb:
            return name_to_hmdb[var]
    return None
"""
def get_hmdb(name, formula, formula_name_to_hmdb, name_to_hmdb):
    formula = clean_formula(formula)  # Ensure formula is cleaned as string
    adjusted_variations = preprocess_name(name)

    for candidate_formula in formula_name_to_hmdb:
        if is_formula_match(formula, candidate_formula):
            for var in adjusted_variations:
                if var in formula_name_to_hmdb[candidate_formula]:
                    return formula_name_to_hmdb[candidate_formula][var]

    for var in adjusted_variations:
        if var in name_to_hmdb:
            return name_to_hmdb[var]
    
    return None

def get_kegg(hmdb: str, hmdb_to_kegg: Dict[str, str]) -> Optional[str]:
    return hmdb_to_kegg.get(hmdb)

def extract_group_prefixes(df):
    """Extract unique group prefixes from dataframe columns."""
    group_prefixes = set()
    for col in df.columns:
        if col not in ['Name', 'KEGG', 'HMDB']:
            parts = col.split('_')
            if len(parts) == 3:
                group_prefix = '_'.join(parts[:2])
                group_prefixes.add(group_prefix)
            elif len(parts) == 2:
                group_prefix = '_'.join(parts[:1])
                group_prefixes.add(group_prefix)
            elif len(parts) > 3:
                group_prefix = '_'.join(parts[:3])
                group_prefixes.add(group_prefix)
    return group_prefixes

def parse_group(group):
    return group.split('_')

def generate_group_pairs(groups, variables):
    """Generate group pairs keeping all but one variable the same."""
    pairs = []
    num_variables = len(variables)

    if num_variables == 1:
        sorted_groups = sorted(groups, key=lambda x: variables[0].index(x.split('_')[0]))
        return list(combinations(sorted_groups, 2))
    else:
        for var_index in range(num_variables):
            for group1, group2 in combinations(groups, 2):
                parts1 = group1.split('_')
                parts2 = group2.split('_')

                # Check if groups have at least as many parts as variables
                if len(parts1) < num_variables or len(parts2) < num_variables:
                    continue

                try:
                    # Check if only the current variable is different
                    if parts1[var_index] != parts2[var_index] and all(parts1[i] == parts2[i] for i in range(num_variables) if i != var_index):
                        # Ensure correct order based on the variables list
                        if variables[var_index].index(parts1[var_index]) > variables[var_index].index(parts2[var_index]):
                            group1, group2 = group2, group1
                        pairs.append((group1, group2))
                except IndexError:
                    # If IndexError occurs, skip this pair
                    continue

    # Custom sorting function
    def custom_sort(pair):
        g1, g2 = pair
        p1 = g1.split('_')
        p2 = g2.split('_')
        return tuple(variables[i].index(p[i]) if i < len(p) else -1 for p in [p1, p2] for i in range(num_variables))

    return sorted(pairs, key=custom_sort)

def generate_group_cols(df, group_prefixes):
    return {
        prefix: [col for col in df.columns if col.startswith(prefix) and col not in ['Name', 'KEGG', 'HMDB']]
        for prefix in group_prefixes
    }



def perform_anova_tukey_fold_change(row, group_cols, comparisons):
    try:
         # Convert any set objects to lists and handle missing data
        processed_groups = {}
        for group, cols in group_cols.items():
            cols = list(cols) if isinstance(cols, set) else cols
            # Convert to float and remove NaN values
            values = row[cols].astype(float)
            valid_values = values[~np.isnan(values)].values
            if len(valid_values) > 0:  # Only include groups with valid data
                processed_groups[group] = valid_values


        # Skip if less than 2 groups have valid data
        if len(processed_groups) < 2:
            return []

        # Create DataFrame with available data
        values_list = []
        groups_list = []
        for group, values in processed_groups.items():
            values_list.extend(values)
            groups_list.extend([group] * len(values))

        anova_data = pd.DataFrame({
            'value': values_list,
            'group': groups_list
        })

        # Perform ANOVA if there are at least two groups
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                model = ols('value ~ group', data=anova_data).fit()
                anova_result = sm.stats.anova_lm(model, typ=2)
                anova_pvalue = anova_result['PR(>F)']['group']
                tukey_result = pairwise_tukeyhsd(anova_data['value'], anova_data['group'])
                tukey_summary = tukey_result.summary().data[1:]
                tukey_dict = {(row[0], row[1]): float(row[3]) for row in tukey_summary}
            except:
                anova_pvalue = np.nan
                tukey_dict = {}

        # Calculate group means for available data
        group_means = {group: np.mean(values) for group, values in processed_groups.items()}

        results = []
        for group1, group2 in comparisons:
            # Skip comparison if either group is missing
            if group1 not in processed_groups or group2 not in processed_groups:
                results.append({
                    'Group1': group1,
                    'Group2': group2,
                    'Fold_Change': np.nan,
                    'Log2_Fold_Change': np.nan,
                    'ANOVA_P_Value': np.nan,
                    'Tukey_HSD_P_Value': np.nan,
                    'T_test_P_Value': np.nan
                })
                continue

            # Calculate statistics for available data
            tukey_pvalue = tukey_dict.get((group1, group2), tukey_dict.get((group2, group1), np.nan))
            mean1, mean2 = group_means[group1], group_means[group2]

            # Calculate fold change
            if mean1 != 0:
                fold_change = mean2 / mean1
                log2_fold_change = np.log2(fold_change) if fold_change > 0 else np.nan
            else:
                fold_change = np.nan
                log2_fold_change = np.nan

            # Calculate pairwise statistics
            try:
                comp_data = anova_data[anova_data['group'].isin([group1, group2])]
                comp_model = ols('value ~ group', data=comp_data).fit()
                comp_anova_result = sm.stats.anova_lm(comp_model, typ=2)
                individual_anova_pvalue = comp_anova_result['PR(>F)']['group']
            except:
                individual_anova_pvalue = np.nan

            # Perform t-test if possible
            try:
                values_group1 = processed_groups[group1]
                values_group2 = processed_groups[group2]

                # Handle different lengths for paired t-test
                min_length = min(len(values_group1), len(values_group2))
                if min_length > 0:
                    t_test_p_value = stats.ttest_rel(values_group1[:min_length],
                                                     values_group2[:min_length]).pvalue
                else:
                    t_test_p_value = np.nan
            except:
                t_test_p_value = np.nan

            results.append({
                'Group1': group1,
                'Group2': group2,
                'Fold_Change': fold_change,
                'Log2_Fold_Change': log2_fold_change,
                'ANOVA_P_Value': individual_anova_pvalue,
                'Tukey_HSD_P_Value': tukey_pvalue,
                'T_test_P_Value': t_test_p_value
            })

        return results
    except Exception as e:
        print(f"Error processing row: {e}")
        return []



def process_row(args):
    index, row, group_cols, comparisons = args
    name = row['Name']
    results = perform_anova_tukey_fold_change(row, group_cols, comparisons)
    for result in results:
        result['Name'] = name
    return results

def analyze_and_combine_results(filtered_df, group_cols, comparisons):
    print(f"Starting combined analysis for {len(filtered_df)} rows...")

    if filtered_df.empty:
        print("Warning: Input DataFrame is empty. No analysis performed.")
        return pd.DataFrame()  # Return an empty DataFrame

    args_list = [(index, row, group_cols, comparisons) for index, row in filtered_df.iterrows()]

    num_cores = max(1, multiprocessing.cpu_count() - 1)
    chunk_size = max(1, len(filtered_df) // (num_cores * 10))

    # with ProcessPoolExecutor(max_workers=num_cores) as executor:
    #     all_results = list(tqdm(
    #         executor.map(process_row, args_list, chunksize=chunk_size),
    #         total=len(args_list),
    #         desc="Processing rows"
    #     ))
    all_results = []
    for args in tqdm(args_list, desc="Processing rows"):
        all_results.append(process_row(args))
    results = [item for sublist in all_results for item in sublist if item]
    results_df = pd.DataFrame(results)

    print("\nStatistical analysis complete. Combining results...")

    if results_df.empty:
        print("Warning: No results generated. Returning empty DataFrame.")
        return pd.DataFrame()  # Return an empty DataFrame

    filtered_df_relevant = filtered_df
    results_df = results_df.rename(columns={'Name': 'Name'})
    merged_df = results_df.merge(filtered_df_relevant, on='Name')

    pivot_df = merged_df.pivot_table(
        index=['Name'],
        columns=['Group1', 'Group2'],
        values=['Fold_Change', 'Log2_Fold_Change', 'ANOVA_P_Value', 'Tukey_HSD_P_Value', 'T_test_P_Value'],
        aggfunc='first'
    )

    # Flatten column names
    pivot_df.columns = [f'{col[0]}({col[1]} vs {col[2]})' for col in pivot_df.columns]
    pivot_df.reset_index(inplace=True)

    # Get additional columns, including KEGG
    additional_cols = ['Name', 'KEGG', 'HMDB'] + [col for col in filtered_df.columns if
                                                  col not in ['Name', 'KEGG', 'HMDB', 'FORMULA']]

    # Merge additional columns and filter for non-null KEGG
    metadata_df = filtered_df[additional_cols].drop_duplicates()
    pivot_df = metadata_df[metadata_df['KEGG'].notna()].merge(pivot_df, on='Name', how='left')

    print(f"Final number of rows in statistic result after removing duplicates: {len(pivot_df)}")
    print("Analysis and combination complete.")

    return pivot_df


def create_volcano_plot(df, group1, group2, tukey_p_value_threshold, upregulated_fc_threshold, downregulated_fc_threshold, output_dir='.'):
    def get_color(log2_fc, tukey_p_value):
        if tukey_p_value < tukey_p_value_threshold and log2_fc > math.log2(upregulated_fc_threshold):
            return '#FF4136'
        elif tukey_p_value < tukey_p_value_threshold and log2_fc < math.log2(downregulated_fc_threshold):
            return '#2ECC40'
        else:
            return '#AAAAAA'

    comparison_columns_1 = [
        f'Log2_Fold_Change({group1} vs {group2})',
        f'Tukey_HSD_P_Value({group1} vs {group2})',
        'Name'
    ]

    comparison_columns_2 = [
        f'Log2_Fold_Change({group2} vs {group1})',
        f'Tukey_HSD_P_Value({group2} vs {group1})',
        'Name'
    ]

    if all(col in df.columns for col in comparison_columns_1):
        comparison_columns = comparison_columns_1
    elif all(col in df.columns for col in comparison_columns_2):
        comparison_columns = comparison_columns_2
        group1, group2 = group2, group1
    else:
        print(f"No comparison data found for {group1} and {group2}.")
        return

    comparison_df = df[comparison_columns]
    comparison_df.columns = ['Log2_Fold_Change', 'Tukey_HSD_P_Value', 'Name']
    comparison_df.loc[comparison_df['Tukey_HSD_P_Value'] == 0, 'Tukey_HSD_P_Value'] = 0.0001

    colors = comparison_df.apply(lambda row: get_color(row['Log2_Fold_Change'], row['Tukey_HSD_P_Value']), axis=1)

    plt.figure(figsize=(16, 12))
    sns.set_style("whitegrid")

    scatter = plt.scatter(comparison_df['Log2_Fold_Change'],
                          -np.log10(comparison_df['Tukey_HSD_P_Value']),
                          c=colors,
                          alpha=0.7,
                          s=50)

    plt.title(f'Volcano Plot: {group1} vs {group2}', fontsize=20, fontweight='bold')
    plt.xlabel('Log2 Fold Change', fontsize=14)
    plt.ylabel('-Log10 Tukey HSD P Value', fontsize=14)

    plt.axhline(y=-np.log10(tukey_p_value_threshold), color='#FF851B', linestyle='--', linewidth=2, alpha=0.5)
    plt.axvline(x=math.log2(downregulated_fc_threshold), color='#B10DC9', linestyle='--', linewidth=2, alpha=0.5)
    plt.axvline(x=math.log2(upregulated_fc_threshold), color='#B10DC9', linestyle='--', linewidth=2, alpha=0.5)

    x_min, x_max = comparison_df['Log2_Fold_Change'].min(), comparison_df['Log2_Fold_Change'].max()
    y_min, y_max = (-np.log10(comparison_df['Tukey_HSD_P_Value'])).min(), (-np.log10(comparison_df['Tukey_HSD_P_Value'])).max()
    x_padding = 0.1 * (x_max - x_min)
    y_padding = 0.1 * (y_max - y_min)
    plt.xlim(x_min - x_padding, x_max + x_padding)
    plt.ylim(y_min - y_padding, y_max + y_padding)

    plt.tick_params(axis='both', which='major', labelsize=12)

    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label=f'Significant (FC > {upregulated_fc_threshold})', markerfacecolor='#FF4136', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label=f'Significant (FC < {downregulated_fc_threshold})', markerfacecolor='#2ECC40', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Not Significant', markerfacecolor='#AAAAAA', markersize=10)
    ]
    plt.legend(handles=legend_elements, loc='lower right', fontsize=12)

    texts = []
    for i, row in comparison_df.iterrows():
        if row['Tukey_HSD_P_Value'] < tukey_p_value_threshold and (row['Log2_Fold_Change'] <= math.log2(downregulated_fc_threshold) or row['Log2_Fold_Change'] >= math.log2(upregulated_fc_threshold)):
            texts.append(plt.text(row['Log2_Fold_Change'],
                                  -np.log10(row['Tukey_HSD_P_Value']),
                                  row['Name'],
                                  fontsize=5,
                                  ha='center',
                                  va='center',
                                  fontweight='bold'))

    adjust_text(texts,
                arrowprops=dict(arrowstyle='-', color='#0074D9', lw=0.5),
                expand_points=(1.5, 1.5),
                expand_text=(1.5, 1.5),
                force_text=(0.5, 0.5),
                force_points=(0.2, 0.2),
                autoalign='xy',
                only_move={'points':'y', 'texts':'xy'})

    plt.gca().set_facecolor('#F8F8F8')
    plt.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout()
    output_file = os.path.join(output_dir, f'volcano_plot_{group1}_vs_{group2}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def generate_volcano_csv(df, group1, group2, tukey_p_value_threshold, upregulated_fc_threshold, downregulated_fc_threshold, output_dir='.'):
    comparison_columns_1 = [
        f'Log2_Fold_Change({group1} vs {group2})',
        f'T_test_P_Value({group1} vs {group2})',
        f'Tukey_HSD_P_Value({group1} vs {group2})',
        'Name', 'KEGG', 'HMDB'
    ]

    comparison_columns_2 = [
        f'Log2_Fold_Change({group2} vs {group1})',
        f'T_test_P_Value({group2} vs {group1})',
        f'Tukey_HSD_P_Value({group2} vs {group1})',
        'Name', 'KEGG', 'HMDB'
    ]

    if all(col in df.columns for col in comparison_columns_1):
        comparison_columns = comparison_columns_1
    elif all(col in df.columns for col in comparison_columns_2):
        comparison_columns = comparison_columns_2
        group1, group2 = group2, group1
    else:
        print(f"No comparison data found for {group1} and {group2}.")
        return

    comparison_df = df[comparison_columns]
    comparison_df.columns = ['Log2_Fold_Change', 'T_test_P_Value',  'Tukey_HSD_P_Value',  'Name', 'KEGG', 'HMDB']

    significant = comparison_df[comparison_df['Tukey_HSD_P_Value'] < tukey_p_value_threshold].copy()
    upregulated = significant[significant['Log2_Fold_Change'] > math.log2(upregulated_fc_threshold)].copy()
    downregulated = significant[significant['Log2_Fold_Change'] < math.log2(downregulated_fc_threshold)].copy()

    def get_color(log2_fc):
        if log2_fc > math.log2(upregulated_fc_threshold):
            return '%23ff0000'
        elif log2_fc < math.log2(downregulated_fc_threshold):
            return '%2300ff00'

    upregulated.loc[:, 'Color'] = upregulated['Log2_Fold_Change'].apply(get_color)
    downregulated.loc[:, 'Color'] = downregulated['Log2_Fold_Change'].apply(get_color)

    result_df = pd.DataFrame({
        'Upregulated_Name': upregulated['Name'].tolist() + [''] * (len(downregulated) - len(upregulated)),
        'Upregulated_KEGG': upregulated['KEGG'].tolist() + [''] * (len(downregulated) - len(upregulated)),
        'Upregulated_HMDB': upregulated['HMDB'].tolist() + [''] * (len(downregulated) - len(upregulated)),
        'Upregulated_Color': upregulated['Color'].tolist() + [''] * (len(downregulated) - len(upregulated)),
        'Upregulated_Log2_Fold_Change': upregulated['Log2_Fold_Change'].tolist() + [''] * (
                    len(downregulated) - len(upregulated)),
        'Upregulated_Tukey_HSD': upregulated['Tukey_HSD_P_Value'].tolist() + [''] * (
                    len(downregulated) - len(upregulated)),
        'Downregulated_Name': downregulated['Name'].tolist() + [''] * (len(upregulated) - len(downregulated)),
        'Downregulated_KEGG': downregulated['KEGG'].tolist() + [''] * (len(upregulated) - len(downregulated)),
        'Downregulated_HMDB': downregulated['HMDB'].tolist() + [''] * (len(upregulated) - len(downregulated)),
        'Downregulated_Color': downregulated['Color'].tolist() + [''] * (len(upregulated) - len(downregulated)),
        'Downregulated_Log2_Fold_Change': downregulated['Log2_Fold_Change'].tolist() + [''] * (
                    len(upregulated) - len(downregulated)),
        'Downregulated_Tukey_HSD': downregulated['Tukey_HSD_P_Value'].tolist() + [''] * (
                    len(upregulated) - len(downregulated))
    })
    output_file = os.path.join(output_dir, f'significant_metabolites_{group1}_vs_{group2}.csv')
    result_df.to_csv(output_file, index=False)
    

def generate_all_volcano_plots_and_csvs(statistic_results, group_pairs, tukey_p_value_threshold, upregulated_fc_threshold, downregulated_fc_threshold, base_output_dir='.'):

    print("Gererating volcano plots and CSVs...")
    for group1, group2 in group_pairs:
        comparison_dir = os.path.join(base_output_dir, f'{group1}_vs_{group2}')
        os.makedirs(comparison_dir, exist_ok=True)
        create_volcano_plot(statistic_results, group1, group2, tukey_p_value_threshold, upregulated_fc_threshold, downregulated_fc_threshold, comparison_dir)
        generate_volcano_csv(statistic_results, group1, group2, tukey_p_value_threshold, upregulated_fc_threshold, downregulated_fc_threshold, comparison_dir)
    print("Volcano plots and CSV generation complete.")

def get_hmdb_ids_from_csv(csv_path):
    df = pd.read_csv(csv_path)
    hmdb_ids = pd.concat([df['Upregulated_HMDB'], df['Downregulated_HMDB']]).dropna().unique()
    return [hmdb_id for hmdb_id in hmdb_ids if hmdb_id != '']

def run_enrichment_analysis(hmdb_ids, output_dir, config, group1, group2):
    ea_config = config['enrichment_analysis']
    os.makedirs(output_dir, exist_ok=True)
    # Generate the summary table
    summary_file_path = generate_summary_table(hmdb_ids, id_db, f"{group1}_vs_{group2}", output_dir)

    categories = ['All Categories', 'Metabolic', 'Disease', 'Drug Action']
    for category in categories:
        category_dir = os.path.join(output_dir, category)
        os.makedirs(category_dir, exist_ok=True)

    # Generate enrichment table for all categories
    all_categories_dir = os.path.join(output_dir, 'All Categories')
    enrichment_df = generate_enrichment_table(pathway_db, id_db, hmdb_ids, all_categories_dir)

    # Create visualizations for all categories
    create_bar_chart(
        enrichment_df,
        all_categories_dir,
        ea_config['metabolite_sets_num_for_plots'],
        category_styles
    )
    create_dot_plot(
        enrichment_df,
        all_categories_dir,
        ea_config['metabolite_sets_num_for_plots'],
        category_styles
    )
    create_network_view(
        enrichment_df,
        pathway_db,
        all_categories_dir,
        ea_config['metabolite_sets_num_for_network_view'],
        ea_config['threshold_shared_hmdb'],
        category_styles
    )

    # Generate enrichment tables and visualizations for specific categories
    for category in ['Metabolic', 'Disease', 'Drug Action']:
        category_dir = os.path.join(output_dir, category)
        os.makedirs(category_dir, exist_ok=True)  # Ensure category directory exists

        filtered_pathway_db = {k: v for k, v in pathway_db.items() if pathway_categories.get(k) == category}

        category_enrichment_df = generate_enrichment_table(filtered_pathway_db, id_db, hmdb_ids, category_dir, category)

        create_bar_chart(
            category_enrichment_df,
            category_dir,
            ea_config['metabolite_sets_num_for_plots'],
            category_styles,
            category
        )
        create_dot_plot(
            category_enrichment_df,
            category_dir,
            ea_config['metabolite_sets_num_for_plots'],
            category_styles,
            category
        )
        create_network_view(
            category_enrichment_df,
            filtered_pathway_db,
            category_dir,
            ea_config['metabolite_sets_num_for_network_view'],
            ea_config['threshold_shared_hmdb'],
            category_styles,
            category
        )

    print(f"Enrichment analysis completed for {group1} vs {group2}")
    
def read_csv_with_encoding(data_path, encoding='utf-8'):
    try:
        return pd.read_csv(data_path, encoding=encoding)
    except UnicodeDecodeError:
        print(f"Failed to read with {encoding} encoding.")
        return None
def process_file(data_path):
    encodings_to_try = ['utf-8', 'iso-8859-1', 'latin1', 'ascii', 'cp1252']
    df = None
    successful_encoding = None

    for encoding in encodings_to_try:
        try:
            df = read_csv_with_encoding(data_path, encoding)
            if df is not None:
                successful_encoding = encoding
                print(f"Successfully read file with {encoding} encoding.")
                break
        except Exception as e:
            print(f"Error with {encoding} encoding: {str(e)}")

    if df is None:
        print("Failed to read the file with any encoding.")
        return None

    if successful_encoding.lower() != 'utf-8':
        try:
            df = df.astype(str).apply(lambda x: x.str.encode('utf-8', errors='ignore').str.decode('utf-8'))
            print("Converted data to UTF-8 encoding.")
        except Exception as e:
            print(f"Error converting to UTF-8: {str(e)}")
            return None

    return df

def main():

    # Load configuration from JSON file
    with open('config.json', 'r') as config_file:
        config = json.load(config_file)

    data_dir = config['generate_groups']['output_for_data']
    data_path = os.path.join(data_dir, "groups_data.csv")
    id_db_path = config['global']['id_translation_path']
    df = process_file(data_path)

    # Remove rows with empty names or string 'nan'
    df['Name'] = df['Name'].astype(str).str.strip()
    df = df[~((df['Name'] == '') | (df['Name'].str.lower() == 'nan'))]

    dtype_dict = {f'NAME{i}': str for i in range(1, 282)}
    dtype_dict.update({'KEGG': str, 'HMDB': str, 'CID': str, 'SID': str, 'FORMULA': str, 'SMILES': str})
    columns_to_use = ['KEGG', 'HMDB', 'CID', 'SID', 'FORMULA', 'SMILES'] + [f'NAME{i}' for i in range(1, 282)]

    id_db = pd.read_csv(id_db_path, usecols=columns_to_use, dtype=dtype_dict)
    id_db['FORMULA'] = id_db['FORMULA'].apply(clean_formula)

    try:
        formula_name_to_hmdb, name_to_hmdb = create_formula_name_to_hmdb_dict(id_db)
        hmdb_to_kegg = create_hmdb_to_kegg_dict(id_db)
        print('Dictionaries created successfully.')
    except Exception as e:
        print(f'Error creating dictionaries: {e}')
        return

    # name_to_hmdb = create_name_to_hmdb_dict(id_db)
    # hmdb_to_kegg = create_hmdb_to_kegg_dict(id_db)
    df['FORMULA'] = df['FORMULA'].apply(clean_formula).fillna("")
    # df['HMDB'] = df['Name'].apply(lambda x: get_hmdb(x, df['FORMULA'], formula_name_to_hmdb, name_to_hmdb))
    df['HMDB'] = df.apply(lambda row: get_hmdb(row['Name'], row['FORMULA'], formula_name_to_hmdb, name_to_hmdb), axis=1)
    df['KEGG'] = df['HMDB'].apply(lambda x: get_kegg(x, hmdb_to_kegg))

    df_with_IDs = df.copy()

    # Drop rows where HMDB is NaN
    df_with_IDs = df_with_IDs.dropna(subset=['HMDB'])

    cols = ['Name', 'ID', 'HMDB', 'KEGG'] + [col for col in df_with_IDs.columns if col not in ['Name', 'ID', 'HMDB', 'KEGG']]
    df_with_IDs = df_with_IDs.reindex(columns=cols)

    group_prefixes = extract_group_prefixes(df_with_IDs)

    group_pairs = generate_group_pairs(group_prefixes, config['global']['group_variables'])

    group_cols = generate_group_cols(df_with_IDs, group_prefixes)

    statistic_results = analyze_and_combine_results(df_with_IDs, group_cols, group_pairs)
    cols_order = ['Name', 'ID', 'HMDB', 'KEGG'] + [col for col in statistic_results.columns if col not in ['Name', 'ID', 'HMDB', 'KEGG']]
    statistic_results = statistic_results.reindex(columns=cols_order)

    if not statistic_results.empty:
        statistic_output_path = config['statistic_analysis']['statistic_analysis_output_path']
        statistic_results.to_csv(statistic_output_path, index=False)
        print(f"Saved {len(statistic_results)} rows to {statistic_output_path}")

        # Generate volcano plots
        figs_output_path = config['statistic_analysis']['figs_output_path']
        tukey_p_value_threshold = config['generate_groups']['TukeyHSD_Adjusted_p_value_threshold']
        upregulated_fc_threshold = config['generate_groups']['upregulated_fold_change_threshold']
        downregulated_fc_threshold = config['generate_groups']['downregulated_fold_change_threshold']
        generate_all_volcano_plots_and_csvs(statistic_results, group_pairs, tukey_p_value_threshold, upregulated_fc_threshold, downregulated_fc_threshold, base_output_dir=figs_output_path)
        
        print(f"Gererating enrichment analysis for each group pair...")
        for group1, group2 in group_pairs:
            comparison_dir = os.path.join(figs_output_path, f'{group1}_vs_{group2}')
            enrichment_dir = os.path.join(comparison_dir, 'enrichment_analysis')
            os.makedirs(enrichment_dir, exist_ok=True)

            csv_path = os.path.join(comparison_dir, f'significant_metabolites_{group1}_vs_{group2}.csv')
            hmdb_ids = get_hmdb_ids_from_csv(csv_path)
            
            if hmdb_ids:
                run_enrichment_analysis(hmdb_ids, enrichment_dir, config, group1, group2)
            else:
                print(f"No significant metabolites found for {group1} vs {group2}. Skipping enrichment analysis.")

        print("Statistic analysis, volcano plots, and enrichment analysis completed.")
    else:
        print("No results to save or plot. Skipping output generation.")


if __name__ == '__main__':
    multiprocessing.freeze_support()  
    main()