# %%
"""
enrichment_analysis.py
language: Python3
author: L. Xie <xiel4@vcu.edu>
"""

## packages
import os
import pandas as pd
import numpy as np
import json
import scipy.stats as stats
import statsmodels.stats.multitest as multitest

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.text import Annotation
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors

import networkx as nx
from adjustText import adjust_text

with open('config.json', 'r') as config_file:
    config = json.load(config_file)


pathways = pd.read_csv(config['global']['pathway_path'])
pathway_db = {row['pathway names']: row['HMDB IDs'].split(', ') for _, row in pathways.iterrows()}
pathway_categories = dict(zip(pathways['pathway names'], pathways['category']))

unique_categories = pathways['category'].unique()
# Generate a color map
n_colors = len(unique_categories)
color_map = plt.cm.get_cmap('tab10')  # You can change 'tab10' to any other colormap
colors_for_category = [mcolors.rgb2hex(color_map(i)[:3]) for i in range(n_colors)]

# Create category styles dynamically
category_styles = {
    category: ('normal', 'normal', color)
    for category, color in zip(unique_categories, colors_for_category)
}

dtype_dict = {f'NAME{i}': str for i in range(1, 282)}
dtype_dict.update({'KEGG': str, 'HMDB': str, 'CID': str, 'SID': str})
columns_to_use = ['KEGG', 'HMDB', 'CID', 'SID'] + [f'NAME{i}' for i in range(1, 282)]

id_db = pd.read_csv(config['global']['id_translation_path'], usecols=columns_to_use, dtype=dtype_dict)
    
def get_input_hmdb_ids():
    group_variables = config['global']['group_variables']

    if len(group_variables) == 3:
        # 3 variables case
        group_variable = group_variables[2]
        print(f"The group variable is: {group_variable}")
        
        if 'group_variable_values' in config['global']:
            values = config['global']['group_variable_values'].get(group_variable, group_variable.split('/'))
        else:
            values = group_variable.split('/')
        
        while True:
            chosen_value = input(f"Choose a value for {group_variable} ({'/'.join(values)}): ").strip()
            if chosen_value not in values:
                print(f"Invalid value. Please choose from {', '.join(values)}.")
                continue

            file_name = f"factorial_output_{chosen_value}.csv"
            file_path = os.path.join(config['data_selection']['selection_output_path'], file_name)
            
            try:
                df = pd.read_csv(file_path)
                break
            except FileNotFoundError:
                print(f"File {file_path} not found. Please make sure it exists in the output directory.")
                continue
    else:
        # 2 variables case
        file_name = "factorial_output.csv"
        file_path = os.path.join(config['data_selection']['selection_output_path'], file_name)
        
        try:
            df = pd.read_csv(file_path)
            chosen_value = None
        except FileNotFoundError:
            print(f"File {file_path} not found. Please make sure it exists in the output directory.")
            return None, None

        # print("Available groups:")
        # for col in df.columns:
        #     print(col)
    hmdb_columns = [col for col in df.columns if col.endswith('_HMDB')]

    input_hmdb_ids_dict = {}
    for col in hmdb_columns:
        group = col.replace('_HMDB', '')
        input_hmdb_ids_dict[group] = df[col].dropna().tolist()

    return input_hmdb_ids_dict, chosen_value
        
def generate_summary_table(hmdb_ids, id_db, group, output_path):
    # Create a dictionary to store the results
    summary_dict = {'HMDB': [], 'KEGG': [], 'Name': []}
    
    # For each HMDB ID, find the corresponding KEGG ID and Name
    for hmdb_id in hmdb_ids:
        if pd.isna(hmdb_id) or hmdb_id == '':
            continue  # Skip empty or NaN HMDB IDs
        
        summary_dict['HMDB'].append(hmdb_id)
        
        # Find the corresponding row in id_db
        id_row = id_db[id_db['HMDB'] == hmdb_id]
        
        if not id_row.empty:
            summary_dict['KEGG'].append(id_row['KEGG'].iloc[0])
            summary_dict['Name'].append(id_row['NAME1'].iloc[0])
        else:
            summary_dict['KEGG'].append('N/A')
            summary_dict['Name'].append('N/A')
    
    # Create a DataFrame from the dictionary
    summary_df = pd.DataFrame(summary_dict)
    
    # Sort the DataFrame by HMDB ID
    summary_df = summary_df.sort_values('HMDB').reset_index(drop=True)
    
    # Remove any rows where all values are 'N/A'
    summary_df = summary_df[~(summary_df == 'N/A').all(axis=1)]
    
    # Save the summary table
    summary_output_path = os.path.join(output_path, f'{group}_names.csv')
    summary_df.to_csv(summary_output_path, index=False)
    
    return summary_output_path

def generate_enrichment_table(pathway_db, id_db, input_hmdb_ids, output_path, category=None):
    hmdb_to_kegg = dict(zip(id_db['HMDB'], id_db['KEGG']))
    hmdb_to_name = dict(zip(id_db['HMDB'], id_db['NAME1']))

    unique_metabolites = set([met for sublist in pathway_db.values() for met in sublist])
    N = len(unique_metabolites)
    n = len(input_hmdb_ids)

    results = []
    for pathway, metabolites in pathway_db.items():
        K = len(metabolites)
        if K >= config['enrichment_analysis']['metabolite_sets_num']:
            observed_metabolites_set = set(metabolites) & set(input_hmdb_ids)
            observed_metabolites = len(observed_metabolites_set)
            if observed_metabolites > 0:
                expected_value = (K * n) / N
                p_value = 1 - stats.hypergeom.cdf(observed_metabolites - 1, N, K, n)
                hits_detail_kegg = ', '.join(hmdb_to_kegg.get(met, 'N/A') for met in observed_metabolites_set)
                hits_detail_name = ', '.join(hmdb_to_name.get(met, 'N/A') for met in observed_metabolites_set)
                pathway_detail = ', '.join(metabolites)
                results.append({
                    'Pathway': pathway,
                    'Category': pathway_categories.get(pathway, 'N/A'),
                    'Total': K,
                    'Hits': observed_metabolites,
                    'Expect': expected_value,
                    'Enrichment Ratio': observed_metabolites / expected_value,
                    'P Value': p_value,
                    'Hits Detail (HMDB)': ', '.join(observed_metabolites_set),
                    'Hits Detail (KEGG)': hits_detail_kegg,
                    'Hits Detail (Name)': hits_detail_name,
                    'Pathway Metabolites': pathway_detail
                })
    if not results:
        return pd.DataFrame()
        
    results_df = pd.DataFrame(results)
    df_hits_greater_than_1 = results_df[results_df['Hits'] > 1].sort_values(by='P Value').reset_index(drop=True)
    df_hits_less_or_equal_1 = results_df[results_df['Hits'] <= 1].sort_values(by='P Value').reset_index(drop=True)
    sorted_filtered_results_df = pd.concat([df_hits_greater_than_1, df_hits_less_or_equal_1]).reset_index(drop=True)

    holm_p = multitest.multipletests(sorted_filtered_results_df['P Value'], alpha=0.05, method='holm')[1]
    sorted_filtered_results_df['Holm P'] = holm_p

    fdr_p = multitest.multipletests(sorted_filtered_results_df['P Value'], alpha=0.05, method='fdr_bh')[1]
    sorted_filtered_results_df['FDR'] = fdr_p

    cols = ['Pathway', 'Category'] + [col for col in sorted_filtered_results_df.columns if col not in ['Pathway', 'Category'] and not col.startswith('Hits Detail') and col != 'Pathway Metabolites']
    cols += [col for col in sorted_filtered_results_df.columns if col.startswith('Hits Detail')] + ['Pathway Metabolites']
    sorted_filtered_results_df = sorted_filtered_results_df[cols]

    sorted_filtered_results_df['Pathway Metabolites'] = sorted_filtered_results_df.apply(
        lambda row: highlight_hits(row['Pathway Metabolites'], row['Hits Detail (HMDB)']),
        axis=1
    )

    if category:
        csv_file_path = os.path.join(output_path, f'enrichment_analysis_results_{category}.csv')
    else:
        csv_file_path = os.path.join(output_path, 'enrichment_analysis_results.csv')

    sorted_filtered_results_df.to_csv(csv_file_path, index=False)
    return sorted_filtered_results_df

def highlight_hits(pathway_metabolites, hits_detail_hmdb):
    hits_set = set(hits_detail_hmdb.split(', '))
    highlighted = [f'{met}*' if met in hits_set else met for met in pathway_metabolites.split(', ')]
    return ', '.join(highlighted)


def create_bar_chart(df, output_path, metabolite_sets_num_for_plots, category_styles, suffix=''):
    if df.empty:
        print("Warning: No data available to plot.")
        return
    df['Enrichment Ratio'] = df['Hits'] / df['Expect']
    df_to_plot = df.sort_values(by='P Value', ascending=True).head(metabolite_sets_num_for_plots)
    df_to_plot = df_to_plot.sort_values(by='P Value', ascending=False)
    num_metabolite_sets = df_to_plot.shape[0]
    fig_height = max(10, num_metabolite_sets * 0.3)  # Increased height factor

    fig, ax = plt.subplots(figsize=(18, fig_height))  # Increased width for legend
    colors = cm.autumn_r(df_to_plot['Enrichment Ratio'] / df_to_plot['Enrichment Ratio'].max())
    bars = ax.barh(range(len(df_to_plot)), df_to_plot['Enrichment Ratio'], color=colors)

    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    ax.set_xlabel('Enrichment Ratio', fontsize=12)
    ax.set_title(f"Overview of Enriched Metabolite Sets (Top {metabolite_sets_num_for_plots})",
                 fontsize=14, fontweight='bold')

    # Adjust plot to make room for labels and legend
    plt.subplots_adjust(left=0.4, right=0.75)  # Adjusted right margin for legend

    # Add labels to the right of bars
    for i, (bar, hits, total) in enumerate(zip(bars, df_to_plot['Hits'], df_to_plot['Total'])):
        ax.text(bar.get_width(), i,
                f'{bar.get_width():.2f} ({hits}/{total})',
                ha='left', va='center', fontweight='bold', fontsize=10)

    # Change font based on category
    for i, (pathway, category) in enumerate(zip(df_to_plot['Pathway'], df_to_plot['Category'])):
        style = category_styles.get(category, ('normal', 'normal', 'black'))
        ax.text(-0.05, i, pathway, ha='right', va='center', fontsize=10,
                fontweight=style[0], fontstyle=style[1], color=style[2],
                transform=ax.get_yaxis_transform())

    # Remove y-axis ticks and labels
    ax.set_yticks([])
    ax.set_yticklabels([])

    # Remove the y-axis label that's too close to the axis
    ax.set_ylabel('')
    
    # Add "Metabolite Sets" label to the left of pathway names
    # Using a very small x position to ensure it's far to the left and won't overlap
    fig.text(0.02, 0.5, 'Metabolite Sets', va='center', rotation='vertical', 
             fontsize=12)  # No bold formatting

    # Create a separate axes for the legend
    legend_ax = fig.add_axes([0.78, 0.1, 0.2, 0.8])  # [left, bottom, width, height]
    legend_ax.axis('off')

    # Add legend for category styles
    legend_elements = [Line2D([0], [0], color=style[2], lw=3, label=category)
                       for category, style in category_styles.items()]
    legend = legend_ax.legend(handles=legend_elements, loc='center left', title='Categories',
                              fontsize=10, title_fontsize=12)
    legend.get_frame().set_alpha(0.8)  # Semi-transparent background
    legend.get_frame().set_edgecolor('lightgray')  # Light border

    # Ensure all bars and labels are within the figure
    ax.set_xlim(right=ax.get_xlim()[1] * 1.2)  # Extend x-axis limit by 20%

    bar_chart_output_path = os.path.join(output_path, f'bar_chart{suffix}.png')
    plt.savefig(bar_chart_output_path, bbox_inches='tight', dpi=300)
    plt.close()


def create_dot_plot(df, output_path, metabolite_sets_num_for_plots, category_styles, suffix=''):
    if df.empty:
        print("Warning: No data available to plot.")
        return
    df_to_plot = df.sort_values(by='P Value', ascending=True).head(metabolite_sets_num_for_plots)
    df_to_plot = df_to_plot.sort_values(by='P Value', ascending=False)
    df_to_plot['-log10(p-value)'] = -np.log10(df_to_plot['P Value'])

    # Calculate appropriate figure height and width
    num_metabolite_sets = df_to_plot.shape[0]
    # Increase height factor to ensure enough space
    fig_height = max(12, num_metabolite_sets * 0.25)
    # Increase width to accommodate long pathway names
    fig_width = 18
    
    # Create figure with increased width
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Calculate left margin based on maximum pathway name length
    max_pathway_length = max(len(pathway) for pathway in df_to_plot['Pathway'])
    # Ensure left margin is proportional to the longest pathway name
    left_margin = min(0.5, max(0.3, max_pathway_length * 0.015))  # Increased minimum margin
    
    sizes = df_to_plot['Enrichment Ratio'] * 5
    colors = df_to_plot['P Value']

    scatter = ax.scatter(
        df_to_plot['-log10(p-value)'],
        df_to_plot['Pathway'],
        s=sizes,
        c=colors,
        cmap='autumn',
        edgecolor='k',
        alpha=0.6,
        marker='o'
    )

    ax.grid(True, linestyle='--', linewidth=0.5)
    ax.set_xlabel('-log10(p-value)')
    ax.set_title(f"Overview of Enriched Metabolite Sets (Top {metabolite_sets_num_for_plots})", fontsize=14, fontweight='bold')

    # Colorbar for P-value
    cbar = plt.colorbar(scatter)
    cbar.set_label('P-value')

    # Adjust y-axis labels with category colors
    ax.set_yticks(range(len(df_to_plot)))
    ax.set_yticklabels([])
    
    # Remove the y-axis label that's too close to the axis
    ax.set_ylabel('')
    
    # Use the same font size as in bar chart (10) for pathway names
    # This places text at a fixed distance to the left of the y-axis
    for i, (pathway, category) in enumerate(zip(df_to_plot['Pathway'], df_to_plot['Category'])):
        style = category_styles.get(category, ('normal', 'normal', 'black'))
        ax.text(-0.05, i, pathway, ha='right', va='center',
                fontsize=10,  # Fixed font size to match bar chart
                fontweight=style[0], fontstyle=style[1], color=style[2],
                transform=ax.get_yaxis_transform())  # Use y-axis transform for consistent placement
    
    # Add "Metabolite Sets" label to the left of pathway names
    # Calculate label position: closer to the middle of pathway names
    # For dot plots, we use a fixed left margin of 0.4 (from subplots_adjust)
    # The 0.03 offset ensures there's just enough space between the label and pathway names
    label_x_pos = max(0.05, min(0.2, 0.4 - 0.03))
    
    # Position the label vertically in the middle of the figure
    fig.text(label_x_pos, 0.5, 'Metabolite Sets', va='center', rotation='vertical', 
             fontsize=12)  # Removed bold formatting

    # Create separate legends for Categories and Enrichment Ratio
    category_elements = [Line2D([0], [0], color=style[2], lw=2, label=category)
                         for category, style in category_styles.items()]

    enrichment_elements = [Line2D([0], [0], marker='o', color='w',
                                  label=f'Enrichment Ratio: {size}',
                                  markerfacecolor='gray',
                                  markersize=np.sqrt(size * 5))
                           for size in [10, 20, 30]]

    # Add Categories legend
    category_legend = ax.legend(handles=category_elements,
                                loc='lower right',
                                title='Categories',
                                fontsize=8,
                                title_fontsize=10,
                                bbox_to_anchor=(1, 0.1),
                                bbox_transform=ax.transAxes)
    ax.add_artist(category_legend)

    # Add Enrichment Ratio legend below Categories legend
    enrichment_legend = ax.legend(handles=enrichment_elements,
                                  loc='lower right',
                                  title='Enrichment Ratio',
                                  fontsize=8,
                                  title_fontsize=10,
                                  bbox_to_anchor=(1, 0),
                                  bbox_transform=ax.transAxes)
    ax.add_artist(enrichment_legend)

    # Use subplots_adjust instead of tight_layout for more precise margin control
    plt.subplots_adjust(left=left_margin, right=0.85, top=0.95, bottom=0.1)
    
    # Save figure with increased margins
    dot_plot_output_path = os.path.join(output_path, f'dot_plot{suffix}.png')
    # Use a larger pad_inches value to ensure text is not cut off
    plt.savefig(dot_plot_output_path, bbox_inches='tight', dpi=300, pad_inches=0.8)
    plt.close()


def create_network_view(df, pathway_db, output_path, metabolite_sets_num_for_network_view, threshold_shared_hmdb,
                        category_styles, suffix=''):
    if df.empty:
        print("Warning: No data available to plot.")
        return
    G = nx.Graph()
    df_to_network = df.head(metabolite_sets_num_for_network_view)
    for i, row in df_to_network.iterrows():
        G.add_node(row['Pathway'], size=row['Enrichment Ratio'] * 10, pvalue=row['P Value'], category=row['Category'])
    for i, row1 in df_to_network.iterrows():
        pathway1 = row1['Pathway']
        hmdb_ids1 = set(pathway_db[pathway1])
        for j, row2 in df_to_network.iterrows():
            if i < j:
                pathway2 = row2['Pathway']
                hmdb_ids2 = set(pathway_db[pathway2])
                shared_hmdb = hmdb_ids1.intersection(hmdb_ids2)
                if len(shared_hmdb) > threshold_shared_hmdb:
                    G.add_edge(pathway1, pathway2)

    fig, ax = plt.subplots(figsize=(20, 15))

    pos = nx.spring_layout(G, k=1.5, iterations=50)
    sizes = [G.nodes[node]['size'] for node in G.nodes()]
    colors = [G.nodes[node]['pvalue'] for node in G.nodes()]

    nx.draw_networkx_edges(G, pos, alpha=0.4, ax=ax, edge_color='gray', width=0.5)

    nodes = nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color=colors, cmap='autumn', alpha=0.8)

    # Create text labels with colors based on categories
    texts = []
    for node, (x, y) in pos.items():
        category = G.nodes[node]['category']
        style = category_styles.get(category, ('normal', 'normal', 'black'))
        texts.append(ax.text(x, y, node, fontsize=10, fontweight=style[0],
                             fontstyle=style[1], color=style[2],
                             ha='center', va='center'))

    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', alpha=0.3),
                expand_points=(1.2, 1.2), force_points=(0.5, 0.5))

    sm = plt.cm.ScalarMappable(cmap='autumn', norm=plt.Normalize(vmin=min(colors), vmax=max(colors)))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label='P-value')

    # Add Enrichment Ratio legend
    for size in [10, 20, 30]:
        plt.scatter([], [], c='k', alpha=0.3, s=size * 10, label=f'Enrichment Ratio: {size}')
    enrichment_legend = plt.legend(scatterpoints=1, frameon=False, labelspacing=1, loc='lower right',
                                   title='Enrichment Ratio')
    plt.gca().add_artist(enrichment_legend)

    # Add Category legend
    category_elements = [Line2D([0], [0], color=style[2], lw=2, label=category)
                         for category, style in category_styles.items()]
    category_legend = plt.legend(handles=category_elements, title='Categories',
                                 loc='upper right', fontsize=8, title_fontsize=10)
    plt.gca().add_artist(category_legend)

    plt.title('Network of Enriched Metabolite Sets', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.margins(0.1)

    network_view_output_path = os.path.join(output_path, f'network_view{suffix}.png')
    plt.savefig(network_view_output_path, bbox_inches='tight', dpi=300)
    plt.close()

def filter_pathways_by_category(pathway_db, pathway_categories, category):
    filtered_db = {
        pathway: metabolites
        for pathway, metabolites in pathway_db.items()
        if pathway_categories.get(pathway) == category
    }
    num_pathways = len(filtered_db)
    
    return filtered_db, num_pathways

def create_visualizations(df, pathway_db, output_path, config, category_styles, category=None):
    suffix = f'_{category}' if category else ''
    create_bar_chart(df, output_path, config['enrichment_analysis']['metabolite_sets_num_for_plots'], category_styles, suffix)
    create_dot_plot(df, output_path, config['enrichment_analysis']['metabolite_sets_num_for_plots'], category_styles, suffix)
    create_network_view(df, pathway_db, output_path,
                        config['enrichment_analysis']['metabolite_sets_num_for_network_view'],
                        config['enrichment_analysis']['threshold_shared_hmdb'], category_styles, suffix)

def main():
    print("Getting input HMDB IDs...")
    input_hmdb_ids_dict, chosen_value = get_input_hmdb_ids()

    if input_hmdb_ids_dict is None:
        print("Error: Could not load input HMDB IDs. Exiting.")
        return

    enrichment_analysis_path = config['enrichment_analysis']['output_path']
    os.makedirs(enrichment_analysis_path, exist_ok=True)

    categories_to_analyze = ['Metabolic', 'Disease', 'Drug Action']
    for group, input_hmdb_ids in input_hmdb_ids_dict.items():
        if chosen_value:
            output_folder = f"enrichment_analysis_{chosen_value}_{group}"
        else:
            output_folder = f"enrichment_analysis_{group}"

        output_path = os.path.join(config['enrichment_analysis']['output_path'], output_folder)
        os.makedirs(output_path, exist_ok=True)

        print(f"Processing enrichment analysis for {group}:")
        print(f"Results will be saved in: {output_path}")
  
        # Generate summary table
        summary_file_path = generate_summary_table(input_hmdb_ids, id_db, group, output_path)
        all_categories_folder = os.path.join(output_path, "All Categories")
        os.makedirs(all_categories_folder, exist_ok=True)

        print("Performing analysis for all categories...")
        enrichment_df = generate_enrichment_table(pathway_db, id_db, input_hmdb_ids, all_categories_folder)
        create_visualizations(enrichment_df, pathway_db, all_categories_folder, config, category_styles)

        for category in categories_to_analyze:
            print(f'Performing analysis for {category} category...')
            category_output_path = os.path.join(output_path, category)
            os.makedirs(category_output_path, exist_ok=True)

            filtered_pathway_db, num_pathways = filter_pathways_by_category(pathway_db, pathway_categories, category)

            if num_pathways < 5:
                print(f"Less than 5 pathways found for {category} category. Skipping analysis.")
                continue


            enrichment_df = generate_enrichment_table(filtered_pathway_db, id_db, input_hmdb_ids, category_output_path,
                                                      category)

            if len(enrichment_df) < 5:
                print(f"Less than 5 enriched pathways found for {category} category. Skipping visualization.")
                continue
            create_visualizations(enrichment_df, filtered_pathway_db, category_output_path, config, category_styles,
                                  category)

    print("All enrichment analyses completed!")


if __name__ == "__main__":
    main()