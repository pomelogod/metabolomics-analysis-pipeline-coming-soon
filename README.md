# Metabolomic Analysis Pipeline - Metaboenrich


This repository contains a set of Python scripts for performing metabolomics data analysis. The pipeline includes data preprocessing, statistical analysis, heatmap visualization, and enrichment analysis.

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Setup](#setup)
3. [Configuration](#configuration)
4. [Scripts](#scripts)
5. [Usage](#usage)
6. [Output](#output)
7. [Examples](#examples)
8. [Notes](#notes)
8. [Author](#author)

## Prerequisites

- Anaconda or Miniconda
- Git (optional, for cloning the repository)
- Python 3.8

## Setup

1. Clone the repository (or download and extract the ZIP file):
   ```sh
   git clone https://github.com/yourusername/metabolomics-analysis.git
   cd metabolomics-analysis
   ```

2. Create a Conda environment with Python 3.8:
   ```sh
   conda create --name myenv python=3.8
   ```

3. Activate the environment:
   ```sh
   conda activate myenv
   ```

4. Install the required packages:
   ```sh
   pip install -r requirements.txt
   ```

## Configuration

Before running the scripts, review the `config.json` file which contains all the necessary parameters for the analysis. Ensure you update the required fields if needed. All results will be saved by default in the `output` folder inside the workflow directory. The `output` folder will be automatically created when you run the scripts:

- `global`: Set the paths for pathways database and ID translation database, and define the group variables for analysis
- `generate_groups`: Configure the number of values for each group, specify the output directory, and set thresholds for significance and fold change
- `statistic_analysis`: Set the output paths for statistical analysis results and figures from volcano plots
- `heatmap`: Configure clustering options, output paths for the heatmap visualization, and datasets derived from the heatmap
- `data_selection`: Configure the selection criteria and output path for significant metabolites
- `enrichment_analysis`: Set the parameters for enrichment analysis, including threshold values, visualization options, and output paths

## Scripts

1. `generate_groups.py`: Generates group data framework based on the configuration and creates an explanation of selection criteria
2. `statistic_analysis.py`: Performs statistical analysis on the group data including ANOVA, Tukey HSD, fold change calculations, generates volcano plots, and runs enrichment analysis for significant metabolites from each group pair
3. `heatmap.py`: Creates hierarchical clustered heatmaps of fold changes across comparisons, identifies patterns of significant metabolites, and performs enrichment analysis for each pattern group
4. `data_selection.py`: Selects significant metabolites based on statistical results and selection criteria for further analysis
5. `enrichment_analysis.py`: Performs pathway enrichment analysis on selected metabolites and generates visualizations including bar charts, dot plots, and network views
6. `genes_ranked.py`: Ranks genes based on statistical significance and fold change, and provides a prioritized list for targeted analysis

## Usage

Run the scripts in the following order:

1. Generate group data:
   ```sh
   python generate_groups.py
   ```
   You need to copy your data into the generated `groups_data.csv` file and ensure that your data is in the correct columns.


2. Perform statistical analysis, generate volcano plots and perform enrichment analysis for dataset selected from volcano plots:
   ```sh
   python statistic_analysis.py
   ```

3. Generate heatmap and perform enrichment analysis for dataset selected from heatmap:
   ```sh
   python heatmap.py
   ```
   This script only supports two variables with two unique values each; if your dataset does not meet this requirement, please modify it or skip using this script.

4. Select significant data based on specific criteria:
   ```sh
   python data_selection.py
   ```

5. Perform enrichment analysis on selected datasets:
   ```sh
   python enrichment_analysis.py
   ```
   During this step, you may be prompted to choose specific data for analysis if you have three group variables.

6. Perform gene ranking analysis:
   ```sh
   python genes_ranked.py
   ```
   During this step, you may be prompted to enter the CSV file path for the enrichment analysis result—ensure you provide the file path, not the folder path.
## Output

- `generate_groups.py`: Creates a CSV file with group data structure and a text file explaining selection criteria
- `statistic_analysis.py`: Produces a CSV file with statistical analysis results, volcano plots, and enrichment analysis results for each group pair
- `heatmap.py`: Generates a clustered heatmap visualization showing fold changes across comparisons, identifies patterns of metabolites with similar expression profiles, saves separate CSV files for each pattern group, and performs enrichment analysis for each identified group
- `data_selection.py`: Generates CSV files with selected significant metabolites based on defined criteria
- `enrichment_analysis.py`: Creates CSV files with enrichment analysis results and various plots (bar chart, dot plot, and network view) for different pathway categories
- `genes_ranked.py`: Produces ranked lists of genes based on statistical significance and fold change magnitude for targeted analysis

Check the output directories specified in the `config.json` file for the results of each step.

## Examples

Here are examples of running each script and their expected output:

1. Generate group data:
   ```
   User @ MacBook-Pro ➜  workflow  python generate_groups.py
   DataFrame saved to ./output/groups_data.csv
   Explanation saved to ./output/selection_criteria_explanation.txt
   ```

2. Perform statistical analysis:
   ```
   User @ MacBook-Pro ➜  workflow  python statistic_analysis.py
   Starting combined analysis for 1000 rows...
   Processing rows: 100%|██████████| 1000/1000 [00:30<00:00, 33.33it/s]
   Statistical analysis complete. Combining results...
   Final number of rows in statistic result after removing duplicates: 239
   Analysis and combination complete.
   Saved 239 rows to ./output/statistic_analysis.csv
   Gererating volcano plots and CSVs...
   Volcano plots and CSV generation complete.
   Gererating enrichment analysis for each group pair...
   Enrichment analysis completed for group_pair_1 vs group_pair_2
   ...
   Statistic analysis, volcano plots, and enrichment analysis completed.
   ```

3. Generate heatmap and perform enrichment analysis:
   ```
   User @ MacBook-Pro ➜  workflow  python heatmap.py
   Clustered gene table saved as './output/clustered_genes_log2fc_significance.csv'
   
   Processing Group A...
   Completed enrichment analysis for Group A
   
   Processing Group B...
   Completed enrichment analysis for Group B
   
   Processing Group C...
   No valid HMDB IDs found for Group C
   
   All analyses completed: Heatmap and enrichment analyses have been saved
   ```

4. Select significant data:
   ```
   User @ MacBook-Pro ➜  workflow  python data_selection.py
   Performing data selection...
   Data selection complete.
   Saving results to CSV files in ./output/data_selection/...
   Results saved to CSV files.
   ```

5. Perform enrichment analysis:
   ```
   User @ MacBook-Pro ➜  workflow  python enrichment_analysis.py
   Getting input HMDB IDs...
   The group variable is: treatment/control
   Choose a value for treatment/control (treatment/control): treatment
   Processing enrichment analysis for group1:
   Results will be saved in: ./output/enrichment_analysis/enrichment_analysis_treatment_group1
   Performing analysis for all categories...
   Performing analysis for Metabolic category...
   Performing analysis for Disease category...
   Performing analysis for Drug Action category...
   Processing enrichment analysis for group2:
   Results will be saved in: ./output/enrichment_analysis/enrichment_analysis_treatment_group2
   Performing analysis for all categories...
   Performing analysis for Metabolic category...
   Performing analysis for Disease category...
   Performing analysis for Drug Action category...
   Processing enrichment analysis for group3:
   Results will be saved in: ./output/enrichment_analysis/enrichment_analysis_treatment_group3
   Performing analysis for all categories...
   Performing analysis for Metabolic category...
   Less than 5 enriched pathways found for Metabolic category. Skipping visualization.
   Performing analysis for Disease category...
   Performing analysis for Drug Action category...
   Less than 5 enriched pathways found for Drug Action category. Skipping visualization.
   Processing enrichment analysis for group4:
   Results will be saved in: ./output/enrichment_analysis/enrichment_analysis_treatment_group4
   Performing analysis for all categories...
   Performing analysis for Metabolic category...
   Performing analysis for Disease category...
   Performing analysis for Drug Action category...
   All enrichment analyses completed!
   Enrichment analysis done!
   ```

6. Perform gene ranking analysis:
   ```
   User @ MacBook-Pro ➜  workflow  python genes_ranked.py
   Please enter the path to your enrichment analysis results CSV file: 'PATH/OF/YOUR/DATA'

   ✅ Successfully loaded all input files.
   ✅ Gene ranking saved to: 'PATH/OF/YOUR/DATA'
   ```

## Notes

- Make sure to activate the Conda environment (`conda activate myenv`) before running any scripts.
- This pipeline supports analyses involving 1, 2, or 3 variables. The heatmap.py script is designed to handle 2 variables with 2 values per variable.
- When specifying group variables in your `config.json` file, avoid using underscores (`_`) in the variable names. 

**Correct Format**: Use concatenated names without separators. Example: `["MCM/CerS5cKO", "CD/HFD"]`  
**Incorrect Format**: Avoid underscores. Example: `["MCM/CerS5_cKO", "CD/HFD"]`

- Review and update the `config.json` file according to your specific dataset and analysis requirements.
- After finishing your work, don't forget to deactivate the virtual environment:
   ```sh
   conda deactivate
   ```
- If you encounter any issues or have questions, please open an issue in the GitHub repository.

## Author

    •   Longsheng Xie
    •   Email: longsheng.xie41@gmail.com
