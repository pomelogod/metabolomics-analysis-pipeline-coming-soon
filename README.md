# Metabolomics Analysis Pipeline

This repository contains a set of Python scripts for performing metabolomics data analysis. The pipeline includes data preprocessing, statistical analysis, and enrichment analysis.

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
   pip install -r requirments.txt
   ```

## Configuration

Before running the scripts, review the `config.json` file which contains all the necessary parameters for the analysis. Ensure you update the required fields if you needed. All the result will saved by default in `output`folder inside the workflow directory. The `output` folder will be automatically created when you run the scripts:

- `global`: Set the paths for pathways database and id translation database and the data variables
- `generate_groups`: Configure the number of values for each group, specify the output directory, and set the threshold for selection criteria explanation
- `statistic_analysis`: Set the output paths for statistical analysis results and figures
- `data_selection`: Configure the selection criteria and output path
- `enrichment_analysis`: Set the parameters for enrichment analysis and output path

## Scripts

1. `generate_groups.py`: Generates group data based on the configuration
2. `statistic_analysis.py`: Performs statistical analysis on the generated group data, generates volcano plots and run enrichment analysis for each group pairs
3. `data_selection.py`: Selects significant data based on statistical results
4. `enrichment_analysis.py`: Performs enrichment analysis on the selected data

## Usage

Run the scripts in the following order:

1. Generate group data:
   ```sh
   python generate_groups.py
   ```
   You need to copy your data into the generated `groups_data.csv` file and ensure that your data is in the correct columns.


2. Perform statistical analysis:
   ```sh
   python statistic_analysis.py
   ```

3. Select significant data:
   ```sh
   python data_selection.py
   ```

4. Perform enrichment analysis:
   ```sh
   python enrichment_analysis.py
   ```
   During this step, you may be prompted to choose specific data for analysis if you have three group variables.

## Output

- `generate_groups.py`: Creates a CSV file with group data and a text file explaining selection criteria
- `statistic_analysis.py`: Produces a CSV file with statistical analysis results, volcano plots and enrichment analysis results for each group pairs
- `data_selection.py`: Generates CSV files with selected significant data
- `enrichment_analysis.py`: Creates CSV files with enrichment analysis results and various plots (bar chart, dot plot, and network view)

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

3. Select significant data:
   ```
   User @ MacBook-Pro ➜  workflow  python data_selection.py
   Performing data selection...
   Data selection complete.
   Saving results to CSV files in ./output/data_selection/...
   Results saved to CSV files.
   ```

4. Perform enrichment analysis:
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

## Notes

- Make sure to activate the Conda environment (`conda activate myenv`) before running any scripts.
- This pipeline only supports analyses involving 1, 2 or 3 variables. Using more than 3 variables will result in errors or inaccurate results.
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
    •   Email: xiel4@vcu.edu/longsheng.xie41@gmail.com
