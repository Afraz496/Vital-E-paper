# Vital-E-paper

This markdown supports the Vital-E paper and is structured as follows:

1. Python: Includes all the Python scripts used for the ML portion of the manuscript (Generates Figure 3 of the manuscript)
2. R: Includes all the R Scripts used for creating the simulated Ct Data (used in Figure 3) and for all the SEIR analysis on the provincial and outbreak data (Figure 4). For a full guide on how to run the software in R please refer to the [virosolver](https://github.com/jameshay218/virosolver) and [virosolver_paper](https://github.com/jameshay218/virosolver_paper) repositories by James Hay.
3. Simulated Data: Includes 3 different sizes of simulated ct values with seed 0.

The remaining repository is structured as follows:

1. Setup
2. Generating Simulated Ct Data
3. Usage
4. Discussion

# Setup

In order to run the software in this repository there are dependencies which are required to run both the pipeline and the compartmental models. These portion of the README divides the dependencies into the Python related/R related dependencies and expands on their installation.

## Python

The environment comes with an accompanying `environment.yml` file to create a virtual environment to run the Machine Learning envrionment.

In project directory run the following:
```bash
conda env create -f environment.yml
```
This creates an environment named `mlenv` (note you can change the environment name by editing the `environment.yml` file). To activate the virtual environment run:

```bash
conda activate mlenv
```

Note that in osx you may instead need to run:
```bash
source activate mlenv
```

for more information on managing virtual environments see this [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment).

## R

For R we require the installation of the virosolver package and all of its dependencies. Please refer to the guide available on the [virosolver repository](https://github.com/jameshay218/virosolver)

# Generating Simulated Ct Data

To generate simulated Ct data you will need to run the R script `generate_ct_for_ml.R`. This script makes use of bash scripting to automate various seeds (of various sizes) and can in turn generate simulated Ct data for use with the Python ML pipeline. Before running the script it would be best to run the Markdown `generate_ct_for_machine_learning.Rmd` to gain some understanding on how the Ct simulation is created. 

Bash script example for Ct simulations:

```bash
vals=($(seq 0 100))
sizes=(100 1000 10000)
for val in ${vals[@]}; do
  for size in ${sizes[@]}; do
    Rscript generate_simulated_ct_for_ml.R ${val} ${size} 
  done
done
```

This script uses the seeds [0, 1, 2, 3, ..., 100] with sample sizes of [100,1000,10000] to create Ct values with an arbitrary time and a generated Rt value. This is then fed into the Machine Learning Pipeline to generate First moments on the Ct values as features and separate the predictor `Rt` into a separate file. 

To help you get started with the usage section there are 3 files of seed 0 with different sample sizes to use in the Machine Learning pipeline.

# Usage

This section describes how to operate the Simulated Ct data in the Machine Learning pipeline and gives some tips on creating compartmental models using the associated R scripts.

## Python

The ML pipeline proceeds in the following order:

1. Fix the Data (rename columns and prepare it)
2. Generate the First Moments (Feature creation)
3. Train on the Data
4. Cross Validate on the data using `TimeSeriesSplit()`
5. Generate Test scores (we generate Rt plots in this pipeline)

Steps 3,4 and 5 are controlled by a single script.

### Fix the Data

If you ran `generate_simulated_ct_for_ml.R` it would generate a file `simulated_ct_output_{SEED}_{SIZE}.csv`. This needs to be supplied as the first argument of the following command:

```bash
python fix_ct_columns.py -i <path_to_generated_ct_file> -f <path_to_output_file> -rt <path_to_rt_file>
```

This then generates 2 files: `simulated_ct_output_{SEED}_{SIZE}.csv` and `rt_simulated_{SEED}_{SIZE}.csv` if you provided these names to the command. This step creates an artificial collection date for the Ct value and segregates the predictor (Rt) from the dataset. This is particularly important to conceptualize because if you are provided real world data it would look like the output of this command.

### Generate the First Moments

```bash
python create_ct_data.py -i <path_to_ct_with_dates_file> -f <path_to_output_file> -rt <path_to_rt_file>
```

### Pipeline

This pipeline also comes equipped with a bootstrapping method (in case your dataset is of a sufficiently small training size). For the simulated data the bootstrapper is unnecessary.

To run the pipeline you will need to use the dataset generated in the last step and run the following command:

```bash
python train_test.py -d <path_to_ct_dataset> -o <path_to_output> -t <train_size> -b <BOOTSTRAP>
```

An example of this would be:

```bash
python train_test.py -d "Data/ct_data.csv" -o "Output" -t 0.8 
```

This command reads the file `Data/ct_data.csv` and sends the pipeline output to a folder called `Output` (if there is no folder, the pipeline will make one for you) and specifies a train size of 0.8 which makes an 80:20 Train Test Split.

## R

Apart from the Ct simulation R scripts, there are accompanying scripts which generate the figures used for the compartmental models in the paper. There are 2 files of interest `outbreak_SEIR.R` and `multiple_cross_sections_SEIR_horizons.R`. 

The `outbreak_SEIR.R` file was used to perform a validation piece on an outbreak in BC. 

The `mutliple_cross_sections_SEIR_horizons.R` file was run on a random population of individuals in British Columbia to predict an Omicron wave.

Both these scripts leverage the [vignette](https://jameshay218.github.io/virosolver/articles/vignette.html) used in the virosolver paper.  

# Discussion

It is important to generate sizeable pseudo-random Ct populations to gauge signficant predictive power in the scope. This method works well on sizeable pseudorandom data and leaves the scope open for the interpretation of how Ct Values epidemic trends. 
