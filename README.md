# PyWhy HIV Simulation

This repository demonstrates causal inference methods from the PyWhy framework applied to simulated HIV treatment data generated using the WhyNot library.

## Project Overview

This project uses simulated HIV treatment data to explore and test causal inference methods from the PyWhy framework. By using synthetic data with known ground truth causal effects, we can evaluate the performance of different causal inference techniques under various conditions.

The simulation creates datasets that model HIV treatment scenarios with the following variations:
- Full vs. partial treatment compliance
- Presence vs. absence of confounding factors

## Datasets

The simulation generates four different datasets, each representing a different experimental scenario:

1. **Full Treatment Compliance**: All patients in the treatment group receive the treatment, and no patients in the control group receive it.
2. **Partial Treatment Compliance**: Some patients in the treatment group do not take the treatment, and some in the control group may get access to it.
3. **Full Treatment Compliance with Confounding**: Treatment compliance is affected by confounding factors like immune response and virus levels.
4. **Partial Treatment Compliance with Confounding**: Combines partial compliance with confounding effects.

Each dataset contains information about:
- Patient covariates (immune response, virus levels, etc.)
- Treatment assignment
- Treatment outcome (infected macrophages at the end of observation)
- True causal effect (known from the simulation)

## Setup

### Environment Setup

This project uses conda for environment management. There are two separate environments with one for general analysis
and one for creating the simulated datasets. The reason for this is that WhyNot package relies on old version of certain
libraries (e.g. mesa) that would force us to use old version of the PyWhy packages

```bash
# Create the analysis environment
conda env create -f environment.yml

# Activate the environment
conda activate pywhy-hiv-simulation
```

```bash
# Create the data creation environment
conda env create -f environment.yml

# Activate the environment
conda activate whynot

```

### Running the Simulation

The simulated data is created by single script. If you want to replicate the analysis keep the parameters in this
file as they are, but you can also play around with them to see what happens to the results

```bash
# Activate the environment
conda activate whynot

# Run the data creation script
python create_simulated_data.py
```

## Analysis

Analyses using PyWhy are in 4 different notebook, with one for each type of experiment setup.
