# PyWhy HIV Simulation

This repository demonstrates causal inference methods from the PyWhy framework applied to simulated HIV treatment data generated using the WhyNot library.

## Project Overview

This project uses simulated HIV treatment data to explore and test causal inference methods from the PyWhy framework. By using synthetic data with known ground truth causal effects, we can evaluate the performance of different causal inference techniques under various conditions.

The simulation creates datasets that model HIV treatment scenarios with the following variations:
- Full vs. partial treatment compliance
- Presence vs. absence of confounding factors

## Datasets

We are simulation HIV treatment experiments, where people with HIV infection randomized into treatment or control groups and then
we follow them for X days. We care about the numbers at the end of treatment period.

# How HIV infection and treatment works

We can think of HIV infection as an interplay of five factions:

- T1 – Helper T‑cells are the main target for the HIV virus wants to infect. Spawn at a steady rate, die naturally, or get turned into zombies if hit by virus.
- T2 – Macrophages are additional targets for the HIV that are a bit harder to protect with treatment.
- Infected T1 / T2 – HIV infected cells T-cells and macrophages. Once infected they start pumping out new virus particles, then eventually die on their own or when attacked by the immune killer T-cells
- V – Free virus flowing around the body looking for cells to infect. Disappears when it finds and infects a healthy cell, but every infected cell creates more free virus
- E – Killer T‑cells. Police force that tries to kill infected cells and multiplies faster the more infected cells there are

What treatment does is to gives targets cells better protection against the free virus. There is a twist though, that explains why HIV is so hard to treat: treatment works better on T1 than on T2 , which means macrophages can hide a low‑level infection even when most T1 are safe, making the virus almost impossible to eradicate.

# The Experiments

The simulation generates four different datasets, each representing a different experimental scenario:

1. **Full Treatment Compliance**: All patients in the treatment group receive the treatment, and no patients in the control group receive it.
2. **Partial Treatment Compliance**: Some patients in the treatment group do not take the treatment, and some in the control group may get access to it.
3. **Full Treatment Compliance with Confounding**: Treatment compliance is affected by confounding factors like immune response and virus levels.
4. **Partial Treatment Compliance with Confounding**: Combines partial compliance with confounding effects.

# Datasets

Covariates at the start of experiment:
- uninfected_t1 = Uninfected T-cells (cells/ml)
- infected_t1 = Infected T-cells (cells/ml)
- uninfected_t2 = Uninfected macrophages (cells/ml)
- infected_t2 = Infected macrophages (cells/ml)
- free_virus = Free virus in bloodstream (copies/ml)
- immune_response = Killer T-cells (cells/ml)

The information about confounders comes from simulated lab-measurements which are not
entirely accurate. This means that there is noise in all these measurements

Other:
- enrolled = was the patient enrolled to the treatment group
- treatment = did the patient take the treatment
- outcome = The change in Infected macrophages (infected_t2) from the start to the end of experiment
- true_effect = what would have been the real effect of the treatment for the patient if they would have taken it
- experiment_number = identifier for the simulated experiment

Outcome is based on the same 'lab-measurements' as the covariates so it is not entirely accurate. True effect is not a measurement, but full actual ground truth.

When analyzing the experiments we act as we would not know the treatment compliance and the true effect of the treatment as in a real setting we would not have this information

## Project Setup

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
