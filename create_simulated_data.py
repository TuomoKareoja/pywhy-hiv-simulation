# %%

import whynot as wn
import pandas as pd
import numpy as np
import os

# %%

seed = 123

# Experiment structure parameters

patients = 500
observation_days = 150
treatment_assignment_day = 100

# Covariate and treatment assignment parameters

# NOTE: uniform jitter with (min, max) around the default state
uninfected_t1_jitter = (0.45, 2.15)
infected_t1_jitter = (0.45, 2.15)
uninfected_t2_jitter = (0.45, 2.15)
infected_t2_jitter = (0.45, 2.15)
free_virus_jitter = (0.45, 2.15)
immune_response_jitter = (0.45, 2.15)

study_enrollment_prob = 0.5

# Treatment efficacy parameters
# NOTE: no_treatment_efficacy is the natural course of the disease
# with standard care, while treatment_drug_efficacy is the course of the disease
# with the treatment drug
no_treatment_efficacy = 0.1
treatment_drug_efficacy = 0.5

# imperfect treatment compliance parameters
# NOTE: with imperfect compliance, some of the patients in the treatment group
# will not take the treatment, and some of the patients in the control group
# somehow get access to it

# probability of taking the treatment if enrolled in the study
study_treatment_prob = 0.9
# change to take treatment if not enrolled in the study
control_treatment_prob = 0.05

# Confounder parameters
# NOTE: these parameters are only used if the treatment probability is confounded by the
# severity of the symptoms. They change how likely a patient in the treatment group
# is to take the treatment based on their immune response and free virus levels
confounder_immune_response_threshold = 10
confounder_free_virus_threshold = 1
confounder_treatment_prob_threshold = 0.8
confounder_treatment_prob = 0.2

# %%


def initial_covariate_distribution(rng):
    """Create patients by scrambling the default parameters.

    Parameters
    ----------
        rng: numpy random number generator (
        controlled by the seed parameter).

    Return
    ------
        wn.hiv.State: Initial state of the simulator.
    """
    state = wn.hiv.State()

    state.uninfected_T1 *= rng.uniform(*uninfected_t1_jitter)
    state.infected_T1 *= rng.uniform(*infected_t1_jitter)
    state.uninfected_T2 *= rng.uniform(*uninfected_t2_jitter)
    state.infected_T2 *= rng.uniform(*infected_t2_jitter)
    state.free_virus *= rng.uniform(*free_virus_jitter)
    state.immune_response *= rng.uniform(*immune_response_jitter)

    # Whether or not the unit is "enrolled in the study"
    # NOTE: custom attribute
    state.enrolled = int(rng.rand() < study_enrollment_prob)

    return state


# %%


def create_treatment_propensity(full_compliance, confounding):
    """Create a propensity function with full or imperfect compliance and confounding or no confounding.

    Parameters
    ----------
    full_compliance : bool
        If True, all no patients in the control group take the treatment and all patients
        in the treatment group take the treatment unless there is confounding.

    Returns
    -------
    function
        The propensity scorer function
    """

    def treatment_propensity(intervention, untreated_run):
        run = untreated_run
        if run.initial_state.enrolled > 0:
            # if the there is confounding some of the patients in the treatment group
            # are more likely to take the treatment based on confounders at the time when
            # the treatment starts
            if confounding:
                if (
                    run[intervention.time].immune_response > confounder_immune_response_threshold
                    and run[intervention.time].free_virus > confounder_free_virus_threshold
                ):
                    return confounder_treatment_prob_threshold
                else:
                    return confounder_treatment_prob
            # if there is no confounding and full compliance, all patients in the treatment group
            # take the treatment.
            elif full_compliance:
                return 1.0
            # if there is no confounding and imperfect compliance, not all patients in the treatment group
            # take the treatment.
            else:
                return study_treatment_prob
        return 0.0 if full_compliance else control_treatment_prob

    return treatment_propensity


# %%

# Create four experiments with different settings

experiment_full_compliance = wn.DynamicsExperiment(
    name="hiv_full_compliance",
    description="Full treatment compliance without confounding",
    simulator=wn.hiv,
    simulator_config=wn.hiv.Config(epsilon_1=no_treatment_efficacy, end_time=observation_days),
    intervention=wn.hiv.Intervention(time=treatment_assignment_day, epsilon_1=treatment_drug_efficacy),
    state_sampler=initial_covariate_distribution,
    propensity_scorer=create_treatment_propensity(full_compliance=True, confounding=False),
    outcome_extractor=lambda run: run[observation_days - 1].infected_T2,
    covariate_builder=lambda intervention, run: np.append(
        run[treatment_assignment_day].values(), run.initial_state.enrolled
    ),
)

experiment_partial_compliance = wn.DynamicsExperiment(
    name="hiv_partial_compliance",
    description="Partial treatment compliance without confounding",
    simulator=wn.hiv,
    simulator_config=wn.hiv.Config(epsilon_1=no_treatment_efficacy, end_time=observation_days),
    intervention=wn.hiv.Intervention(time=treatment_assignment_day, epsilon_1=treatment_drug_efficacy),
    state_sampler=initial_covariate_distribution,
    propensity_scorer=create_treatment_propensity(full_compliance=False, confounding=False),
    outcome_extractor=lambda run: run[observation_days - 1].infected_T2,
    covariate_builder=lambda intervention, run: np.append(
        run[treatment_assignment_day].values(), run.initial_state.enrolled
    ),
)

experiment_full_compliance_confounding = wn.DynamicsExperiment(
    name="hiv_full_compliance_confounding",
    description="Full treatment compliance with confounding",
    simulator=wn.hiv,
    simulator_config=wn.hiv.Config(epsilon_1=no_treatment_efficacy, end_time=observation_days),
    intervention=wn.hiv.Intervention(time=treatment_assignment_day, epsilon_1=treatment_drug_efficacy),
    state_sampler=initial_covariate_distribution,
    propensity_scorer=create_treatment_propensity(full_compliance=True, confounding=True),
    outcome_extractor=lambda run: run[observation_days - 1].infected_T2,
    covariate_builder=lambda intervention, run: np.append(
        run[treatment_assignment_day].values(), run.initial_state.enrolled
    ),
)

experiment_partial_compliance_confounding = wn.DynamicsExperiment(
    name="hiv_partial_compliance_confounding",
    description="Partial treatment compliance with confounding",
    simulator=wn.hiv,
    simulator_config=wn.hiv.Config(epsilon_1=no_treatment_efficacy, end_time=observation_days),
    intervention=wn.hiv.Intervention(time=treatment_assignment_day, epsilon_1=treatment_drug_efficacy),
    state_sampler=initial_covariate_distribution,
    propensity_scorer=create_treatment_propensity(full_compliance=False, confounding=True),
    outcome_extractor=lambda run: run[observation_days - 1].infected_T2,
    covariate_builder=lambda intervention, run: np.append(
        run[treatment_assignment_day].values(), run.initial_state.enrolled
    ),
)

# Run all experiments with the same seed
datasets = {
    "hiv_full_compliance": experiment_full_compliance.run(num_samples=patients, seed=seed),
    "hiv_partial_compliance": experiment_partial_compliance.run(num_samples=patients, seed=seed),
    "hiv_full_compliance_confounding": experiment_full_compliance_confounding.run(num_samples=patients, seed=seed),
    "hiv_partial_compliance_confounding": experiment_partial_compliance_confounding.run(
        num_samples=patients, seed=seed
    ),
}

# Process and save all datasets
for dataset_name, dset in datasets.items():
    # Turn into pandas dataframe
    df = pd.DataFrame(
        dset.covariates,
        columns=[
            "uninfected_t1",
            "infected_t1",
            "uninfected_t2",
            "infected_t2",
            "free_virus",
            "immune_response",
            "instrument",
        ],
    )
    df["treatment"] = dset.treatments
    df["outcome"] = dset.outcomes
    df["true_effect"] = dset.true_effects

    # Save to csv
    df.to_csv(os.path.join("data", f"{dataset_name}.csv"), index=False)
    print(f"Created dataset: {dataset_name}.csv")

# %%
