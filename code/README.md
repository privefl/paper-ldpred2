## Structure of the code

#### Prepare data (for both simulations and real data)

1. `prepare-genotypes`: UKBB Europeans with imputed dosages for HapMap3 variants

2. `prepare-corr`: LD matrices for LDpred2 based on the 10K individuals in the validation set


### Simulations

#### Prepare data

1. `prepare-simu-phenotypes`: Simulation of phenotypes using a liability threshold model

2. `prepare-simu-sumstats`: Logistic GWAS on simulated phenotypes

#### Run methods

- `run-simu-ldpred1` + `process-betas-simu-ldpred1`

- `run-simu-ldpred2` (per chromosome)

- `run-simu-ldpred2-gwide`

#### Visualize results

- `plot-res-simus-ldpred`: Figure 1

- `plot-res-simus2-ldpred`: Figures 2 and S2


### Real data

#### Prepare data

- `prepare-phenotypes`: Define phenotypes in the UK Biobank

- `prepare-sumstats`: Quality control and formatting of external published summary statistics

#### Run methods

- `run-ldpred1` + `process-betas-ldpred1`

- `run-ldpred2` (all LDpred2 models, per chromosome)

- `run-ldpred2-gwide` (all LDpred2 models, genome-wide)

- `run-lassosum` (and lassosum-auto)

- `run-prscs` (and PRS-CS-auto) + `process-betas-prscs`

- `run-sbayesr` + `process-betas-sbayesr`

- `run-sct` (C+T and SCT)

#### Visualize results

- `plot-res-ldpred`: Figures 3 and S3

- `plot-res-all`: Figures 4 and S1


### Provide an LD reference

- `provide-ld-ref`: Make the LD references provided + figure S8

- `example-with-provided-ldref`: Example script of using the provided LD reference for running LDpred2-auto directly


### Misc

- `plot-params`: Figures S4 and S5

- `plot-approx-sd`: Figures S6-S9
