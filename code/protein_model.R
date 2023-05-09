# TODO:
# Create animation
# Create diagram of inputs/outputs and weights
# Research propensity function equations and parameters for biological accuracy
# Write down all equations and justify them
# Thoroughly comment code
# RMD report
# Presentation

library(adaptivetau)
library(tidyverse)
library(gganimate)
library(lhs)

# Model parameters
params <- c(
  k_folding_unmutated = 1,         # Folding rate constant for proteins from unmutated RNA
  k_folding_mutated = 0.9,         # Folding rate constant for proteins from mutated RNA
  k_misfolding_unmutated = 0.1,    # Misfolding rate constant for proteins from unmutated RNA
  k_misfolding_mutated = 0.3,      # Misfolding rate constant for proteins from mutated RNA
  k_chaperone_folding = 0.05,      # Chaperone-assisted folding rate constant
  k_degradation = 0.01,            # Protein degradation rate constant
  feedback_strength = 0.0001       # Feedback strength for negative feedback
)

# Initial state of the system
init <- c(
  RNA_unmutated = 1000,            # Number of unmutated RNA molecules present
  RNA_mutated = 100,               # Number of mutated RNA molecules present
  Protein_success = 0,             # Number of successfully folded proteins
  Protein_misfolded = 0,           # Number of misfolded proteins
  Chaperone = 500,                 # Number of chaperone molecules available
  Degradation_signal = 50          # Initial level of the degradation signal
)


# Propensity function
propensity_function <- function(state, params, t) {
  return(c(
    # Folding reaction for proteins translated from unmutated RNA
    params[["k_folding_unmutated"]] * state[["RNA_unmutated"]],
    # Folding reaction for proteins translated from mutated RNA
    params[["k_folding_mutated"]] * state[["RNA_mutated"]],
    # Misfolding reaction for proteins translated from unmutated RNA
    params[["k_misfolding_unmutated"]] * state[["RNA_unmutated"]],
    # Misfolding reaction for proteins translated from mutated RNA
    params[["k_misfolding_mutated"]] * state[["RNA_mutated"]],
    # Chaperone-assisted folding reaction for misfolded proteins
    params[["k_chaperone_folding"]] * state[["Protein_misfolded"]] * state[["Chaperone"]],
    # Protein degradation reaction for misfolded proteins through lysosome-like process
    params[["k_degradation"]] * state[["Protein_misfolded"]] * state[["Degradation_signal"]],
    # Negative feedback reaction that increases the degradation signal when
    # there is a high concentration of successfully folded proteins
    params[["feedback_strength"]] * state[["Protein_success"]]^2
  ))
}

transitions <- list(
  # Folding reaction for proteins translated from unmutated RNA
  c(RNA_unmutated = -1, Protein_success = +1),
  # Folding reaction for proteins translated from mutated RNA
  c(RNA_mutated = -1, Protein_success = +1),
  # Misfolding reaction for proteins translated from unmutated RNA
  c(RNA_unmutated = -1, Protein_misfolded = +1),
  # Misfolding reaction for proteins translated from mutated RNA
  c(RNA_mutated = -1, Protein_misfolded = +1),
  # Chaperone-assisted folding reaction for misfolded proteins
  c(Protein_misfolded = -1, Protein_success = +1),
  # Protein degradation reaction for misfolded proteins through lysosome-like process
  c(Protein_misfolded = -1),
  # Negative feedback reaction that increases the degradation signal when
  # there is a high concentration of successfully folded proteins
  c(Degradation_signal = +1)
)

t_final <- 1000   # Final time for the simulation
n_steps <- 1000  # Number of time steps for storing the output

result <- ssa.adaptivetau(init.values = init, 
                          transitions = transitions, 
                          rateFunc = propensity_function, 
                          params = params, 
                          tf = 25)
result

result <- as.data.frame(result)


# test LHS code - unsure if correct approach

# Define the range for each parameter to test
param_ranges <- list(
  k_folding_unmutated = c(0.1, 10),
  k_folding_mutated = c(0.1, 10),
  k_misfolding_unmutated = c(0.01, 1),
  k_misfolding_mutated = c(0.01, 1),
  k_chaperone_folding = c(0.01, 1),
  k_degradation = c(0.001, 0.1),
  feedback_strength = c(0.00001, 0.001)
)

n_samples <- 100  # Number of samples (adjust as needed)

# Generate normalized LHS samples
lhs_samples <- randomLHS(n_samples, length(param_ranges))

# Scale the samples to the defined parameter ranges
param_samples <- matrix(NA, nrow = n_samples, ncol = length(param_ranges),
                        dimnames = list(NULL, names(param_ranges)))
for (i in seq_along(param_ranges)) {
  range_i <- param_ranges[[i]]
  param_samples[, i] <- range_i[1] + lhs_samples[, i] * (range_i[2] - range_i[1])
}

t_max <- 100
results <- vector("list", length = n_samples)
for (i in 1:n_samples) {
  params_i <- as.list(param_samples[i, ])
  results[[i]] <- ssa.adaptivetau(init, transitions, propensity_function, params_i, t_max)
}