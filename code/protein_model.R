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
  RNA_unmutated = 1000,            # 
  RNA_mutated = 100,
  Protein_success = 0,
  Protein_misfolded = 0,
  Chaperone = 500,
  Degradation_signal = 50
)

# Propensity function
propensity_function <- function(state, params, t) {
  return(c(
    # Folding reaction for proteins translated from unmutated RNA
    params["k_folding_unmutated"] * state["RNA_unmutated"],
    # Folding reaction for proteins translated from mutated RNA
    params["k_folding_mutated"] * state["RNA_mutated"],
    # Misfolding reaction for proteins translated from unmutated RNA
    params["k_misfolding_unmutated"] * state["RNA_unmutated"],
    # Misfolding reaction for proteins translated from mutated RNA
    params["k_misfolding_mutated"] * state["RNA_mutated"],
    # Chaperone-assisted folding reaction for misfolded proteins
    params["k_chaperone_folding"] * state["Protein_misfolded"] * state["Chaperone"],
    # Protein degradation reaction for misfolded proteins through lysosome-like process
    params["k_degradation"] * state["Protein_misfolded"] * state["Degradation_signal"],
    # Negative feedback reaction that increases the degradation signal when
    # there is a high concentration of successfully folded proteins
    params["feedback_strength"] * state["Protein_success"]^2
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