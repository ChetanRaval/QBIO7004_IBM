library(adaptivetau)
library(tidyverse)
library(lhs)

# more interesting base parameter values
params <- c(
  k_folding_unmutated = 0.3,         # Folding rate constant for proteins from unmutated RNA
  k_folding_mutated = 0.5,           # Folding rate constant for proteins from mutated RNA
  k_misfolding_unmutated = 0.4,      # Misfolding rate constant for proteins from unmutated RNA
  k_misfolding_mutated = 0.9,        # Misfolding rate constant for proteins from mutated RNA
  k_chaperone_folding = 0.02,        # Chaperone-assisted folding rate constant
  k_degradation = 0.02,              # Protein degradation rate constant
  feedback_strength = 0.0001         # Feedback strength for negative feedback
)

# Initial state of the system
init <- c(
  RNA_unmutated = 1000,            # Number of unmutated RNA molecules present
  RNA_mutated = 1000,              # Number of mutated RNA molecules present
  Protein_success = 0,             # Number of successfully folded proteins
  Protein_misfolded = 0,           # Number of misfolded proteins
  Chaperone = 15,                  # Number of chaperones available
  Degradation_signal = 5         # Initial level of the degradation signal
)

# Calculate the rate at which each reaction in the system can occur dependent on parameters and initial state
# This function will allow for calculations of rates and probabilities to change over time as we move through time steps
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
  c(RNA_mutated = -1, Protein_misfolded = +1),
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

# Set the maximum number of tries
max_tries <- 10

# Initialize the maximum misfolded fraction
max_misfold <- 0

# Initialize the number of tries
tries <- 0

while (max_misfold <= 0) {
  # Your parameter ranges, number of samples, and LHS sampling code here...
  param_ranges <- list(
    k_folding_unmutated = c(0.3, 20),
    k_folding_mutated = c(0.9, 30),
    k_misfolding_unmutated = c(0.4, 7),
    k_misfolding_mutated = c(1, 10),
    k_chaperone_folding = c(0.001, 0.05),
    k_degradation = c(0.0001, 0.1),
    feedback_strength = c(0.000001, 0.01)
  )
  
  n_samples <- 1000
  lhs_samples <- randomLHS(n_samples, length(param_ranges))
  
  param_samples <- matrix(NA, nrow = n_samples, ncol = length(param_ranges),
                          dimnames = list(NULL, names(param_ranges)))
  
  for (i in seq_along(param_ranges)) {
    range_i <- param_ranges[[i]]
    param_samples[, i] <- runif(n_samples, min = range_i[1], max = range_i[2])
  }
  
  t_max <- 1000
  
  # Your simulation code here...
  results <- vector("list", length = n_samples)
  
  for (i in 1:n_samples) {
    params_i <- as.list(param_samples[i, ])
    results[[i]] <- ssa.adaptivetau(init, transitions, propensity_function, params_i, t_max)
  }
  
  fractions_misfolded <- numeric(length(results))
  
  for (i in seq_along(results)) {
    result <- results[[i]]
    final_state <- result[nrow(result), 4:5] 
    fraction_misfolded <- final_state[2] / sum(final_state)
    fractions_misfolded[i] <- fraction_misfolded
  }
  
  max_index <- which.max(fractions_misfolded)
  max_misfold <- fractions_misfolded[max_index]
  
  print(paste("Try: ", tries, " - Current max misfold:", max_misfold))
  
  tries <- tries + 1
  if (tries > max_tries) {
    print(paste("Stopped after", max_tries, "tries."))
    break
  }
  
  if (max_misfold > 0) {
    optimal_params <- as.list(param_samples[max_index, ])
    df <- data.frame(max_misfold = max_misfold, unlist(optimal_params))
    write_csv(df, "optimal_params.csv")
    print("CSV file written.")
  }
}














