# TODO:
# Create diagram of inputs/outputs and weights
# Research propensity function equations and parameters for biological accuracy
# Write down all equations and justify them - not entirely required
# this model can be quite general and then parameters can be changed to fit to biological systems that you are modelling
# Thoroughly comment code
# RMD report
# Presentation

# output summary statistics from the LHS - such as what combination of parameters is causing the highest fraction of misfolding
# how to find this?? output summary statistics

# show the full dynamics at the start
# then narrow down to a specific question that you are summarising the model for
# show different instances of the same parameter in one plot
# dont have to choose all parameters - choose one parameter that is of interest - fix one parameter e.g. k-folding-mutated and change all of the other parameters on the y axis and see how this affects it
# possibly remove chaperones from the diagram and from plots etc. not interesting

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
  k_chaperone_folding = 0.01,      # Chaperone-assisted folding rate constant
  k_degradation = 0.01,            # Protein degradation rate constant
  feedback_strength = 0.0001       # Feedback strength for negative feedback
)

# more interesting base parameter values
params <- c(
  k_folding_unmutated = 0.3,         # Folding rate constant for proteins from unmutated RNA
  k_folding_mutated = 0.9,           # Folding rate constant for proteins from mutated RNA
  k_misfolding_unmutated = 0.4,      # Misfolding rate constant for proteins from unmutated RNA
  k_misfolding_mutated = 0.9,        # Misfolding rate constant for proteins from mutated RNA
  k_chaperone_folding = 0.02,        # Chaperone-assisted folding rate constant
  k_degradation = 0,              # Protein degradation rate constant
  feedback_strength = 0.0001         # Feedback strength for negative feedback
)


# Initial state of the system
init <- c(
  RNA_unmutated = 1000,            # Number of unmutated RNA molecules present
  RNA_mutated = 1000,               # Number of mutated RNA molecules present
  Protein_success = 0,             # Number of successfully folded proteins
  Protein_misfolded = 0,           # Number of misfolded proteins
  Chaperone = 2,                 # Number of chaperones available
  Degradation_signal = 2          # Initial level of the degradation signal
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

t_final <- 5   # Final time for the simulation
n_steps <- 50  # Number of time steps for storing the output

result <- ssa.adaptivetau(init.values = init, 
                          transitions = transitions, 
                          rateFunc = propensity_function, 
                          params = params, 
                          tf = t_final)
result

result <- as.data.frame(result)


# Plot the time-course data for protein_success and protein_misfolded
ggplot(result, aes(x = time)) +
  geom_line(aes(y = Protein_success, color = "Protein_success")) +
  geom_line(aes(y = Protein_misfolded, color = "Protein_misfolded")) +
  labs(
    title = "Time-course data for protein_success and protein_misfolded",
    x = "Time",
    y = "Concentration",
    color = "Species"
  ) +
  theme_minimal()



# more interesting base parameter values
params <- c(
  k_folding_unmutated = 0.3,         # Folding rate constant for proteins from unmutated RNA
  k_folding_mutated = 0.9,           # Folding rate constant for proteins from mutated RNA
  k_misfolding_unmutated = 0.4,      # Misfolding rate constant for proteins from unmutated RNA
  k_misfolding_mutated = 0.9,        # Misfolding rate constant for proteins from mutated RNA
  k_chaperone_folding = 0.02,        # Chaperone-assisted folding rate constant
  k_degradation = 0,              # Protein degradation rate constant
  feedback_strength = 0.0001         # Feedback strength for negative feedback
)


# test LHS code - unsure if correct approach

param_ranges <- list(
  k_folding_unmutated = c(0.3, 20),
  k_folding_mutated = c(0.9, 30),
  k_misfolding_unmutated = c(0.4, 7),
  k_misfolding_mutated = c(1, 10),
  k_chaperone_folding = c(0.001, 0.05),
  k_degradation = c(0.0001, 0.1),
  feedback_strength = c(0.000001, 0.01)
)



n_samples <- 1000  # Number of samples (adjust as needed)

# Generate normalized LHS samples
lhs_samples <- randomLHS(n_samples, length(param_ranges))

# Scale the samples to the defined parameter ranges
# these param samples are n_samples * the parameter ranges
param_samples <- matrix(NA, nrow = n_samples, ncol = length(param_ranges),
                        dimnames = list(NULL, names(param_ranges)))

for (i in seq_along(param_ranges)) {
  range_i <- param_ranges[[i]]
  param_samples[, i] <- runif(n_samples, min = range_i[1], max = range_i[2])
  }

t_max <- 10

results <- vector("list", length = n_samples)
for (i in 1:n_samples) {
  params_i <- as.list(param_samples[i, ])
  results[[i]] <- ssa.adaptivetau(init, transitions, propensity_function, params_i, t_max)
}


fractions_misfolded <- numeric(length(results))

for (i in seq_along(results)) {
  result <- results[[i]]
  final_state <- result[nrow(result), 4:5] # Get the final state of Protein_success (4th column) and Protein_misfolded (5th column)
  # print(final_state)
  fraction_misfolded <- final_state[2] / sum(final_state) # Calculate the fraction of misfolded proteins
  fractions_misfolded[i] <- fraction_misfolded
}

fractions_misfolded

# Find the index of the simulation with the highest fraction of misfolded proteins
max_index <- which.max(fractions_misfolded)

# Retrieve the parameter combination that corresponds to the index
optimal_params <- param_samples[max_index, ]

# Combine the input parameters and fractions_misfolded into a single data frame
input_parameters <- data.frame(param_samples)
colnames(input_parameters) <- names(param_ranges)
input_parameters$fraction_misfolded <- fractions_misfolded

# Create a scatter plot matrix
pairs(input_parameters)














library(ggplot2)

all_data <- do.call(rbind, lapply(seq_along(results), function(i) {
  result <- results[[i]]
  data_i <- data.frame(time = result[, "time"],
                       RNA_unmutated = result[, "RNA_unmutated"],
                       RNA_mutated = result[, "RNA_mutated"],
                       Protein_success = result[, "Protein_success"],
                       Protein_misfolded = result[, "Protein_misfolded"],
                       Chaperone = result[, "Chaperone"],
                       Degradation_signal = result[, "Degradation_signal"],
                       sample_id = i)
  return(data_i)
}))





# Helper function to plot a single simulation
plot_simulation <- function(data, simulation_id) {
  ggplot(data[data$sample_id == simulation_id,], aes(x = time)) +
    geom_line(aes(y = Protein_success, color = "Protein_success")) +
    geom_line(aes(y = Protein_misfolded, color = "Protein_misfolded")) +
    labs(title = paste("Simulation", simulation_id),
         x = "Time",
         y = "Protein count",
         color = "Species") +
    theme_minimal()
}

# Plot a few simulations (e.g., first 3)
simulation_ids <- unique(all_data$sample_id)[1:3]
plots <- lapply(simulation_ids, function(simulation_id) {
  plot_simulation(all_data, simulation_id)
})

# Display the plots
gridExtra::grid.arrange(grobs = plots, ncol = 1)


