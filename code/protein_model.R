# QBIO7004 - COMPUTATIONAL METHODS IN BIOLOGY
# Stochastic Event-Based Model
# Protein folding/misfolding and aggregation implications in neurodegenerative diseases
# Chetan Raval
# c.raval@uqconnect.edu.au

library(adaptivetau)
library(tidyverse)
library(lhs)
library(gganimate)
library(ppcor)

# load .RData for simulation reproducibility
# load("../output/lhs_simulation.RData")

# if .RData file not available, set the following seed
set.seed(7004)

# Model parameters
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
  Chaperone = 50,                  # Number of chaperones available
  Degradation_signal = 150         # Initial level of the degradation signal
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

# Function to calculate transition rates, given the propensity function and parameters
# init cell state values are updated based on the above propensity functions
# each transition increments the mentioned variables in the same order they appear in propensities function
# e.g. transitions[1] updates after the 1st calculation in propensity function is completed and so on
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

# Simulation time
t_final <- 10   

result <- ssa.adaptivetau(init.values = init, 
                          transitions = transitions, 
                          rateFunc = propensity_function, 
                          params = params, 
                          tf = t_final)

# convert adaptive tau results into dataframe for plotting
result <- as.data.frame(result)
# check head of dataframe
head(result)

# convert result dataframe into long format for ggplot output
result_long <- result %>% 
  pivot_longer(cols = -time, names_to = "variable", values_to = "value")

# plot the dynamics of the model for the above simulation
ggplot(result_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "Time", y = "Value", color = "Variable", 
       title = "Protein Misfolding Model") +
  theme_minimal()


# Plot the time-course data for protein_success and protein_misfolded
ggplot(result, aes(x = time)) +
  geom_line(aes(y = Protein_success, color = "Protein_success")) +
  geom_line(aes(y = Protein_misfolded, color = "Protein_misfolded")) +
  labs(
    title = "Model dynamics for successfully folded and misfolded proteins",
    x = "Time",
    y = "Concentration",
    color = "Species"
  ) +
  theme_minimal()


# -------- Latin Hypercube Sampling -------- #

# create parameter ranges to sample
param_ranges <- list(
  k_folding_unmutated = c(0.3, 20),
  k_folding_mutated = c(0.9, 30),
  k_misfolding_unmutated = c(0.4, 7),
  k_misfolding_mutated = c(1, 10),
  k_chaperone_folding = c(0.001, 0.05),
  k_degradation = c(0.0001, 0.1),
  feedback_strength = c(0.000001, 0.01)
)

# Number of samples
n_samples <- 10000

# Generate normalized LHS samples
lhs_samples <- randomLHS(n_samples, length(param_ranges))

# Scale the samples to the defined parameter ranges
# these param samples are n_samples * the parameter ranges
param_samples <- matrix(NA, nrow = n_samples, ncol = length(param_ranges),
                        dimnames = list(NULL, names(param_ranges)))


# Loop through the indices of the list of parameter ranges
for (i in seq_along(param_ranges)) {
  # Get the range (min and max) for each parameter
  range_i <- param_ranges[[i]]
  # Generate a random sample of size n_samples from a uniform distribution 
  # with the minimum and maximum values specified by range_i.
  # These random numbers represent sampled values of each parameter, 
  # which are stored in the i-th column of the param_samples matrix.
  param_samples[, i] <- runif(n_samples, min = range_i[1], max = range_i[2])
}

# max simulation time
t_max <- 10

# Initialize a list of length n_samples to store the results of each simulation
results <- vector("list", length = n_samples)
# Loop over the number of samples
for (i in 1:n_samples) {
  # For each sample, extract the i-th row from the param_samples matrix,
  # which corresponds to the parameter values for the i-th simulation.
  # Convert this row to a list, because the ssa.adaptivetau function expects its parameters in this format.
  params_i <- as.list(param_samples[i, ])
  # Run the stochastic simulation algorithm with the Adaptive Tau method for the current set of parameters.
  # The initial state of the system is given by "init", the set of reactions by "transitions",
  # the propensity function by "propensity_function", the parameters by "params_i", 
  # and the simulation runs until time "t_max".
  # Store the result of the simulation in the i-th element of the results list.
  results[[i]] <- ssa.adaptivetau(init, transitions, propensity_function, params_i, t_max)
}

# create a list to store misfolded fraction results in
fractions_misfolded <- numeric(length(results))

# 
for (i in seq_along(results)) {
  result <- results[[i]]
  # Get the final state of Protein_success (4th column) and Protein_misfolded (5th column)
  final_state <- result[nrow(result), 4:5] 
  # print(final_state) - testing loop
  fraction_misfolded <- final_state[2] / sum(final_state) # Calculate the fraction of misfolded proteins
  # append to i-th index of empty list
  fractions_misfolded[i] <- fraction_misfolded
}

# Find the index of the simulation with the highest fraction of misfolded proteins
max_index <- which.max(fractions_misfolded)

# Maximum fraction/percentage misfolded value
max_misfold <- fractions_misfolded[max_index]
# show max_misfold list
max_misfold

# Retrieve the parameter combination that corresponds to the index
optimal_params <- as.data.frame(param_samples[max_index, ])
optimal_params$max_misfold <- max_misfold

# Combine the input parameters and fractions_misfolded into a single data frame
input_parameters <- as.data.frame(param_samples)
colnames(input_parameters) <- names(param_ranges)
input_parameters$fraction_misfolded <- fractions_misfolded



# extract the simulation conditions for the maximum misfolded event
simulation_data <- as.data.frame(results[[max_index]])
# restructure data for plotting
simulation_data_long <- simulation_data %>% 
  pivot_longer(cols = -time, names_to = "variable", values_to = "value")
# plot this simulation
ggplot(simulation_data_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "Time", y = "Value", color = "Variable", 
       title = "Misfolding Model - Maximum Misfolded Fraction") +
  theme_minimal()


# export sessionInfo to .txt file for reproducibility
writeLines(capture.output(sessionInfo()), "../sessionInfo.txt")

