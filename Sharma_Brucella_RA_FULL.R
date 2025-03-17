#####################################################################
####      Annotated Brucella abortus Risk Assessment Model       ####
#####################################################################

# Load libraries
library(ggplot2)
library(epiR)
library(MASS)
library(sensitivity)
library(RColorBrewer)
library(reshape2)
options(scipen = 999) # Avoid scientific notation

# Set seed for reproducibility
set.seed(12345)

# Define output directory (use relative path or allow user to specify)
output_dir <- "output_figures"
dir.create(output_dir, showWarnings = FALSE)

# Number of Monte Carlo iterations
iter <- 10000

#---------------------------------------------------------------------------
#### Step 1: Function to estimate apparent seroprevalence (Beta Distribution) ####

calculate_apparent_seroprevalence <- function(positive, total, label) {
  app_seroprevalence <- rbeta(iter, positive + 1, total - positive + 1)
  
  # Create histogram
  hist(app_seroprevalence, 
       breaks = 'fd', 
       col = 'red3', 
       main = "", 
       xlim = c(0, 0.3), 
       xlab = paste("Apparent seroprevalence -", label), 
       ylim = c(0, 1000))
  
  # Return summary statistics
  return(list(
    seroprevalence = app_seroprevalence,
    summary = summary(app_seroprevalence),
    ci95 = quantile(app_seroprevalence, c(0.025, 0.975))
  ))
}

# Calculate seroprevalence for different populations
results_overall <- calculate_apparent_seroprevalence(41, 261, "overall")
results_cattle <- calculate_apparent_seroprevalence(33, 181, "cattle")
results_buffalo <- calculate_apparent_seroprevalence(8, 80, "buffalo")

App_Seroprev <- results_overall$seroprevalence
App_SeroprevC <- results_cattle$seroprevalence
App_SeroprevB <- results_buffalo$seroprevalence

#---------------------------------------------------------------------------
#### Step 2: Adjust for Sensitivity & Specificity ####

# Diagnostic test characteristics
Se1 <- 0.874 
Se2 <- 0.846 
Sp1 <- 0.994 
Sp2 <- 0.996 

# Serial Testing (Higher Specificity)
Se_series <- Se1 * Se2
Sp_series <- 1 - (1 - Sp1) * (1 - Sp2)

# Function to compute true prevalence
calculate_true_prevalence <- function(apparent_prevalence, sensitivity, specificity) {
  true_prevalence <- (apparent_prevalence + specificity - 1) / (sensitivity + specificity - 1)
  return(true_prevalence)
}

# Compute True Prevalence for each population
TPrevO <- calculate_true_prevalence(App_Seroprev, Se_series, Sp_series)
TPrevC <- calculate_true_prevalence(App_SeroprevC, Se_series, Sp_series)
TPrevB <- calculate_true_prevalence(App_SeroprevB, Se_series, Sp_series)

# Create boxplot of true prevalence values
png(file.path(output_dir, "TruePrev.png"), width = 5, height = 4, units = 'in', res = 300)
par(mfrow = c(1, 1), mar = c(3, 4, 0.3, 0.2))
boxplot(TPrevC, TPrevB, TPrevO, 
        names = c("Cattle", "Buffalo", "Combined"), 
        ylab = "True prevalence")
dev.off()

#---------------------------------------------------------------------------
#### Step 3: Probability of Shedding by a seropositive animal ####

# Data from studies
shedders_iran <- 6      # Cows in Iran
seropos_iran <- 260
shedders_italy <- 101   # Buffalo in Italy
seropos_italy <- 337

# Calculate shedding prevalence distributions
calculate_shedding_prevalence <- function(shedders, seropos, label) {
  shed_prev <- rbeta(iter, shedders + 1, seropos - shedders + 1)
  return(list(
    prevalence = shed_prev,
    summary = summary(shed_prev),
    ci95 = quantile(shed_prev, c(0.025, 0.975)),
    label = label
  ))
}

# Most conservative (buffalo in Italy)
shed_results_italy <- calculate_shedding_prevalence(shedders_italy, seropos_italy, "Italy (buffalo)")
Shed_PrevItaly <- shed_results_italy$prevalence

# Least conservative (cows in Iran)
shed_results_iran <- calculate_shedding_prevalence(shedders_iran, seropos_iran, "Iran (cattle)")
Shed_PrevIran <- shed_results_iran$prevalence

# Plot combined histogram
png(file.path(output_dir, "Shedding_Distribution.png"), width = 5, height = 4, units = 'in', res = 300)
par(mfrow = c(1, 1), mar = c(3, 4, 0.3, 0.2))
hist(Shed_PrevItaly, breaks = 50, col = 'grey50', main = "", 
     xlim = c(0, 0.5), ylim = c(0, 1000), 
     xlab = "Shedding proportion (of seropositive animals)")
hist(add = TRUE, Shed_PrevIran, breaks = 50, col = 'grey70')
legend("topright", legend = c("Buffalo (Italy)", "Cattle (Iran)"), 
       fill = c("grey50", "grey70"))
dev.off()

#---------------------------------------------------------------------------
#### Step 4: Combine prevalence and shedding to get overall shedding prevalence ####

# Calculate shedding prevalence for each population
TPrevShedB <- Shed_PrevItaly * TPrevB
TPrevShedC <- Shed_PrevIran * TPrevC

# Calculate overall shedding prevalence based on proportion of cattle and buffalo
PropCows <- 181/(181+80)  # Based on study data

# Function to calculate combined prevalence
calculate_combined_prevalence <- function(prop_cows, prev_cattle, prev_buffalo, n_simulations = 1000) {
  results <- numeric(iter)
  
  for (i in 1:iter) {
    # Sample population
    count_cattle <- rbinom(1, n_simulations, prop_cows)
    count_buff <- n_simulations - count_cattle
    
    # Calculate shedders in each group
    shed_count_cattle <- round(count_cattle * prev_cattle[i], 0)
    shed_count_buff <- round(count_buff * prev_buffalo[i], 0)
    
    # Calculate overall proportion
    results[i] <- (shed_count_cattle + shed_count_buff) / n_simulations
  }
  
  return(results)
}

PropOverallShed <- calculate_combined_prevalence(PropCows, TPrevShedC, TPrevShedB)

# Create boxplot comparing shedding prevalence
png(file.path(output_dir, "Shedding.png"), width = 5, height = 4, units = 'in', res = 300)
par(mfrow = c(1, 1), mar = c(3, 4, 0.3, 0.2))
boxplot(TPrevShedC, TPrevShedB, PropOverallShed, 
        names = c("Cattle", "Buffalo", "Combined"), 
        ylab = "Shedding prevalence")
dev.off()

#---------------------------------------------------------------------------
#### Step 5: Brucella CFU/mL in Milk ####

# Distribution of shedding levels
# 73 Low shedder 30-999 CFU/mL
# 16 high shedder 10000-40000 CFU/mL
# 12 intermittent low shedder

# Calculate probabilities of low vs high shedding
LShed <- 73 + 12
HShed <- 16
TotShed <- LShed + HShed

# Generate distributions
PLowShed <- rbeta(iter, LShed + 1, TotShed - LShed + 1)
PHighShed <- 1 - PLowShed

# Generate CFU distributions
RLowShed <- runif(iter, 0, 999)
RHighShed <- runif(iter, 10000, 40000)

# Function to calculate CFU counts per mL of infected milk
generate_cfu_counts <- function(n_samples, prob_low_shed, low_shed_range, high_shed_range) {
  cfu_counts <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    p_low_shed <- rbeta(1, LShed + 1, TotShed - LShed + 1)
    is_low_shedder <- runif(1) < p_low_shed
    
    if (is_low_shedder) {
      cfu_counts[i] <- sample(low_shed_range, 1)
    } else {
      cfu_counts[i] <- sample(high_shed_range, 1)
    }
  }
  
  return(cfu_counts)
}

# Generate CFU counts
CFU_counts <- generate_cfu_counts(iter, PLowShed, RLowShed, RHighShed)

# Plot distribution of CFU counts
png(file.path(output_dir, "CFU_Distribution.png"), width = 5, height = 4, units = 'in', res = 300)
hist(CFU_counts, col = 'yellow', main = '', breaks = 100, 
     xlab = 'CFU/mL infected milk')
dev.off()

#---------------------------------------------------------------------------
#### Step 6: Calculate CFU in 100mL of bulk milk ####

# Calculate percentage of milk that contains Brucella (assuming 10% reduction from animal prevalence)
calculate_milk_infection_percent <- function(shedding_prevalence) {
  total_milk_shedding <- shedding_prevalence * 0.9
  total_milk_non_shed <- 1 - shedding_prevalence
  percent_shedding <- total_milk_shedding / (total_milk_non_shed + total_milk_shedding)
  return(percent_shedding)
}

# Calculate for each population
PercentSheddingCattle <- calculate_milk_infection_percent(TPrevShedC)
PercentSheddingBuff <- calculate_milk_infection_percent(TPrevShedB)
PercentSheddingComb <- calculate_milk_infection_percent(PropOverallShed)

# Function to simulate CFU in 100mL of milk
simulate_cfu_in_bulk_milk <- function(percent_shedding, cfu_distribution, volume_ml = 100) {
  results <- numeric(length(percent_shedding))
  
  for (i in 1:length(percent_shedding)) {
    infected_ml <- round(percent_shedding[i] * volume_ml, 0)
    if (infected_ml > 0) {
      results[i] <- sum(sample(cfu_distribution, infected_ml, replace = TRUE))
    } else {
      results[i] <- 0
    }
  }
  
  return(results)
}

# Simulate CFU in 100mL for each population
CFU_countCattle100 <- simulate_cfu_in_bulk_milk(PercentSheddingCattle, CFU_counts)
CFU_countBuff100 <- simulate_cfu_in_bulk_milk(PercentSheddingBuff, CFU_counts)
CFU_countComb100 <- simulate_cfu_in_bulk_milk(PercentSheddingComb, CFU_counts)

# Calculate percentage above safety threshold (1000 CFU)
percent_cattle <- mean(CFU_countCattle100 > 1000) * 100
percent_buffalo <- mean(CFU_countBuff100 > 1000) * 100
percent_combined <- mean(CFU_countComb100 > 1000) * 100

# Create summary table
cfu_summary <- data.frame(
  Population = c("Cattle", "Buffalo", "Combined"),
  Mean_CFU = c(mean(CFU_countCattle100), mean(CFU_countBuff100), mean(CFU_countComb100)),
  Median_CFU = c(median(CFU_countCattle100), median(CFU_countBuff100), median(CFU_countComb100)),
  Percent_Above_Threshold = c(percent_cattle, percent_buffalo, percent_combined)
)

print(cfu_summary)

# Create boxplot
png(file.path(output_dir, "CFU100ml.png"), width = 5, height = 4, units = 'in', res = 300)
par(mfrow = c(1, 1), mar = c(3, 4, 1, 0.3))
boxplot(CFU_countCattle100, CFU_countBuff100, CFU_countComb100, 
        names = c("Cattle", "Buffalo", "Combined"), 
        ylab = "CFU/100mL")
abline(h = 1000, lty = 2, col = 2, lwd = 2)
text(2.5, 1200, "Safety threshold", col = 2)
dev.off()

#---------------------------------------------------------------------------
#### Step 7: Dose-Response Model ####

# Infectious dose (ID50)
ID50 <- 10^4

# Calculate infection probability
calculate_infection_probability <- function(cfu_count, id50) {
  1 - exp(-cfu_count / id50)
}

# Calculate for each population
Prob_Infection_Cattle <- calculate_infection_probability(CFU_countCattle100, ID50)
Prob_Infection_Buffalo <- calculate_infection_probability(CFU_countBuff100, ID50)
Prob_Infection_Combined <- calculate_infection_probability(CFU_countComb100, ID50)

# Create summary table
infection_summary <- data.frame(
  Population = c("Cattle", "Buffalo", "Combined"),
  Mean_Probability = c(mean(Prob_Infection_Cattle), mean(Prob_Infection_Buffalo), mean(Prob_Infection_Combined)),
  Median_Probability = c(median(Prob_Infection_Cattle), median(Prob_Infection_Buffalo), median(Prob_Infection_Combined)),
  CI_Lower = c(quantile(Prob_Infection_Cattle, 0.025), quantile(Prob_Infection_Buffalo, 0.025), quantile(Prob_Infection_Combined, 0.025)),
  CI_Upper = c(quantile(Prob_Infection_Cattle, 0.975), quantile(Prob_Infection_Buffalo, 0.975), quantile(Prob_Infection_Combined, 0.975))
)

print(infection_summary)

# Create visualization
png(file.path(output_dir, "InfectionProbability.png"), width = 5, height = 4, units = 'in', res = 300)
boxplot(Prob_Infection_Cattle, Prob_Infection_Buffalo, Prob_Infection_Combined,
        names = c("Cattle", "Buffalo", "Combined"),
        ylab = "Probability of Infection per 100mL Consumption")
dev.off()

#---------------------------------------------------------------------------
#### Step 8: Sensitivity Analysis (Sobol method) ####

# Define improved Sobol function
improved_sobol_function <- function(X) {
  # Calculate shedding prevalence
  tprev_sh_c <- X$PropSheddingC * X$Seroprev_cattle
  tprev_sh_b <- X$PropSheddingB * X$Seroprev_buffalo
  
  # Calculate overall prevalence
  prop_cows <- X$PropCows / 1000  # Convert to proportion
  count_cattle_s <- prop_cows * tprev_sh_c 
  count_buffalo_s <- (1 - prop_cows) * tprev_sh_b
  prop_overall_shed <- count_cattle_s + count_buffalo_s
  
  # Calculate CFU in milk
  prob_low_shedding <- X$ProbLowShedding
  cfu_shedders <- (prob_low_shedding * X$AmountLowShed) + 
    ((1 - prob_low_shedding) * X$AmountHighShed)
  
  # Final CFU concentration
  cfu_all <- cfu_shedders * prop_overall_shed
  return(cfu_all)
}

# Create parameter samples for sensitivity analysis
create_parameter_samples <- function(n) {
  data.frame(
    Seroprev_cattle = sample(TPrevC, n, replace = TRUE),
    Seroprev_buffalo = sample(TPrevB, n, replace = TRUE),
    PropSheddingB = sample(Shed_PrevItaly, n, replace = TRUE),
    PropSheddingC = sample(Shed_PrevIran, n, replace = TRUE),
    PropCows = sample(seq(0, 1000, by = 1), n, replace = TRUE),  # Values from 0-1000
    ProbLowShedding = sample(PLowShed, n, replace = TRUE),
    AmountLowShed = sample(RLowShed, n, replace = TRUE),
    AmountHighShed = sample(RHighShed, n, replace = TRUE)
  )
}

# Run sensitivity analysis
n_samples <- 1000  # Reduced for computation time, increase for final analysis
X1 <- create_parameter_samples(n_samples)
X2 <- create_parameter_samples(n_samples)

# Perform Sobol analysis
sobol_results <- sobol2007(model = improved_sobol_function, X1 = X1, X2 = X2, nboot = 100)

# Create improved visualization with ggplot2
sobol_data <- data.frame(
  Parameter = rownames(sobol_results$S),
  MainEffect = sobol_results$S[, "original"],
  TotalEffect = sobol_results$T[, "original"]
)
rownames(sobol_data) <- NULL

# Define better parameter labels
custom_labels <- c(
  "Seroprev_cattle" = "Seroprevalence cattle",
  "Seroprev_buffalo" = "Seroprevalence buffalo",
  "PropSheddingB" = "Shedding proportion buffalo",
  "PropSheddingC" = "Shedding proportion cattle",
  "PropCows" = "Proportion of cattle",
  "ProbLowShedding" = "Probability of low shedding",
  "AmountLowShed" = "Amount of low shedding",
  "AmountHighShed" = "Amount of high shedding"
)

# Reshape data for plotting
sobol_data_long <- reshape2::melt(sobol_data, id.vars = "Parameter", 
                                  variable.name = "EffectType", 
                                  value.name = "Value")

# Create plot
png(file.path(output_dir, "Sobol.png"), width = 6, height = 5, units = 'in', res = 300)
ggplot(sobol_data_long, aes(x = Parameter, y = Value, color = EffectType)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels = custom_labels) +
  scale_color_manual(values = c("MainEffect" = "#E41A1C", "TotalEffect" = "#377EB8"),
                     labels = c("MainEffect" = "Main Effect", "TotalEffect" = "Total Effect")) +
  labs(title = "",
       x = "Parameter", 
       y = "Sensitivity Index", 
       color = "Effect Type") +
  theme_minimal() +
  ylim(0, 0.8) +
  geom_hline(yintercept = c(0.25, 0.5), linetype = "dashed", color = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 6, r = 6, b = 6, l = 20)
  )
dev.off()

# Save results for further analysis
save(TPrevC, TPrevB, Shed_PrevItaly, Shed_PrevIran, CFU_counts, 
     CFU_countCattle100, CFU_countBuff100, CFU_countComb100,
     Prob_Infection_Cattle, Prob_Infection_Buffalo, Prob_Infection_Combined,
     sobol_results, 
     file = file.path(output_dir, "brucella_model_results.RData"))

# Generate summary report
cat("=== Brucella abortus Risk Assessment Model Results ===\n\n")
cat("True Seroprevalence:\n")
cat("  Cattle: ", mean(TPrevC), " (95% CI: ", quantile(TPrevC, 0.025), "-", quantile(TPrevC, 0.975), ")\n")
cat("  Buffalo: ", mean(TPrevB), " (95% CI: ", quantile(TPrevB, 0.025), "-", quantile(TPrevB, 0.975), ")\n\n")

cat("Shedding Prevalence:\n")
cat("  Cattle: ", mean(TPrevShedC), " (95% CI: ", quantile(TPrevShedC, 0.025), "-", quantile(TPrevShedC, 0.975), ")\n")
cat("  Buffalo: ", mean(TPrevShedB), " (95% CI: ", quantile(TPrevShedB, 0.025), "-", quantile(TPrevShedB, 0.975), ")\n")
cat("  Combined: ", mean(PropOverallShed), " (95% CI: ", quantile(PropOverallShed, 0.025), "-", quantile(PropOverallShed, 0.975), ")\n\n")

cat("CFU in 100mL milk:\n")
cat("  Cattle: ", mean(CFU_countCattle100), " (95% CI: ", quantile(CFU_countCattle100, 0.025), "-", quantile(CFU_countCattle100, 0.975), ")\n")
cat("  Buffalo: ", mean(CFU_countBuff100), " (95% CI: ", quantile(CFU_countBuff100, 0.025), "-", quantile(CFU_countBuff100, 0.975), ")\n")
cat("  Combined: ", mean(CFU_countComb100), " (95% CI: ", quantile(CFU_countComb100, 0.025), "-", quantile(CFU_countComb100, 0.975), ")\n\n")

cat("Probability of infection per 100mL consumption:\n")
cat("  Cattle: ", mean(Prob_Infection_Cattle), " (95% CI: ", quantile(Prob_Infection_Cattle, 0.025), "-", quantile(Prob_Infection_Cattle, 0.975), ")\n")
cat("  Buffalo: ", mean(Prob_Infection_Buffalo), " (95% CI: ", quantile(Prob_Infection_Buffalo, 0.025), "-", quantile(Prob_Infection_Buffalo, 0.975), ")\n")
cat("  Combined: ", mean(Prob_Infection_Combined), " (95% CI: ", quantile(Prob_Infection_Combined, 0.025), "-", quantile(Prob_Infection_Combined, 0.975), ")\n\n")

cat("Most influential parameters (from sensitivity analysis):\n")
print(sobol_data[order(-sobol_data$TotalEffect), ][1:3, ])


#---------------------------------------------------------------------------
#### Step 9: Scenario Analysis - No High Shedders ####

# This section calculates the probability of infection when high shedders are removed,
# which could represent a scenario where management practices identify and remove high shedding animals

# Function to simulate CFU counts with only low shedders
generate_low_shedder_cfu_counts <- function(n_samples, low_shed_range) {
  # All animals are low shedders in this scenario
  sample(low_shed_range, n_samples, replace = TRUE)
}

# Generate CFU counts for low shedders only
CFU_counts_low_only <- generate_low_shedder_cfu_counts(iter, RLowShed)

# Plot distribution of CFU counts for comparison
png(file.path(output_dir, "CFU_Distribution_LowOnly.png"), width = 6, height = 4, units = 'in', res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
hist(CFU_counts, col = 'orange', main = 'All Shedders', breaks = 100, 
     xlab = 'CFU/mL infected milk', xlim = c(0, 40000))
hist(CFU_counts_low_only, col = 'lightgreen', main = 'Low Shedders Only', breaks = 50, 
     xlab = 'CFU/mL infected milk', xlim = c(0, 1000))
dev.off()

# Recalculate CFU in 100mL milk with only low shedders
CFU_countCattle100_low <- simulate_cfu_in_bulk_milk(PercentSheddingCattle, CFU_counts_low_only)
CFU_countBuff100_low <- simulate_cfu_in_bulk_milk(PercentSheddingBuff, CFU_counts_low_only)
CFU_countComb100_low <- simulate_cfu_in_bulk_milk(PercentSheddingComb, CFU_counts_low_only)

# Calculate infection probabilities with only low shedders
Prob_Infection_Cattle_low <- calculate_infection_probability(CFU_countCattle100_low, ID50)
Prob_Infection_Buffalo_low <- calculate_infection_probability(CFU_countBuff100_low, ID50)
Prob_Infection_Combined_low <- calculate_infection_probability(CFU_countComb100_low, ID50)

# Create comparison table
infection_comparison <- data.frame(
  Population = rep(c("Cattle", "Buffalo", "Combined"), 2),
  Scenario = c(rep("With High Shedders", 3), rep("Low Shedders Only", 3)),
  Mean_Probability = c(
    mean(Prob_Infection_Cattle), mean(Prob_Infection_Buffalo), mean(Prob_Infection_Combined),
    mean(Prob_Infection_Cattle_low), mean(Prob_Infection_Buffalo_low), mean(Prob_Infection_Combined_low)
  ),
  Median_Probability = c(
    median(Prob_Infection_Cattle), median(Prob_Infection_Buffalo), median(Prob_Infection_Combined),
    median(Prob_Infection_Cattle_low), median(Prob_Infection_Buffalo_low), median(Prob_Infection_Combined_low)
  ),
  CI_Lower = c(
    quantile(Prob_Infection_Cattle, 0.025), quantile(Prob_Infection_Buffalo, 0.025), quantile(Prob_Infection_Combined, 0.025),
    quantile(Prob_Infection_Cattle_low, 0.025), quantile(Prob_Infection_Buffalo_low, 0.025), quantile(Prob_Infection_Combined_low, 0.025)
  ),
  CI_Upper = c(
    quantile(Prob_Infection_Cattle, 0.975), quantile(Prob_Infection_Buffalo, 0.975), quantile(Prob_Infection_Combined, 0.975),
    quantile(Prob_Infection_Cattle_low, 0.975), quantile(Prob_Infection_Buffalo_low, 0.975), quantile(Prob_Infection_Combined_low, 0.975)
  ),
  Risk_Reduction_Percent = c(
    NA, NA, NA,
    (1 - mean(Prob_Infection_Cattle_low)/mean(Prob_Infection_Cattle)) * 100,
    (1 - mean(Prob_Infection_Buffalo_low)/mean(Prob_Infection_Buffalo)) * 100,
    (1 - mean(Prob_Infection_Combined_low)/mean(Prob_Infection_Combined)) * 100
  )
)

print(infection_comparison)

# Create visualization comparing both scenarios
# Reshape data for plotting
plot_data <- data.frame(
  Population = rep(c("Cattle", "Buffalo", "Combined"), each = 2),
  Scenario = rep(c("With High Shedders", "Low Shedders Only"), 3),
  Probability = c(
    Prob_Infection_Cattle, Prob_Infection_Cattle_low,
    Prob_Infection_Buffalo, Prob_Infection_Buffalo_low,
    Prob_Infection_Combined, Prob_Infection_Combined_low
  )
)

# Calculate risk reduction for a typical meal
# Assuming a typical meal contains 250mL of milk
meal_size_ml <- 250
# Function to calculate meal risk
calculate_meal_risk <- function(prob_per_100ml, meal_size) {
  # Scale probability to reflect the meal size
  scaling_factor <- meal_size / 100
  # Calculate probability for the full meal (compound risk)
  1 - (1 - prob_per_100ml)^scaling_factor
}

# Calculate meal risk for both scenarios
Meal_Risk_Cattle <- calculate_meal_risk(Prob_Infection_Cattle, meal_size_ml)
Meal_Risk_Buffalo <- calculate_meal_risk(Prob_Infection_Buffalo, meal_size_ml)
Meal_Risk_Combined <- calculate_meal_risk(Prob_Infection_Combined, meal_size_ml)

Meal_Risk_Cattle_low <- calculate_meal_risk(Prob_Infection_Cattle_low, meal_size_ml)
Meal_Risk_Buffalo_low <- calculate_meal_risk(Prob_Infection_Buffalo_low, meal_size_ml)
Meal_Risk_Combined_low <- calculate_meal_risk(Prob_Infection_Combined_low, meal_size_ml)

# Create meal risk comparison plot
# Define colors
shedding_colors <- c("pink3", "green4")

png(file.path(output_dir, "Meal_Risk_Comparison.png"), width = 7, height = 5, units = 'in', res = 300)
par(mfrow = c(1, 1), mar = c(3, 4, 1, 2))
# Risk/100ml
boxplot(Prob_Infection_Cattle, Prob_Infection_Cattle_low, 
        Prob_Infection_Buffalo, Prob_Infection_Buffalo_low,
        Prob_Infection_Combined, Prob_Infection_Combined_low,
        col = rep(shedding_colors, 3),
        main = "",
        ylab = "Probability of Infection")
axis(1, at = c(1.5, 3.5, 5.5), labels = c("Cattle", "Buffalo", "Combined"), tick = FALSE, line = 1)
legend("topright", 
       legend = c("With High Shedders", "Low Shedders Only"),
       fill = shedding_colors,
       title = "Shedding",
       cex = 0.8,
       bty = "n")
dev.off()

# Create bar chart showing risk reduction percentages
risk_reduction <- data.frame(
  Population = c("Cattle", "Buffalo", "Combined"),
  Reduction_100ml = c(
    (1 - mean(Prob_Infection_Cattle_low)/mean(Prob_Infection_Cattle)) * 100,
    (1 - mean(Prob_Infection_Buffalo_low)/mean(Prob_Infection_Buffalo)) * 100,
    (1 - mean(Prob_Infection_Combined_low)/mean(Prob_Infection_Combined)) * 100
  ),
  Reduction_Meal = c(
    (1 - mean(Meal_Risk_Cattle_low)/mean(Meal_Risk_Cattle)) * 100,
    (1 - mean(Meal_Risk_Buffalo_low)/mean(Meal_Risk_Buffalo)) * 100,
    (1 - mean(Meal_Risk_Combined_low)/mean(Meal_Risk_Combined)) * 100
  )
)

# Plot risk reduction percentages
png(file.path(output_dir, "Risk_Reduction.png"), width = 6, height = 5, units = 'in', res = 300)
par(mfrow = c(1, 1), mar = c(3, 4, 1, 0.3))
barplot(t(risk_reduction[, 2:3]), beside = TRUE, 
        names.arg = risk_reduction$Population,
        col = c("steelblue", "darkgreen"),
        main = "Risk Reduction by Removing High Shedders",
        ylab = "Percent Reduction (%)",
        ylim = c(0, 100))
legend("topright", legend = c("Per 100mL", paste("Per Meal (", meal_size_ml, "mL)", sep="")), 
       fill = c("steelblue", "darkgreen"))
dev.off()

# Calculate risk reduction for a typical meal
# Assuming a typical meal contains 250mL of milk
meal_size_ml <- 100
# Function to calculate meal risk
calculate_meal_risk <- function(prob_per_100ml, meal_size) {
  # Scale probability to reflect the meal size
  scaling_factor <- meal_size / 100
  # Calculate probability for the full meal (compound risk)
  1 - (1 - prob_per_100ml)^scaling_factor
}

# Calculate meal risk for both scenarios
Meal_Risk_Cattle <- calculate_meal_risk(Prob_Infection_Cattle, meal_size_ml)
Meal_Risk_Buffalo <- calculate_meal_risk(Prob_Infection_Buffalo, meal_size_ml)
Meal_Risk_Combined <- calculate_meal_risk(Prob_Infection_Combined, meal_size_ml)

Meal_Risk_Cattle_low <- calculate_meal_risk(Prob_Infection_Cattle_low, meal_size_ml)
Meal_Risk_Buffalo_low <- calculate_meal_risk(Prob_Infection_Buffalo_low, meal_size_ml)
Meal_Risk_Combined_low <- calculate_meal_risk(Prob_Infection_Combined_low, meal_size_ml)

# Create summary table of meal risk
meal_risk_summary <- data.frame(
  Population = rep(c("Cattle", "Buffalo", "Combined"), 2),
  Scenario = c(rep("With High Shedders", 3), rep("Low Shedders Only", 3)),
  Mean_Meal_Risk = c(
    mean(Meal_Risk_Cattle), mean(Meal_Risk_Buffalo), mean(Meal_Risk_Combined),
    mean(Meal_Risk_Cattle_low), mean(Meal_Risk_Buffalo_low), mean(Meal_Risk_Combined_low)
  ),
  Median_Meal_Risk = c(
    median(Meal_Risk_Cattle), median(Meal_Risk_Buffalo), median(Meal_Risk_Combined),
    median(Meal_Risk_Cattle_low), median(Meal_Risk_Buffalo_low), median(Meal_Risk_Combined_low)
  ),
  CI_Lower = c(
    quantile(Meal_Risk_Cattle, 0.025), quantile(Meal_Risk_Buffalo, 0.025), quantile(Meal_Risk_Combined, 0.025),
    quantile(Meal_Risk_Cattle_low, 0.025), quantile(Meal_Risk_Buffalo_low, 0.025), quantile(Meal_Risk_Combined_low, 0.025)
  ),
  CI_Upper = c(
    quantile(Meal_Risk_Cattle, 0.975), quantile(Meal_Risk_Buffalo, 0.975), quantile(Meal_Risk_Combined, 0.975),
    quantile(Meal_Risk_Cattle_low, 0.975), quantile(Meal_Risk_Buffalo_low, 0.975), quantile(Meal_Risk_Combined_low, 0.975)
  ),
  Risk_Reduction_Percent = c(
    NA, NA, NA,
    (1 - mean(Meal_Risk_Cattle_low)/mean(Meal_Risk_Cattle)) * 100,
    (1 - mean(Meal_Risk_Buffalo_low)/mean(Meal_Risk_Buffalo)) * 100,
    (1 - mean(Meal_Risk_Combined_low)/mean(Meal_Risk_Combined)) * 100
  )
)

# Format the table for better readability
formatted_meal_risk <- meal_risk_summary
formatted_meal_risk$Mean_Meal_Risk <- sprintf("%.4f", formatted_meal_risk$Mean_Meal_Risk)
formatted_meal_risk$Median_Meal_Risk <- sprintf("%.4f", formatted_meal_risk$Median_Meal_Risk)
formatted_meal_risk$CI_Lower <- sprintf("%.4f", formatted_meal_risk$CI_Lower)
formatted_meal_risk$CI_Upper <- sprintf("%.4f", formatted_meal_risk$CI_Upper)
formatted_meal_risk$Risk_Reduction_Percent <- ifelse(is.na(formatted_meal_risk$Risk_Reduction_Percent), 
                                                     NA, 
                                                     sprintf("%.1f%%", formatted_meal_risk$Risk_Reduction_Percent))

# Create a more readable version of the table with CI combined
meal_risk_table <- data.frame(
  Population = formatted_meal_risk$Population,
  Scenario = formatted_meal_risk$Scenario,
  Median_Meal_Risk = formatted_meal_risk$Median_Meal_Risk,
  CI_95 = paste0("(", formatted_meal_risk$CI_Lower, "-", formatted_meal_risk$CI_Upper, ")"),
  Risk_Reduction = formatted_meal_risk$Risk_Reduction_Percent
)

# Display the table
print(paste("Infection Risk from a Typical Meal (", meal_size_ml, "mL)"))
print(meal_risk_table)

# Save the table to CSV
write.csv(meal_risk_table, file = file.path(output_dir, "meal_risk_table.csv"), row.names = FALSE)

#---------------------------------------------------------------------------
#### Step 10: Annual Risk Assessment - Monthly Consumption ####

# Calculate annual infection risk based on monthly consumption of 100ml milk (12 times per year)
consumption_frequency <- 12  # Number of consumption events per year
consumption_volume <- 100    # Volume of milk consumed each time (ml)

# Function to calculate annual risk from periodic consumption
calculate_annual_risk <- function(risk_per_event, num_events) {
  # Formula: P(at least one infection in n events) = 1 - (1 - p)^n
  # where p is the probability of infection per event and n is the number of events
  annual_risk <- 1 - (1 - risk_per_event)^num_events
  return(annual_risk)
}

# Calculate annual risks with high shedders
Annual_Risk_Cattle <- calculate_annual_risk(Prob_Infection_Cattle, consumption_frequency)
Annual_Risk_Buffalo <- calculate_annual_risk(Prob_Infection_Buffalo, consumption_frequency)
Annual_Risk_Combined <- calculate_annual_risk(Prob_Infection_Combined, consumption_frequency)

# Calculate annual risks without high shedders
Annual_Risk_Cattle_low <- calculate_annual_risk(Prob_Infection_Cattle_low, consumption_frequency)
Annual_Risk_Buffalo_low <- calculate_annual_risk(Prob_Infection_Buffalo_low, consumption_frequency)
Annual_Risk_Combined_low <- calculate_annual_risk(Prob_Infection_Combined_low, consumption_frequency)

# Create summary table of annual risk
annual_risk_summary <- data.frame(
  Population = rep(c("Cattle", "Buffalo", "Combined"), 2),
  Scenario = c(rep("With High Shedders", 3), rep("Low Shedders Only", 3)),
  Mean_Annual_Risk = c(
    mean(Annual_Risk_Cattle), mean(Annual_Risk_Buffalo), mean(Annual_Risk_Combined),
    mean(Annual_Risk_Cattle_low), mean(Annual_Risk_Buffalo_low), mean(Annual_Risk_Combined_low)
  ),
  Median_Annual_Risk = c(
    median(Annual_Risk_Cattle), median(Annual_Risk_Buffalo), median(Annual_Risk_Combined),
    median(Annual_Risk_Cattle_low), median(Annual_Risk_Buffalo_low), median(Annual_Risk_Combined_low)
  ),
  CI_Lower = c(
    quantile(Annual_Risk_Cattle, 0.025), quantile(Annual_Risk_Buffalo, 0.025), quantile(Annual_Risk_Combined, 0.025),
    quantile(Annual_Risk_Cattle_low, 0.025), quantile(Annual_Risk_Buffalo_low, 0.025), quantile(Annual_Risk_Combined_low, 0.025)
  ),
  CI_Upper = c(
    quantile(Annual_Risk_Cattle, 0.975), quantile(Annual_Risk_Buffalo, 0.975), quantile(Annual_Risk_Combined, 0.975),
    quantile(Annual_Risk_Cattle_low, 0.975), quantile(Annual_Risk_Buffalo_low, 0.975), quantile(Annual_Risk_Combined_low, 0.975)
  ),
  Risk_Reduction_Percent = c(
    NA, NA, NA,
    (1 - mean(Annual_Risk_Cattle_low)/mean(Annual_Risk_Cattle)) * 100,
    (1 - mean(Annual_Risk_Buffalo_low)/mean(Annual_Risk_Buffalo)) * 100,
    (1 - mean(Annual_Risk_Combined_low)/mean(Annual_Risk_Combined)) * 100
  )
)

# Format the table for better readability
formatted_annual_risk <- annual_risk_summary
formatted_annual_risk$Mean_Annual_Risk <- sprintf("%.4f", formatted_annual_risk$Mean_Annual_Risk)
formatted_annual_risk$Median_Annual_Risk <- sprintf("%.4f", formatted_annual_risk$Median_Annual_Risk)
formatted_annual_risk$CI_Lower <- sprintf("%.4f", formatted_annual_risk$CI_Lower)
formatted_annual_risk$CI_Upper <- sprintf("%.4f", formatted_annual_risk$CI_Upper)
formatted_annual_risk$Risk_Reduction_Percent <- ifelse(is.na(formatted_annual_risk$Risk_Reduction_Percent), 
                                                       NA, 
                                                       sprintf("%.1f%%", formatted_annual_risk$Risk_Reduction_Percent))

# Create a more readable version of the table with CI combined
annual_risk_table <- data.frame(
  Population = formatted_annual_risk$Population,
  Scenario = formatted_annual_risk$Scenario,
  Median_Annual_Risk = formatted_annual_risk$Median_Annual_Risk,
  CI_95 = paste0("(", formatted_annual_risk$CI_Lower, "-", formatted_annual_risk$CI_Upper, ")"),
  Risk_Reduction = formatted_annual_risk$Risk_Reduction_Percent
)

# Display the table
print("Annual Infection Risk from Monthly Consumption of 100ml Milk (12 times per year)")
print(annual_risk_table)

# Save the table to CSV
write.csv(annual_risk_table, file = file.path(output_dir, "annual_risk_table.csv"), row.names = FALSE)

# Create a visualization of annual risks
png(file.path(output_dir, "Annual_Risk_Comparison.png"), width = 7, height = 5, units = 'in', res = 300)
par(mfrow = c(1, 1), mar = c(3, 4, 1, 0.3))

# Define colors and positions
shedding_colors <- c("pink3", "green4")
annual_risk_data <- matrix(c(
  mean(Annual_Risk_Cattle), mean(Annual_Risk_Cattle_low),
  mean(Annual_Risk_Buffalo), mean(Annual_Risk_Buffalo_low),
  mean(Annual_Risk_Combined), mean(Annual_Risk_Combined_low)
), ncol = 2, byrow = TRUE)


# Generate summary report
cat("=== Brucella abortus Risk Assessment Model Results ===\n\n")
cat("True Seroprevalence:\n")
cat("  Cattle: ", mean(TPrevC), " (95% CI: ", quantile(TPrevC, 0.025), "-", quantile(TPrevC, 0.975), ")\n")
cat("  Buffalo: ", mean(TPrevB), " (95% CI: ", quantile(TPrevB, 0.025), "-", quantile(TPrevB, 0.975), ")\n\n")

cat("Shedding Prevalence:\n")
cat("  Cattle: ", mean(TPrevShedC), " (95% CI: ", quantile(TPrevShedC, 0.025), "-", quantile(TPrevShedC, 0.975), ")\n")
cat("  Buffalo: ", mean(TPrevShedB), " (95% CI: ", quantile(TPrevShedB, 0.025), "-", quantile(TPrevShedB, 0.975), ")\n")
cat("  Combined: ", mean(PropOverallShed), " (95% CI: ", quantile(PropOverallShed, 0.025), "-", quantile(PropOverallShed, 0.975), ")\n\n")

cat("CFU in 100mL milk:\n")
cat("  Cattle: ", mean(CFU_countCattle100), " (95% CI: ", quantile(CFU_countCattle100, 0.025), "-", quantile(CFU_countCattle100, 0.975), ")\n")
cat("  Buffalo: ", mean(CFU_countBuff100), " (95% CI: ", quantile(CFU_countBuff100, 0.025), "-", quantile(CFU_countBuff100, 0.975), ")\n")
cat("  Combined: ", mean(CFU_countComb100), " (95% CI: ", quantile(CFU_countComb100, 0.025), "-", quantile(CFU_countComb100, 0.975), ")\n\n")

cat("Probability of infection per 100mL consumption:\n")
cat("  Cattle: ", mean(Prob_Infection_Cattle), " (95% CI: ", quantile(Prob_Infection_Cattle, 0.025), "-", quantile(Prob_Infection_Cattle, 0.975), ")\n")
cat("  Buffalo: ", mean(Prob_Infection_Buffalo), " (95% CI: ", quantile(Prob_Infection_Buffalo, 0.025), "-", quantile(Prob_Infection_Buffalo, 0.975), ")\n")
cat("  Combined: ", mean(Prob_Infection_Combined), " (95% CI: ", quantile(Prob_Infection_Combined, 0.025), "-", quantile(Prob_Infection_Combined, 0.975), ")\n\n")

cat("Probability of infection from a typical meal (", meal_size_ml, "mL) - With High Shedders:\n")
cat("  Cattle: ", mean(Meal_Risk_Cattle), " (95% CI: ", quantile(Meal_Risk_Cattle, 0.025), "-", quantile(Meal_Risk_Cattle, 0.975), ")\n")
cat("  Buffalo: ", mean(Meal_Risk_Buffalo), " (95% CI: ", quantile(Meal_Risk_Buffalo, 0.025), "-", quantile(Meal_Risk_Buffalo, 0.975), ")\n")
cat("  Combined: ", mean(Meal_Risk_Combined), " (95% CI: ", quantile(Meal_Risk_Combined, 0.025), "-", quantile(Meal_Risk_Combined, 0.975), ")\n\n")

cat("Probability of infection from a typical meal (", meal_size_ml, "mL) - Low Shedders Only:\n")
cat("  Cattle: ", mean(Meal_Risk_Cattle_low), " (95% CI: ", quantile(Meal_Risk_Cattle_low, 0.025), "-", quantile(Meal_Risk_Cattle_low, 0.975), ")\n")
cat("  Buffalo: ", mean(Meal_Risk_Buffalo_low), " (95% CI: ", quantile(Meal_Risk_Buffalo_low, 0.025), "-", quantile(Meal_Risk_Buffalo_low, 0.975), ")\n")
cat("  Combined: ", mean(Meal_Risk_Combined_low), " (95% CI: ", quantile(Meal_Risk_Combined_low, 0.025), "-", quantile(Meal_Risk_Combined_low, 0.975), ")\n\n")

cat("Risk reduction by removing high shedders (%):\n")
print(risk_reduction)

cat("\nMost influential parameters (from sensitivity analysis):\n")
print(sobol_data[order(-sobol_data$TotalEffect), ][1:3, ])

cat("\nAnnual Infection Risk (12 x 100ml consumption):\n")
cat("  With High Shedders:\n")
cat("    Cattle: ", mean(Annual_Risk_Cattle), " (95% CI: ", quantile(Annual_Risk_Cattle, 0.025), "-", quantile(Annual_Risk_Cattle, 0.975), ")\n")
cat("    Buffalo: ", mean(Annual_Risk_Buffalo), " (95% CI: ", quantile(Annual_Risk_Buffalo, 0.025), "-", quantile(Annual_Risk_Buffalo, 0.975), ")\n")
cat("    Combined: ", mean(Annual_Risk_Combined), " (95% CI: ", quantile(Annual_Risk_Combined, 0.025), "-", quantile(Annual_Risk_Combined, 0.975), ")\n\n")

cat("  Low Shedders Only:\n")
cat("    Cattle: ", mean(Annual_Risk_Cattle_low), " (95% CI: ", quantile(Annual_Risk_Cattle_low, 0.025), "-", quantile(Annual_Risk_Cattle_low, 0.975), ")\n")
cat("    Buffalo: ", mean(Annual_Risk_Buffalo_low), " (95% CI: ", quantile(Annual_Risk_Buffalo_low, 0.025), "-", quantile(Annual_Risk_Buffalo_low, 0.975), ")\n")
cat("    Combined: ", mean(Annual_Risk_Combined_low), " (95% CI: ", quantile(Annual_Risk_Combined_low, 0.025), "-", quantile(Annual_Risk_Combined_low, 0.975), ")\n\n")

cat("  Risk Reduction:\n")
cat("    Cattle: ", sprintf("%.1f%%", (1 - mean(Annual_Risk_Cattle_low)/mean(Annual_Risk_Cattle)) * 100), "\n")
cat("    Buffalo: ", sprintf("%.1f%%", (1 - mean(Annual_Risk_Buffalo_low)/mean(Annual_Risk_Buffalo)) * 100), "\n")
cat("    Combined: ", sprintf("%.1f%%", (1 - mean(Annual_Risk_Combined_low)/mean(Annual_Risk_Combined)) * 100), "\n")


