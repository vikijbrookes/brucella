#####################################################################
####               Balbir and Vijay - RA Brucella abortus        ####
#####################################################################
library(ggplot2)
library(epiR)
library(MASS)
options(scipen = 999)

iter = 10000

#---------------------------------------------------------------------------
#### Probability that a cow or buffalo is sero-positive
# Data from Ludhiana district

#overall
s = 41
n = 261

App_Seroprev = rbeta(iter, s +1, n - s +1)
hist(App_Seroprev, breaks = 'fd', col = 'red3', main = "", xlim = c(0, 0.3), xlab = "Apparent seroprevalence - overall", ylim = c(0, 1000))
summary(App_Seroprev);quantile(App_Seroprev, c(0.025, 0.975))

#Cattle
sC = 33
nC = 181
App_SeroprevC = rbeta(iter, sC +1, nC - sC +1)
hist(App_SeroprevC, breaks = 'fd', col = 'red3', main = "", xlim = c(0, 0.3), xlab = "Apparent seroprevalence - cattle", ylim = c(0, 1000))
summary(App_SeroprevC);quantile(App_SeroprevC, c(0.025, 0.975))

#
sB = 8
nB = 80
App_SeroprevB = rbeta(iter, sB +1, nB - sB +1)
hist(App_SeroprevB, breaks = 'fd', col = 'red3', main = "", xlim = c(0, 0.3), xlab = "Apparent seroprevalence - buffalo", ylim = c(0, 1000))
summary(App_SeroprevB);quantile(App_SeroprevB, c(0.025, 0.975))


Se1 = 0.874 # RBPT 
Se2 = 0.846 # I-ELISA
Sp1 = 0.994 # RBPT 
Sp2 = 0.996 # I-ELISA

Se_series = Se1 * Se2
Sp_series = 1- (1 - Sp1) * (1 - Sp2)

Se_parallel = 1 - (1 - Se1) * (1 - Se2)
Sp_parallel = Sp1 * Sp2

TPrevO = c()
for (i in 1:iter){
  TruePrevO = (sample(App_Seroprev, 1) + Sp_series - 1)/ (Se_series + Sp_series -1)
  TPrevO = c(TPrevO, TruePrevO)
}
hist(TPrevO, breaks = 'fd', main = "", xlim = c(0, 0.4), xlab = "True prevalence", ylim = c(0, 800))
summary(TPrevO);quantile(TPrevO, c(0.025, 0.975))

TPrevC = c()
for (i in 1:iter){
  TruePrevC = (sample(App_SeroprevC, 1) + Sp_series - 1)/ (Se_series + Sp_series -1)
  TPrevC = c(TPrevC, TruePrevC)
}
hist(TPrevC, breaks = 'fd', main = "", xlim = c(0, 0.4), xlab = "True prevalence", ylim = c(0, 800))
summary(TPrevC);quantile(TPrevC, c(0.025, 0.975))

TPrevB = c()
for (i in 1:iter){
  TruePrevB = (sample(App_SeroprevB, 1) + Sp_series - 1)/ (Se_series + Sp_series -1)
  TPrevB = c(TPrevB, TruePrevB)
}
hist(TPrevB, breaks = 'fd', main = "", xlim = c(0, 0.4), xlab = "True prevalence", ylim = c(0, 800))
summary(TPrevB);quantile(TPrevB, c(0.025, 0.975))

setwd("C:/Users/vbro3295/OneDrive - The University of Sydney (Staff)/People_Places/Balbir/Brucella/Figures_Brucella/")
png("TruePrev.png", width = 5, height = 4, units = 'in', res= 300)
par(mfrow = c(1, 1), mar = c(3, 4, 0.3, 0.2))
boxplot(TPrevC, TPrevB, TPrevO, names = c("Cattle", "Buffalo", "Combined"), ylab = "True prevalence")
dev.off()

# Probability that a sero-positive animal sheds in milk --------------------
# Seronegative animals did not shed in milk
shedders = 6 #101
seropos = 260 #301

sheddersItaly = 101
seroposItaly = 337

# GS so AP is TP

# Most conservative (buffalo in Italy)
Shed_PrevItaly = rbeta(iter, sheddersItaly +1, seroposItaly - sheddersItaly +1)
hist(Shed_PrevItaly, breaks = 50, col = 'grey50', main = "", xlim = c(0, 0.5), ylim = c(0, 1000), xlab = "Shedding proportion (of seropositive animals)")
summary(Shed_PrevItaly);quantile(Shed_PrevItaly, c(0.025, 0.975))
shedItaly = as.data.frame(Shed_PrevItaly)

# Least conservative (cows in Iran)
Shed_PrevIran = rbeta(iter, shedders +1, seropos - shedders +1)
hist(add = T, Shed_PrevIran, breaks = 50, col = 'grey70', main = "", xlim = c(0, 0.5), xlab = "Shedding proportion", ylim = c(0, 800))
summary(Shed_PrevIran);quantile(Shed_PrevIran, c(0.025, 0.975))
shedIran = as.data.frame(Shed_PrevIran)


# Plots
setwd("C:/Users/vbro3295/OneDrive - The University of Sydney (Staff)/People_Places/Balbir/Brucella/Figures_Brucella/")
png("Shedding.png", width = 5, height = 4, units = 'in', res= 300)
par(mfrow = c(1, 1), mar = c(3, 4, 0.3, 0.2))
### Overall prevalence of shedding cows and buffalo in Punjab --------------

# High (Buffalo Italy) BUFFALO
  TPrevShedB = Shed_PrevItaly * TPrevB
  # hist(TPrevShedB, breaks = 50, col = '#7FC97F', main = "", xlim = c(0, 0.12), xlab = "Shedding prevalence", ylim = c(0, 1000))
  summary(TPrevShedB);quantile(TPrevShedB, c(0.025, 0.975))
# Low (Cows in Tehran, Iran) COWS
  TPrevShedC = Shed_PrevIran * TPrevC
  # plot(density(add = T, TPrevShedC, breaks = 50,  main = "", xlim = c(0, 0.15), ylim = c(0, 800), line = "#FFFF99"))
  summary(TPrevShedC);quantile(TPrevShedC, c(0.025, 0.975))

# Overall shedding prev given the proportions f cows and buffalo in the dataset (181 cows and 80 buffalo)
PropCows = 181/(181+80)
  CountCattle = rbinom(1000, 1000, PropCows)
  ShedCountCattle = round(CountCattle * TPrevShedC, 0)
  CountBuff = 1000 - CountCattle
  ShedCountBuff = round(CountBuff * TPrevShedB, 0)
  PropOverallShed = (ShedCountCattle + ShedCountBuff)/1000

  # hist(add = T, PropOverallShed, breaks = 'fd', main = "", xlim = c(0, 0.4), ylim = c(0, 800), col = "#386CB0")
summary(PropOverallShed);quantile(PropOverallShed, c(0.025, 0.975))

boxplot(TPrevShedC, TPrevShedB, PropOverallShed, names = c("Cattle", "Buffalo", "Combined"), ylab = "Shedding prevalence")

dev.off()

library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all()
brewer.pal(n=5,"Accent")



# Brucella CFU/ml milk -----------------------------------------------------
# 73 Low shedder 30-999
# 16 high shedder 10000-40000
# 12 intermittent low shedder - Make assumption that these are sometimes below limit of detection of 30CFU/ml

LShed = 73 + 12
HShed = 16
TotShed = LShed + HShed

PLowShed = rbeta(10000, LShed +1,  TotShed - LShed +1)
hist(PLowShed, breaks = 'fd', main = "", xlim = c(0, 1), xlab = "Shedding prevalence", ylim = c(0, 800))
summary(PLowShed);quantile(PLowShed, c(0.025, 0.975))

RLowShed = runif(iter, 0, 999) 

PHighShed = 1 - PLowShed
hist(PHighShed, breaks = 'fd', main = "", xlim = c(0, 1), xlab = "Shedding prevalence", ylim = c(0, 800))
summary(PHighShed);quantile(PHighShed, c(0.025, 0.975))

RHighShed = runif(iter, 10000, 40000) 

#---- CFU/ml milk from shedding animals 
# set up empty vector for CFU counts per ml dependent on whether infected milk is from a high or low shedder
CFU_counts = c()

# Now run a loop to collect CFU per ml infected milk
for (i in 1:iter){
  PLowShed1 = rbeta(1, LShed +1,  TotShed - LShed +1) # alpha1 = x+1, alpha2 = n-x+1
  random1 = runif(1, 0, 1)
  
  if (random1 < PLowShed1) {
    CFU_count = sample(RLowShed, 1, replace = T)
  } else {
    CFU_count = sample(RHighShed, 1, replace = T)
  }
  CFU_counts = c(CFU_counts, CFU_count)
}

hist(CFU_counts, col = 'yellow', main = '', breaks = 100, xlab = 'CFU/ml infected milk')

# Now run a loop to get CFU/ml milk from unregulated sources
# First based on high prevalence of shedders (Italy), then low prev of shedders (Iran)
# Do not assume that volume of milk produced by individual cows is from the same distribution regardless of infection status.
# Therefore, the proportion of 'infected' milk is 10% less than the prevalence of shedders.

# Cow scenario

TotalMilkSheddingCattle = TPrevShedC * 0.9
TotalMilkNonShedCattle = (1-TPrevShedC)
PrecentSheddingCattle = TotalMilkSheddingCattle/(TotalMilkNonShedCattle + TotalMilkSheddingCattle)
summary(PrecentSheddingCattle);quantile(PrecentSheddingCattle, c(0.025, 0.975))

TotalMilkSheddingBuff = TPrevShedB * 0.9
TotalMilkNonShedBuff = (1-TPrevShedB)
PrecentSheddingBuff = TotalMilkSheddingBuff/(TotalMilkNonShedBuff + TotalMilkSheddingBuff)
summary(PrecentSheddingBuff);quantile(PrecentSheddingBuff, c(0.025, 0.975))

TotalMilkSheddingComb = PropOverallShed * 0.9
TotalMilkNonShedComb = (1- PropOverallShed)
PrecentSheddingComb = TotalMilkSheddingComb/(TotalMilkNonShedComb + TotalMilkSheddingComb)
summary(PrecentSheddingComb);quantile(PrecentSheddingComb, c(0.025, 0.975))

# in 100ml, milk will come from the distribution of CFU/ml
CFU_countCattle100 = c()
for (i in 1:iter){
  MilkCFU = sum(sample(CFU_counts, round((PrecentSheddingCattle[i] *100), 0)))
  CFU_countCattle100 = c(CFU_countCattle100, MilkCFU)
}

CFU_countBuff100 = c()
for (i in 1:iter){
  MilkCFU = sum(sample(CFU_counts, round((PrecentSheddingBuff[i] *100), 0)))
  CFU_countBuff100 = c(CFU_countBuff100, MilkCFU)
}

CFU_countComb100 = c()
for (i in 1:iter){
  MilkCFU = sum(sample(CFU_counts, round((PrecentSheddingComb[i] *100), 0)))
  CFU_countComb100 = c(CFU_countComb100, MilkCFU)
}

summary(CFU_countCattle100);quantile(CFU_countCattle100, c(0.025, 0.975))
summary(CFU_countBuff100);quantile(CFU_countBuff100, c(0.025, 0.975))
summary(CFU_countComb100);quantile(CFU_countComb100, c(0.025, 0.975))

# Calculate percentages for each group
percent_cattle <- mean(CFU_countCattle100 > 1000) * 100
percent_buffalo <- mean(CFU_countBuff100 > 1000) * 100
percent_combined <- mean(CFU_countComb100 > 1000) * 100

setwd("C:/Users/vbro3295/OneDrive - The University of Sydney (Staff)/People_Places/Balbir/Brucella/Figures_Brucella/")
png("CFU100ml.png", width = 5, height = 4, units = 'in', res= 300)
par(mfrow = c(1, 1), mar = c(3, 4, 1, 0.3))

boxplot(CFU_countCattle100, CFU_countBuff100, CFU_countComb100, names = c("Cattle", "Buffalo", "Combined"), ylab = "CFU/100ml")
abline(h = 1000, lty = 2, col = 2, lwd = 2)
dev.off()

######################################################################################################
############                  Model for sensitivity analysis                                  ########
######################################################################################################
library(sensitivity)

n <- 10000

Viki.Sobol <- function (X) {

  # Prevalence of shedding for each species
  TPrevShC = X$PropSheddingC * X$Seroprev_cattle
  #print(TPrevShC)
  TPrevShB = X$PropSheddingB * X$Seroprev_buffalo 
  #print(TPrevShB)
  
  # Prevalence of shedding overall
  CountCattleS = round((X$PropCows * TPrevShC), 0)
  NoBuffalo = 1000 - X$PropCows
  CountBuffaloS = round((NoBuffalo * TPrevShB), 0)
  PropOverallShed = (CountCattleS + CountBuffaloS)/1000
  
  # CFU in milk
  CFUShedders = (X$ProbLowShedding * X$AmountLowShed) + (1-X$ProbLowShedding * X$AmountHighShed) 

  CFUall = CFUShedders * PropOverallShed
  return(CFUall)
}


X1 <- data.frame(Seroprev_cattle = sample(TPrevC, n, replace = T),
                 Seroprev_buffalo = sample(TPrevB, n, replace = T),
                 PropSheddingB = sample(Shed_PrevItaly, n, replace = T),
                 PropSheddingC = sample(Shed_PrevIran, n, replace = T),
                 PropCows = sample(CountCattle, n, replace = T),
                 ProbLowShedding = sample(PLowShed, n, replace = T),
                 AmountLowShed = sample(RLowShed, n, replace = T),
                 AmountHighShed = sample(RHighShed, n, replace = T))

X2 <- data.frame(Seroprev_cattle = sample(TPrevC, n, replace = T),
                 Seroprev_buffalo = sample(TPrevB, n, replace = T),
                 PropSheddingB = sample(Shed_PrevItaly, n, replace = T),
                 PropSheddingC = sample(Shed_PrevIran, n, replace = T),
                 PropCows = sample(CountCattle, n, replace = T),
                 ProbLowShedding = sample(PLowShed, n, replace = T),
                 AmountLowShed = sample(RLowShed, n, replace = T),
                 AmountHighShed = sample(RHighShed, n, replace = T))

x1 <- sobol2007(model = Viki.Sobol, X1 = X1, X2 = X2, nboot = 1000)
print(x1)

plot(x1)
title(xlab = "Parameter", ylab = "Index", main = "")
abline (h = 0.5, col = "red")
abline (h = 0.25, col = "red")

library(ggplot2)
library(sensitivity)

# Extracting main and total effects from the sobol2007 results
sobol_data <- data.frame(
  Parameter = rownames(x1$S),
  MainEffect = x1$S[, "original"],
  TotalEffect = x1$T[, "original"]
)
rownames(sobol_data) <- NULL  # Reset rownames for ggplot compatibility

# Custom labels for parameters
custom_labels <- c(
  "Seroprev_cattle" = "Seroprevalence cattle",
  "Seroprev_buffalo" = "Seroprevalence buffalo",
  "PropSheddingB" = "Shedding proportion buffalo",
  "PropSheddingC" = "Shedding Proportion cattle",
  "PropCows" = "Proportion of cattle",
  "ProbLowShedding" = "Shedding type (low/high)",
  "AmountLowShed" = "Low shedding amount",
  "AmountHighShed" = "High shedding amount"
)

# Reshaping data for ggplot
sobol_data_long <- reshape2::melt(sobol_data, id.vars = "Parameter", 
                                  variable.name = "EffectType", 
                                  value.name = "Value")

setwd("C:/Users/vbro3295/OneDrive - The University of Sydney (Staff)/People_Places/Balbir/Brucella/Figures_Brucella/")
png("Sobol.png", width = 6, height = 5, units = 'in', res= 300)
# Create the plot
ggplot(sobol_data_long, aes(x = Parameter, y = Value, color = EffectType)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +  # Increased dot size
  scale_x_discrete(labels = custom_labels) +
  scale_color_manual(values = c("MainEffect" = "#E41A1C", "TotalEffect" = "#377EB8")) +
  labs(x = "Parameter", y = "Sensitivity Index", color = "Effect Type") +
  theme_minimal() +
  ylim(0, 0.8) +
  geom_hline(yintercept = c(0.25, 0.5), linetype = "dashed", color = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 6, r = 6, b = 6, l = 20)  # Increased left margin
  )
dev.off()


######################################################################################################
### Redo the sensitivity analysis with all low shedding
######################################################################################################

# Initially the same code.
# Brucella CFU/ml milk -----------------------------------------------------
# Here, High shedding becomes 0 probability

#---- CFU/ml milk from shedding animals 
# set up empty vector for CFU counts per ml dependent on whether infected milk is from a high or low shedder
CFU_counts = RLowShed

hist(CFU_counts, col = 'yellow', main = '', breaks = 100, xlab = 'CFU/ml infected milk')

# Now run a loop to get CFU/ml milk from unregulated sources
# First based on high prevalence of shedders (Italy), then low prev of shedders (Iran)
# Do not assume that volume of milk produced by individual cows is from the same distribution regardless of infection status.
# Therefore, the proportion of 'infected' milk is 10% less than the prevalence of shedders.

# Cow scenario

TotalMilkSheddingCattle = TPrevShedC * 0.9
TotalMilkNonShedCattle = (1-TPrevShedC)
PrecentSheddingCattle = TotalMilkSheddingCattle/(TotalMilkNonShedCattle + TotalMilkSheddingCattle)
summary(PrecentSheddingCattle);quantile(PrecentSheddingCattle, c(0.025, 0.975))

TotalMilkSheddingBuff = TPrevShedB * 0.9
TotalMilkNonShedBuff = (1-TPrevShedB)
PrecentSheddingBuff = TotalMilkSheddingBuff/(TotalMilkNonShedBuff + TotalMilkSheddingBuff)
summary(PrecentSheddingBuff);quantile(PrecentSheddingBuff, c(0.025, 0.975))

TotalMilkSheddingComb = PropOverallShed * 0.9
TotalMilkNonShedComb = (1- PropOverallShed)
PrecentSheddingComb = TotalMilkSheddingComb/(TotalMilkNonShedComb + TotalMilkSheddingComb)
summary(PrecentSheddingComb);quantile(PrecentSheddingComb, c(0.025, 0.975))

# in 100ml, milk will come from the distribution of CFU/ml
CFU_countCattle100 = c()
for (i in 1:iter){
  MilkCFU = sum(sample(CFU_counts, round((PrecentSheddingCattle[i] *100), 0)))
  CFU_countCattle100 = c(CFU_countCattle100, MilkCFU)
}

CFU_countBuff100 = c()
for (i in 1:iter){
  MilkCFU = sum(sample(CFU_counts, round((PrecentSheddingBuff[i] *100), 0)))
  CFU_countBuff100 = c(CFU_countBuff100, MilkCFU)
}

CFU_countComb100 = c()
for (i in 1:iter){
  MilkCFU = sum(sample(CFU_counts, round((PrecentSheddingComb[i] *100), 0)))
  CFU_countComb100 = c(CFU_countComb100, MilkCFU)
}

summary(CFU_countCattle100);quantile(CFU_countCattle100, c(0.025, 0.975))
summary(CFU_countBuff100);quantile(CFU_countBuff100, c(0.025, 0.975))
summary(CFU_countComb100);quantile(CFU_countComb100, c(0.025, 0.975))

