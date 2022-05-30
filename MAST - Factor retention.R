#### Packages and functions ####
lapply(c("tidyverse", "stargazer", "EGAnet", "lavaan", 
         "parallel", "foreach", "doParallel", "effects"), 
       library, character.only = TRUE)
numCores <- detectCores()

source("Sampling.R")
source("PCM.R")
source("Factor retention procedures.R")

#### Population correlation and loading matrices ####
#number of variables, major, minor and unique factors
J <- 20; M1 <- 3; M2 <- 200; M3 <- J

#in this study, the "high" scenario has communalities in the range 0.6-0.8
##the "low" scenario has communalities in the range 0.2-0.4 (Tucker et al., 1969)
set.seed(1)
communalities <- list(low = sample(c(0.2, 0.3, 0.4), size = J, replace = TRUE), 
                      high = sample(c(0.6, 0.7, 0.8), size = J, replace = TRUE)) 

#in addition interfactor correlation is set at two values
##the "high" scenario has an interfactor correlation of 0.5 (de Winter & Dodou, 2016)
##the "low" scenario has independent factors (correlation = 0)
ifc <- list(low = 0, high = 0.5)

#four population correlation matrices and loading matrices
##H(igh), I(fc), L(ow), C(ommunalities)
HILC <- PCM(ifc$high, communalities$low)
LILC <- PCM(ifc$low, communalities$low)
LIHC <- PCM(ifc$low, communalities$high)
HIHC <- PCM(ifc$high, communalities$high)

#### Factor retention simulation ####
#save the population correlation matrices in a list
PCM <- list(HIHC = HIHC$R, LIHC = LIHC$R, HILC = HILC$R, LILC = LILC$R)

#sample size of simulated sample
##we will consider two sample sizes: 100 and 1000 
N <- list(low = 100, high = 1000)

#empty file to store results
##for every correlation matrix the amount of factors according to the 
##Eigen Value rule, Scree Optimal Coordinate, Scree Acceleration Factor (Raîche et al., 2013)
write.table(t(c("corMatrix", "n", "nSamp", 
                paste0("Cont", c("EV", "ScreeAF", "ScreeAFEV", "PAM", "PA95", "MAP", "RPA", "EGA")), 
                paste0("Dich50Pearson", c("EV", "ScreeAF", "ScreeAFEV", "PAM", "PA95", "MAP", "RPA", "EGA")), 
                paste0("Dich75Pearson", c("EV", "ScreeAF", "ScreeAFEV", "PAM", "PA95", "MAP", "RPA", "EGA")),
                paste0("Dich50Tetra", c("EV","ScreeAF", "ScreeAFEV", "PAM", "PA95", "MAP", "RPA", "EGA")),
                paste0("Dich75Tetra", c("EV","ScreeAF", "ScreeAFEV", "PAM", "PA95", "MAP", "RPA", "EGA")))), 
            file = "ResultsFRC.csv", row.names = FALSE, col.names = FALSE)

registerDoParallel(numCores - 1)
foreach (corNumber = 1:length(PCM)) %dopar% {
  corMatrix <- PCM[[corNumber]]
  nameCorMatrix <- names(PCM[corNumber])

  for (n in N){

    for (nSamp in 1:1000){
      data <- sampling(N = n, popCorMatrix = corMatrix)
      
      #continuous data 
      corCont <- cor(data$Cont)
      
      #dichotomous data; 50-50 split; Pearson correlation
      corDich50Pearson <- cor(data$Dich50)
      
      #dichotomous data; 75-25 split; Pearson correlation
      corDich75Pearson <- cor(data$Dich75)
      
      #dichotomous data; 50-50 split; tetrachoric correlation
      corDich50Tetra <- lavaan::lavCor(data$Dich50, ordered = names(data$Dich50))
      
      #dichotomous data; 75-25 split; tetrachoric correlation
      corDich75Tetra <- lavaan::lavCor(data$Dich75, ordered = names(data$Dich75))
      
      tmp <- list(nameCorMatrix, n, nSamp, 
               EV(corCont), 
               ScreeAF(corCont)$ScreeAF, ScreeAF(corCont)$ScreeAFEV,
               PA(corCont, ncases = n)$PAM, PA(corCont, ncases = n)$PA95,
               MAP(corCont), revisedPA(corCont, ncases = n),
               EGAerror(corCont, n = n),
               
               EV(corDich50Pearson), 
               ScreeAF(corDich50Pearson)$ScreeAF, ScreeAF(corDich50Pearson)$ScreeAFEV, 
               PA(corDich50Pearson, ncases = n)$PAM, PA(corDich50Pearson, ncases = n)$PA95,
               MAP(corDich50Pearson), revisedPA(corDich50Pearson, ncases = n),
               EGAerror(corDich50Pearson, n = n),
               
               EV(corDich75Pearson),
               ScreeAF(corDich75Pearson)$ScreeAF, ScreeAF(corDich75Pearson)$ScreeAFEV,
               PA(corDich75Pearson, ncases = n)$PAM, PA(corDich75Pearson, ncases = n)$PA95,
               MAP(corDich75Pearson), revisedPA(corDich75Pearson, ncases = n),
               EGAerror(corDich75Pearson, n = n),
               
               EV(corDich50Tetra), 
               ScreeAF(corDich50Tetra)$ScreeAF, ScreeAF(corDich50Tetra)$ScreeAFEV,
               PA(corDich50Tetra, ncases = n)$PAM, PA(corDich50Tetra, ncases = n)$PA95,
               MAP(corDich50Tetra), revisedPA(corDich50Tetra, ncases = n),
               EGAerror(corDich50Tetra, n = n),
               
               EV(corDich75Tetra), 
               ScreeAF(corDich75Tetra)$ScreeAF, ScreeAF(corDich75Tetra)$ScreeAFEV,
               PA(corDich75Tetra, ncases = n)$PAM, PA(corDich75Tetra, ncases = n)$PA95,
               MAP(corDich75Tetra), revisedPA(corDich75Tetra, ncases = n),
               EGAerror(corDich75Tetra, n = n))
               
      write.table(tmp, file = "ResultsFRC.csv", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}


#### Data analysis ####
results <- read.table("ResultsFRC.csv", header = TRUE, sep = " ")
colMeans(sapply(results, is.na))
results <- results[-which(is.na(results$Dich75TetraMAP)), ]
#one row with missing values for Dich75TetraMAP, delete it
colMeans(sapply(results, is.na))
results$n <- as.factor(results$n)
table(results$corMatrix)
table(results$n)

resultsWide <- results %>%
  #correct mistake in scree plot criterion combined with Kaiser criterion
  mutate(ContScreeAFEV = ifelse(ContScreeAF > ContEV, ContScreeAFEV, ContScreeAF - 1), 
         Dich50PearsonScreeAFEV = ifelse(Dich50PearsonScreeAF > Dich50PearsonEV, Dich50PearsonScreeAFEV, Dich50PearsonScreeAF - 1), 
         Dich50TetraScreeAFEV = ifelse(Dich50TetraScreeAF > Dich50TetraEV, Dich50TetraScreeAFEV, Dich50PearsonScreeAF - 1),
         Dich75PearsonScreeAFEV = ifelse(Dich75PearsonScreeAF > Dich75PearsonEV, Dich75PearsonScreeAFEV, Dich75PearsonScreeAF - 1), 
         Dich75TetraScreeAFEV = ifelse(Dich75TetraScreeAF > Dich75TetraEV, Dich75TetraScreeAFEV, Dich75PearsonScreeAF - 1)) %>%
  
  #correct mistake in scree plot criterion
  mutate(ContScreeAF = ContScreeAF - 1, 
         Dich50PearsonScreeAF = Dich50PearsonScreeAF - 1, 
         Dich50TetraScreeAF = Dich50PearsonScreeAF - 1,
         Dich75PearsonScreeAF = Dich75PearsonScreeAF - 1, 
         Dich75TetraScreeAF = Dich75PearsonScreeAF - 1) %>%
  
  #also include RPAEV 
  mutate(ContRPAEV = pmin(ContEV, ContRPA), 
         Dich50PearsonRPAEV = pmin(Dich50PearsonEV, Dich50PearsonRPA),
         Dich50TetraRPAEV = pmin(Dich50TetraEV, Dich50TetraRPA),
         Dich75PearsonRPAEV = pmin(Dich75PearsonEV, Dich75PearsonRPA),
         Dich75TetraRPAEV = pmin(Dich75TetraEV, Dich75TetraRPA))

table(c(as.matrix(resultsWide[, 4:48])))

#all estimates whereby 20 factors are estimated should be converted to 0 
##this is also done for the EGAnet function
resultsWide[resultsWide >= 20] <- 0
table(c(as.matrix(resultsWide[, 4:48])))

#overview of the averages per method
colMeans(resultsWide[, 4:48], na.rm = TRUE)

resultsWide$nSamp <- NULL

#recode to long format for analyses
resultsLong <- gather(resultsWide, key = "method", value = "nFactors", ContEV:Dich75TetraRPAEV) 

#create factors for population correlation/type of data/correlation matrix/factor retention criterion
resultsLong <- resultsLong %>%
  mutate(interfactorCor = as.factor(ifelse(substr(corMatrix, 1, 1) == "H", "High", "Low")),
         communalities = as.factor(ifelse(substr(corMatrix, 3, 3) == "H", "High", "Low")),
         data = as.factor(substr(method, 1, 4)), 
         dataType = as.factor(ifelse(data == "Cont", "Cont", substr(method, 1, 6))),
         correlations = as.factor(ifelse(dataType == "Cont", "", 
                                      ifelse(substr(method, 7, 7) == "P", 
                                      "Pearson", "Tetra"))), 
         criterion = as.factor(substr(method, nchar(paste0(dataType, correlations)) + 1, nchar(method))), 
         correlations = as.factor(ifelse(dataType == "Cont", "Pearson", as.character(correlations))), 
         method = NULL, data = NULL, corMatrix = NULL)
         

table(resultsLong$n)
table(resultsLong$interfactorCor, resultsLong$communalities)
table(resultsLong$correlations, resultsLong$dataType)
table(resultsLong$dataType, resultsLong$criterion)

#bias variable and correctly predicted (Auerswald & Moshagen, 2019)
resultsLong$bias <- resultsLong$nFactors - 3
resultsLong$correct <- ifelse(resultsLong$bias == 0, 1, 0)

#for other chapters, linear regression makes more sense, not in this case
##lots of criteria -> plots with percentage correctly predicted more efficient
##lots of interactions -> effect sizes from ANOVA give more of an overview

#set contrasts to more easily interpretable ones
sapply(resultsLong[, c(1, 3:7)], contrasts)
contrasts(resultsLong$interfactorCor) <- contr.treatment(n = 2, base = 2)
contrasts(resultsLong$communalities) <- contr.treatment(n = 2, base = 2)
contrasts(resultsLong$criterion) <- contr.treatment(n = 9, base = 2)

sapply(resultsLong[, c(1, 3:7)], contrasts)

### LINEAR MODEL
linearFRC <- lm(bias ~ (n + interfactorCor + communalities + dataType + 
                          correlations)*criterion, data = resultsLong)

##create new data for plots
###study effect of n, interfactorCor and communalities for different levels of dataType and correlations
newData <- data.frame(n = c(rep("100", 216), rep("1000", 216)), 
                      interfactorCor = c(rep(c(rep("Low", 108), rep("High", 108)), 2)), 
                      communalities = c(rep(c(rep("Low", 54), rep("High", 54)), 4)), 
                      dataType = c(rep(c(rep("Cont", 18), rep("Dich50", 18), rep("Dich75", 18)), 4)), 
                      correlations = c(rep(c(rep("Pearson", 9), rep("Tetra", 9)), 12)), 
                      criterion = rep(c("EGA", "EV", "MAP", "PA95", "PAM",
                                        "RPA", "RPAEV", "ScreeAF", "ScreeAFEV"), 48)) 

newData$predictions <- predict(linearFRC, newdata = newData, type = "response")

newData$communalities <- factor(newData$communalities, levels=rev(levels(factor(newData$communalities))))
communalitiesLabs <- c("Low communalities", "High communalities")
names(communalitiesLabs) <- c("Low", "High")

newData$interfactorCor <- factor(newData$interfactorCor, levels=rev(levels(factor(newData$interfactorCor))))
interfactorCorLabs <- c("Low interfactor ρ", "High interfactor ρ")
names(interfactorCorLabs) <- c("Low", "High")

#continuous data, Pearson correlations
plot1 <- ggplot(data = filter(newData, dataType == "Cont", correlations == "Pearson"), 
                aes(y = predictions, x = n)) + 
  geom_line(aes(linetype = criterion, group = criterion, colour = criterion)) +
  geom_point(aes(shape = criterion, colour = criterion)) +
  facet_grid(interfactorCor ~ communalities, margins = FALSE, scales = "fixed",
             labeller = labeller(communalities = communalitiesLabs, 
                                 interfactorCor = interfactorCorLabs)) +
  theme_minimal() +
  labs(x = "Sample size", y = "Bias", 
       title = "Continuous data, Pearson correlations", 
       shape = "Criterion", linetype = "Criterion", colour = "Criterion") +
  scale_x_discrete(limits = c("100", "1000"), expand = c(0.05, 0.05)) +
  scale_shape_manual(values = c(1, 2, 3, 16, 17, 8, 0, 5, 6)) +
  scale_colour_manual(values = rep(c("grey13", "grey61", "dodgerblue4"), 3)) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave(plot1, file = "ContPearson.png", height = 4, width = 6)
       
#dichtomous data, 50 - 50 split, Pearson correlation
plot2 <- ggplot(data = filter(newData, dataType == "Dich50", correlations == "Pearson"), 
                aes(y = predictions, x = n)) + 
  geom_line(aes(linetype = criterion, group = criterion, colour = criterion)) +
  geom_point(aes(shape = criterion, colour = criterion)) +
  facet_grid(interfactorCor ~ communalities, margins = FALSE, scales = "fixed",
             labeller = labeller(communalities = communalitiesLabs, 
                                 interfactorCor = interfactorCorLabs)) +
  theme_minimal() +
  labs(x = "Sample size", y = "Bias", 
       title = "Dichotomous data (50-50), Pearson correlations", 
       shape = "Criterion", linetype = "Criterion", colour = "Criterion") +
  scale_x_discrete(limits = c("100", "1000"), expand = c(0.05, 0.05)) +
  scale_shape_manual(values = c(1, 2, 3, 16, 17, 8, 0, 5, 6)) +
  scale_colour_manual(values = rep(c("grey13", "grey61", "dodgerblue4"), 3)) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave(plot2, file = "Dich50Pearson.png", height = 4, width = 6)

#dichtomous data, 50 - 50 split, tetrachoric correlation
plot3 <- ggplot(data = filter(newData, dataType == "Dich50", correlations == "Tetra"), 
                aes(y = predictions, x = n)) + 
  geom_line(aes(linetype = criterion, group = criterion, colour = criterion)) +
  geom_point(aes(shape = criterion, colour = criterion)) +
  facet_grid(interfactorCor ~ communalities, margins = FALSE, scales = "fixed",
             labeller = labeller(communalities = communalitiesLabs, 
                                 interfactorCor = interfactorCorLabs)) +
  theme_minimal() +
  labs(x = "Sample size", y = "Bias", 
       title = "Dichotomous data (50-50), tetrachoric correlations", 
       shape = "Criterion", linetype = "Criterion", colour = "Criterion") +
  scale_x_discrete(limits = c("100", "1000"), expand = c(0.05, 0.05)) +
  scale_shape_manual(values = c(1, 2, 3, 16, 17, 8, 0, 5, 6)) +
  scale_colour_manual(values = rep(c("grey13", "grey61", "dodgerblue4"), 3)) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave(plot3, file = "Dich50Tetra.png", height = 4, width = 6)

#dichtomous data, 75 - 25 split, Pearson correlation
plot4 <- ggplot(data = filter(newData, dataType == "Dich75", correlations == "Pearson"), 
                aes(y = predictions, x = n)) + 
  geom_line(aes(linetype = criterion, group = criterion, colour = criterion)) +
  geom_point(aes(shape = criterion, colour = criterion)) +
  facet_grid(interfactorCor ~ communalities, margins = FALSE, scales = "fixed",
             labeller = labeller(communalities = communalitiesLabs, 
                                 interfactorCor = interfactorCorLabs)) +
  theme_minimal() +
  labs(x = "Sample size", y = "Bias", 
       title = "Dichotomous data (75-25), Pearson correlations", 
       shape = "Criterion", linetype = "Criterion", colour = "Criterion") +
  scale_x_discrete(limits = c("100", "1000"), expand = c(0.05, 0.05)) +
  scale_shape_manual(values = c(1, 2, 3, 16, 17, 8, 0, 5, 6)) +
  scale_colour_manual(values = rep(c("grey13", "grey61", "dodgerblue4"), 3)) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave(plot4, file = "Dich75Pearson.png", height = 4, width = 6)

#dichtomous data, 75 - 25 split, tetrachoric correlation
plot5 <- ggplot(data = filter(newData, dataType == "Dich75", correlations == "Tetra"), 
                aes(y = predictions, x = n)) + 
  geom_line(aes(linetype = criterion, group = criterion, colour = criterion)) +
  geom_point(aes(shape = criterion, colour = criterion)) +
  facet_grid(interfactorCor ~ communalities, margins = FALSE, scales = "fixed",
             labeller = labeller(communalities = communalitiesLabs, 
                                 interfactorCor = interfactorCorLabs)) +
  theme_minimal() +
  labs(x = "Sample size", y = "Bias", 
       title = "Dichotomous data (75-25), tetrachoric correlations", 
       shape = "Criterion", linetype = "Criterion", colour = "Criterion") +
  scale_x_discrete(limits = c("100", "1000"), expand = c(0.05, 0.05)) +
  scale_shape_manual(values = c(1, 2, 3, 16, 17, 8, 0, 5, 6)) +
  scale_colour_manual(values = rep(c("grey13", "grey61", "dodgerblue4"), 3)) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave(plot5, file = "Dich75Tetra.png", height = 4, width = 6)

### LOGISTIC MODEL

logisticFRC <- glm(correct ~ (n + interfactorCor + communalities + dataType + 
                                correlations)*criterion, data = resultsLong, 
                   family = binomial)
anova(logisticFRC, test = "Chisq")

newData$predictionsLog <- predict(logisticFRC, newdata = newData, type = "response")

#continuous data, Pearson correlations
plot6 <- ggplot(data = filter(newData, dataType == "Cont", correlations == "Pearson"), 
                aes(y = predictionsLog, x = n)) + 
  geom_line(aes(linetype = criterion, group = criterion, colour = criterion)) +
  geom_point(aes(shape = criterion, colour = criterion)) +
  facet_grid(interfactorCor ~ communalities, margins = FALSE, scales = "fixed",
             labeller = labeller(communalities = communalitiesLabs, 
                                 interfactorCor = interfactorCorLabs)) +
  theme_minimal() +
  labs(x = "Sample size", y = "P(Correct #factors)", 
       title = "Continuous data, Pearson correlations", 
       shape = "Criterion", linetype = "Criterion", colour = "Criterion") +
  scale_x_discrete(limits = c("100", "1000"), expand = c(0.05, 0.05)) +
  scale_shape_manual(values = c(1, 2, 3, 16, 17, 8, 0, 5, 6)) +
  scale_colour_manual(values = rep(c("grey13", "grey61", "dodgerblue4"), 3)) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave(plot6, file = "ContPearsonLog.png", height = 4, width = 6)

#dichtomous data, 50 - 50 split, Pearson correlation
plot7 <- ggplot(data = filter(newData, dataType == "Dich50", correlations == "Pearson"), 
                aes(y = predictionsLog, x = n)) + 
  geom_line(aes(linetype = criterion, group = criterion, colour = criterion)) +
  geom_point(aes(shape = criterion, colour = criterion)) +
  facet_grid(interfactorCor ~ communalities, margins = FALSE, scales = "fixed",
             labeller = labeller(communalities = communalitiesLabs, 
                                 interfactorCor = interfactorCorLabs)) +
  theme_minimal() +
  labs(x = "Sample size", y = "P(Correct #factors)", 
       title = "Dichotomous data (50-50), Pearson correlations", 
       shape = "Criterion", linetype = "Criterion", colour = "Criterion") +
  scale_x_discrete(limits = c("100", "1000"), expand = c(0.05, 0.05)) +
  scale_shape_manual(values = c(1, 2, 3, 16, 17, 8, 0, 5, 6)) +
  scale_colour_manual(values = rep(c("grey13", "grey61", "dodgerblue4"), 3)) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave(plot7, file = "Dich50PearsonLog.png", height = 4, width = 6)

#dichtomous data, 50 - 50 split, tetrachoric correlation
plot8 <- ggplot(data = filter(newData, dataType == "Dich50", correlations == "Tetra"), 
                aes(y = predictionsLog, x = n)) + 
  geom_line(aes(linetype = criterion, group = criterion, colour = criterion)) +
  geom_point(aes(shape = criterion, colour = criterion)) +
  facet_grid(interfactorCor ~ communalities, margins = FALSE, scales = "fixed",
             labeller = labeller(communalities = communalitiesLabs, 
                                 interfactorCor = interfactorCorLabs)) +
  theme_minimal() +
  labs(x = "Sample size", y = "P(Correct #factors)", 
       title = "Dichotomous data (50-50), tetrachoric correlations", 
       shape = "Criterion", linetype = "Criterion", colour = "Criterion") +
  scale_x_discrete(limits = c("100", "1000"), expand = c(0.05, 0.05)) +
  scale_shape_manual(values = c(1, 2, 3, 16, 17, 8, 0, 5, 6)) +
  scale_colour_manual(values = rep(c("grey13", "grey61", "dodgerblue4"), 3)) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave(plot8, file = "Dich50TetraLog.png", height = 4, width = 6)

#dichtomous data, 75 - 25 split, Pearson correlation
plot9 <- ggplot(data = filter(newData, dataType == "Dich75", correlations == "Pearson"), 
                aes(y = predictionsLog, x = n)) + 
  geom_line(aes(linetype = criterion, group = criterion, colour = criterion)) +
  geom_point(aes(shape = criterion, colour = criterion)) +
  facet_grid(interfactorCor ~ communalities, margins = FALSE, scales = "fixed",
             labeller = labeller(communalities = communalitiesLabs, 
                                 interfactorCor = interfactorCorLabs)) +
  theme_minimal() +
  labs(x = "Sample size", y = "P(Correct #factors)", 
       title = "Dichotomous data (75-25), Pearson correlations", 
       shape = "Criterion", linetype = "Criterion", colour = "Criterion") +
  scale_x_discrete(limits = c("100", "1000"), expand = c(0.05, 0.05)) +
  scale_shape_manual(values = c(1, 2, 3, 16, 17, 8, 0, 5, 6)) +
  scale_colour_manual(values = rep(c("grey13", "grey61", "dodgerblue4"), 3)) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave(plot9, file = "Dich75PearsonLog.png", height = 4, width = 6)

#dichtomous data, 75 - 25 split, tetrachoric correlation
plot10 <- ggplot(data = filter(newData, dataType == "Dich75", correlations == "Tetra"), 
                aes(y = predictionsLog, x = n)) + 
  geom_line(aes(linetype = criterion, group = criterion, colour = criterion)) +
  geom_point(aes(shape = criterion, colour = criterion)) +
  facet_grid(interfactorCor ~ communalities, margins = FALSE, scales = "fixed",
             labeller = labeller(communalities = communalitiesLabs, 
                                 interfactorCor = interfactorCorLabs)) +
  theme_minimal() +
  labs(x = "Sample size", y = "P(Correct #factors)", 
       title = "Dichotomous data (75-25), tetrachoric correlations", 
       shape = "Criterion", linetype = "Criterion", colour = "Criterion") +
  scale_x_discrete(limits = c("100", "1000"), expand = c(0.05, 0.05)) +
  scale_shape_manual(values = c(1, 2, 3, 16, 17, 8, 0, 5, 6)) +
  scale_colour_manual(values = rep(c("grey13", "grey61", "dodgerblue4"), 3)) +
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave(plot10, file = "Dich75TetraLog.png", height = 4, width = 6)



