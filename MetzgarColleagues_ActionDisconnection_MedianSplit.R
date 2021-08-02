## Disconnection median split
library("dplyr")
library("ggplot2")
library("lme4")
library("readxl")
library("MASS")
## MEDIAN SPLIT 

# load data
data <- read_excel("disconnection_data.xlsx")

Factor <- c("GTS_HP", "MGI_HP")
data[Factor] <- lapply(data[Factor], as.numeric)

# find numbers for median split
n <- length(data$SubID)
n <- as.integer(n / 2)

# sort by lesion size
data <- data[order(data$LesionSize),]

# split data
data <- data[1:n,]

## GTS_HP ~ MGI_HP
model <- lm(GTS_HP ~ MGI_HP, na.action = na.exclude, data = data)
# add residuals to data structure
data$msGM_HP_resid <- residuals(model)
# studentized residuals
data$sGM_50 <- studres(model)
# put a 1 for all significant participants
data$GM[data$sGM_50 <= -1.701] <- 1
data$fGM_50 <- fitted(model)

## MGI_HP ~ GTS_HP
model <- lm(MGI_HP ~ GTS_HP, na.action = na.exclude, data = data)
data$msMG_HP_resid <- residuals(model)
data$sMG_50 <- studres(model)
data$MG[data$sMG_50 <= -1.701] <- 1
data$fMG_50 <- fitted(model)

## correlations (pearson is default)
# original behavior
cor.test(data$GTS_HP, data$MGI_HP)

# save
write.csv(data, "disconnection_residuals_mediansplit.csv")
