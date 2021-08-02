## Disconnection Residuals controlling for lesion volume

library("dplyr")
library("ggplot2")
library("lme4")
library("readxl")
library("MASS")

# load data
data <- read_excel("disconnection_data.xlsx")

Factor <- c("GTS_HP", "MGI_HP")
data[Factor] <- lapply(data[Factor], as.numeric)

### Lesion Volume Correction
## control for lesion volume
model <- lm(MGI_HP ~ LesionSize, na.action = na.exclude, data = data)
# add residuals to data structure
data$MGI_LVC <- residuals(model)
# studentize residuals
model <- lm(GTS_HP ~ LesionSize, na.action = na.exclude, data = data)
# put a 1 for all significant participants
data$GTS_LVC <- residuals(model)

# GTS ~ MGI
model <- lm(GTS_LVC ~ MGI_LVC, na.action = na.exclude, data = data)
data$GM_LVC_resid <- residuals(model)
data$sGM_LVC <- studres(model)
data$GM_LVC[data$sGM_LVC <= -1.672] <- 1

# MGI ~ GTS
model <- lm(MGI_LVC ~ GTS_LVC , na.action = na.exclude, data = data)
data$MG_LVC_resid <- residuals(model)
data$sMG_LVC <- studres(model)
data$MG_LVC[data$sMG_LVC <= -1.672] <- 1

## correlations (pearson is default)

# original behavior
cor.test(data$GTS_HP, data$MGI_HP)
# behavior and TLV
cor.test(data$GTS_HP, data$LesionSize)
cor.test(data$MGI_HP, data$LesionSize)
# LVC behavior
cor.test(data$MGI_LVC, data$GTS_LVC)

### SAVE
write.csv(data, "disconnection_residuals_LVC.csv")
