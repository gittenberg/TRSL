setwd("~/git/TRSL/parameters")

# load data from csv
data = read.csv("regression_data.csv", sep=",")
# rename columns
colnames(data) = c("X", "gene", "translation.rate", "transcript.abundance", "initiation.rate", "ORF.length", "CAI")

# unnormalized model
model_unnorm = lm(data$translation.rate ~ data$initiation.rate + data$CAI)
summary(model_unnorm)

IR_norm = scale(data$initiation.rate)
CAI_norm = scale(data$CAI)

# normalized model
model = lm(data$translation.rate ~ IR_norm + CAI_norm)
summary(model)
