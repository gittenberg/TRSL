setwd("~/git/TRSL/parameters")

# load data from csv
data = read.csv("regression_data.csv", sep=",")
# rename columns
colnames(data) = c("X", "gene", "translation.rate", "transcript.abundance", "initiation.rate", "ORF.length", "CAI")

# unnormalized model
model_unnorm = lm(data$translation.rate ~ data$initiation.rate + data$CAI)
summary(model_unnorm)

# scale (standard normalize) variables
IR_norm = scale(data$initiation.rate)
CAI_norm = scale(data$CAI)

# normalized model
model = lm(data$translation.rate ~ IR_norm + CAI_norm)
summary(model)

# matrix of values, only complete cases. we include the dependent variable for later plotting.
factors_norm = cbind(IR_norm, CAI_norm, data$translation.rate)[complete.cases(cbind(IR_norm, CAI_norm)),]
# add vector of ones for score calculation
factors_norm = cbind(rep(1, nrow(factors_norm)), factors_norm)
colnames(factors_norm) = c('const', 'IR_norm', 'CAI_norm', 'translation.rate')

# define score s = t + b1 * x1 + b2 * x2
scores = factors_norm[, c('const', 'IR_norm', 'CAI_norm')] %*% model$coefficients
plot.new()
curve(x^1, from=0, to=20, col="red", xlab="score", ylab="translation rate")
points(scores, factors_norm[, "translation.rate"], ylim=c(0,50), main="", pch=1, col='blue')
# this is obviously not a good model.
# next we plot against log(translation rate):

# normalized model
log_model = lm(log(data$translation.rate) ~ IR_norm + CAI_norm)
summary(log_model)

# matrix of values, only complete cases. we include the dependent variable for later plotting.
factors_norm = cbind(IR_norm, CAI_norm, log(data$translation.rate))[complete.cases(cbind(IR_norm, CAI_norm)),]
# add vector of ones for score calculation
factors_norm = cbind(rep(1, nrow(factors_norm)), factors_norm)
colnames(factors_norm) = c('const', 'IR_norm', 'CAI_norm', 'log.translation.rate')

# define score s = t + b1 * x1 + b2 * x2
scores = factors_norm[, c('const', 'IR_norm', 'CAI_norm')] %*% log_model$coefficients
plot.new()
curve(x^1, from=-5, to=8, col="red", xlab="score", ylab="log translation rate")
points(scores, factors_norm[, "log.translation.rate"], main="", xlab="score", ylab="log(translation rate)", pch=1, col='blue')
# much better fit!
