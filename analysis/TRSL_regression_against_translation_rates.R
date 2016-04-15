setwd("~/git/TRSL/parameters")

# load data from csv
data = read.csv("regression_data.csv", sep=",")
# rename columns
colnames(data) = c("X", "gene", "translation.rate", "transcript.abundance", "initiation.rate", "ORF.length", "CAI")

###################################################################################################################
# unnormalized model
###################################################################################################################
model_unnorm = lm(data$translation.rate ~ data$initiation.rate + data$CAI)
summary(model_unnorm)
# the coefficients are not comparable because the factors have different orders of magnitude

# scale (standard normalize) variables
IR_norm = scale(data$initiation.rate)
CAI_norm = scale(data$CAI)

###################################################################################################################
# normalized model
###################################################################################################################
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

###################################################################################################################
# logarithmic model
###################################################################################################################
log_model = lm(log(data$translation.rate) ~ IR_norm + CAI_norm)
summary(log_model)

# matrix of values, only complete cases. we include the dependent variable for later plotting.
factors_log = cbind(IR_norm, CAI_norm, log(data$translation.rate))[complete.cases(cbind(IR_norm, CAI_norm)),]
# add vector of ones for score calculation
factors_log = cbind(rep(1, nrow(factors_log)), factors_log)
colnames(factors_log) = c('const', 'IR_norm', 'CAI_norm', 'log.translation.rate')

# define score s = t + b1 * x1 + b2 * x2
scores = factors_log[, c('const', 'IR_norm', 'CAI_norm')] %*% log_model$coefficients
plot.new()
curve(x^1, from=-5, to=8, col="red", xlab="score", ylab="log translation rate")
points(scores, factors_log[, "log.translation.rate"], main="", pch=1, col='blue')
# much better fit!

###################################################################################################################
# transformed model
###################################################################################################################
transf_data = (log(data$translation.rate)-min(log(data$translation.rate), na.rm=TRUE))^2
transf_model = lm(transf_data ~ IR_norm + CAI_norm)
summary(transf_model)

# matrix of values, only complete cases. we include the dependent variable for later plotting.
factors_transf = cbind(IR_norm, CAI_norm, transf_data)[complete.cases(cbind(IR_norm, CAI_norm)),]
# add vector of ones for score calculation
factors_transf = cbind(rep(1, nrow(factors_transf)), factors_transf)
colnames(factors_transf) = c('const', 'IR_norm', 'CAI_norm', 'transformed.translation.rate')

# define score s = t + b1 * x1 + b2 * x2
scores = factors_transf[, c('const', 'IR_norm', 'CAI_norm')] %*% transf_model$coefficients
plot.new()
curve(x^1, from=0, to=150, col="red", xlab="score", ylab="transformed translation rate")
points(scores, factors_transf[, "transformed.translation.rate"], main="", pch=1, col='blue')
# looks good!
