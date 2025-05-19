library("astsa")
library(stats)
library(nnet)

data <- bnrf1ebv

x <- data[1:1000]
y <- data[2:1001]

# create a multinomial logit model
m <- multinom(as.factor(y) ~ as.factor(x), Hess = TRUE, trace = FALSE)

# print summary of model results
summary(m)

n <- length(y)
k <- length(coef(m))
log_likelihood <- logLik(m)
aic <- -2 * log_likelihood + 2 * k
bic <- -2 * log_likelihood + k * log(n)

# print AIC and BIC values
cat("AIC:", aic, "\n")
cat("BIC:", bic, "\n")


new_x <- data[991:1000]
predictions <- predict(m, newdata = data.frame(x = new_x), type = "probs")
next_nucleotides <- apply(predictions, 1, function(x) sample(c("G", "A", "C", "T"), size = 1, prob = x))


actual <- data[1001:1010]
actual
next_nucleotides
