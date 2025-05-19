library(readr)
library(stats)
library(nnet)

# Download the FASTA file
link <- "https://www.ebi.ac.uk/ena/browser/api/fasta/CP009973.1?lineLimit=1000"
data <- readLines(link)

# Remove the header line from the data
data <- data[-1]

# Concatenate all lines into a single string
data <- paste(data, collapse = "")

# Replace nucleotide characters with corresponding numbers
data <- gsub("G", "1", data)
data <- gsub("A", "2", data)
data <- gsub("C", "3", data)
data <- gsub("T", "4", data)

# Separate the DNA sequence into individual nucleotides
data <- unlist(strsplit(data, ""))

# Cast the character vector to an integer vector
data <- as.integer(data)

x <- data[1:1000]
y <- data[2:1001]
m <- multinom(as.factor(y) ~ as.factor(x), Hess = TRUE, trace = FALSE)

# print summary of model results
summary(m, Wald.ratios = TRUE)

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
