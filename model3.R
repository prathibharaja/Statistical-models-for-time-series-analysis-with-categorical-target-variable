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

x <- character(length(data) - 3)
for (i in 1:(length(data) - 3)) {
  x[i] <- paste0(data[i], data[i+1], data[i+2])
}
x <- as.integer(x)
x
x <- x[1:1000]
data <- as.integer(data)
y <- data[4:1003]
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

