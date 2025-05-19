library(readr)

# Download the FASTA file
link <- "https://www.ebi.ac.uk/ena/browser/api/fasta/CP009973.1?lineLimit=1000"
data <- readLines(link)

# Remove the header line from the data
data <- data[-1]

# Concatenate all lines into a single string
data <- paste(data, collapse = "")

# Separate the DNA sequence into individual nucleotides
data <- unlist(strsplit(data, ""))

x <- character(length(data) - 10)
for (i in 1:(length(data) - 10)) {
  x[i] <- paste0(data[i], data[i+1], data[i+2], data[i+3], data[i+4], data[i+5], data[i+6],data[i+7],data[i+8],data[i+9])
}

data <- x[1:100]
data
noquote(data)
