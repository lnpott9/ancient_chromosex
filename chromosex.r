# chromosex.r
# Written by Laura Pott
# December 14 2023

# Input: idxstats file
# Output: estimate for the individual's chromosomal sex
# Usage: This script is designed to be run as part of a shell script, but can be run on an idxstats file individually.
######## Rscript chromosex.r <idxstats file>

# 1. Allow script to accept command-line arguments
if (length(commandArgs(trailingOnly = TRUE)) != 1) {
  stop("Usage: Rscript chromosex.r <input_file>")
}

# 2. Get the input file from the command-line argument
input_file <- commandArgs(trailingOnly = TRUE)[1]

# 3. Remove file extension, then extract the base file name
file_name_noextension <- tools::file_path_sans_ext(basename(input_file))
base_file_name <- sub("_rmdup", "", file_name_noextension)

# 4. Read idxstats file into a data frame, with no header, 24 rows for 24 chromosomes, and the chromosome names as the row names
reads_per_chromo <- read.table(input_file, header = FALSE, nrows = 24, row.names = 1)

# 5. Extract values from column 1 (sequence length of ref sequence) and column 2 (number of mapped reads) of the idxstats file
ref_lengths <- c(as.numeric(reads_per_chromo[,1]))
mapped_reads <- c(as.numeric(reads_per_chromo[,2]))

# 6. Perform a linear regression to test if the numbers of sequenced and mapped reads on each chromosome are correlated with the number of reference reads
LM <- lm(ref_lengths~mapped_reads)
p_value <- summary(LM)$coefficients[2, "Pr(>|t|)"]

# 7. Calculate sum of all reference sequence lengths and total number of mapped reads
total_ref_mapped <- sum(ref_lengths)
total_sample_mapped <- sum(mapped_reads)

# 8. Normalized R_x calculation
# a. Initialize a vector to hold the normalized ratios for each chromosome
normalized_ratios <- numeric(24)

# b. Divide the number of reads that map to each chromosome by the total number of reads that map
# Then, do that same calculation for the reference. Divide the number for the sample by that of the reference
for (i in 1:24) {
  normalized_ratios[i] <- (reads_per_chromo[i, 2]/total_sample_mapped) / (reads_per_chromo[i, 1]/total_ref_mapped)
}

# c. Divide the ratio for the X chromosome by the ratio of every other autosome
normalizedchromoX <- normalized_ratios[23] / normalized_ratios[1:22]

# d. Calculate R_x and 95% confidence interval upper/lower bounds if the number of reads is above 1,000 and the p value is significant
rx <- NA
rx_chromosex <- NA
rx_confidence_interval <- NA
rx_lowerbound <- NA
rx_upperbound <- NA

if (total_sample_mapped > 1000 && p_value < 0.05) {
  rx <- mean(normalizedchromoX)
  rx_confidence_interval <- 1.96 * (sd(normalizedchromoX)/sqrt(22))
  rx_lowerbound <- rx - rx_confidence_interval
  rx_upperbound <- rx + rx_confidence_interval
}

# 9. R_y calculation
ry <- NA
ry_chromosex <- NA
ry_confidence_interval <- NA
ry_lowerbound <- NA
ry_upperbound <- NA

# a. Calculate number of reads reads mapping to the X and Y chromosomes
reads_mapping_X <- reads_per_chromo[23, 2]
reads_mapping_Y <- reads_per_chromo[24, 2]
reads_mapping_sex_chromos <- reads_mapping_X + reads_mapping_Y

# b. Calculate R_y and 95% confidence interval upper/lower bounds if the number of reads is above 100,000
if (total_sample_mapped > 100000) {
  ry <- reads_mapping_Y / reads_mapping_sex_chromos
  
  # b. calculate standard error and 95% confidence interval
  ry_se <- sqrt((ry * (1.0 - ry)) / reads_mapping_sex_chromos)
  ry_confidence_interval <- 1.96 * ry_se
  ry_lowerbound <- ry - ry_confidence_interval
  ry_upperbound <- ry + ry_confidence_interval
}

# 10. Use upper and lower bound cutoffs for R_x to produce a chromosomal sex assignment for the individual based on R_x
# a. assign XX if the Rx CI lower and upper bounds are >0.8
if (!is.na(rx_lowerbound) && !is.na(rx_upperbound) && rx_lowerbound > 0.8 && rx_upperbound > 0.8) {
  rx_chromosex <- 'XX'
}
# b. assign XY if the Rx CI lower and upper bounds are <0.6
if (!is.na(rx_lowerbound) && !is.na(rx_upperbound) && rx_lowerbound < 0.6 && rx_upperbound < 0.6) {
  rx_chromosex <- 'XY'
}
# c. assign "consistent with XX but not XY" when the CI upper bound is >0.8 but the lower bound is <0.8
if (!is.na(rx_lowerbound) && !is.na(rx_upperbound) && rx_lowerbound < 0.8 && rx_upperbound > 0.8) {
  rx_chromosex <- 'consistent with XX but not XY'
}
# d. assign "consistent with XY but not XX" when the CI lower bound is <0.6 but the upper bound is >0.6
if (!is.na(rx_lowerbound) && !is.na(rx_upperbound) && rx_lowerbound < 0.6 && rx_upperbound > 0.6 && rx_upperbound < 0.8) {
  rx_chromosex <- 'consistent with XY but not XX'
}


# 11. Use upper and lower bound cutoffs for R_y to produce a chromosomal sex assignment for the individual based on R_y
# a. assign XX if the R_y CI lower and upper bounds are <0.016
if (!is.na(ry_lowerbound) && !is.na(ry_upperbound) && ry_lowerbound < 0.016 && ry_upperbound < 0.016) {
  ry_chromosex <- 'XX'
}
# b. assign XY if the R_y CI lower and upper bounds are >0.075
if (!is.na(ry_lowerbound) && !is.na(ry_upperbound) && ry_lowerbound > 0.075 && ry_upperbound > 0.075) {
  ry_chromosex <- 'XY'
}
# c. assign "consistent with XX but not XY" when the CI lower bound is <0.016 but the upper bound is >0.016
if (!is.na(ry_lowerbound) && !is.na(ry_upperbound) && ry_lowerbound < 0.016 && ry_upperbound > 0.016 && ry_upperbound < 0.075) {
  ry_chromosex <- 'consistent with XX but not XY'
}
# d. assign "consistent with XY but not XX" when the CI upper bound is >0.075 but the lower bound is <0.075
if (!is.na(ry_lowerbound) && !is.na(ry_upperbound) && ry_lowerbound < 0.075 && ry_upperbound > 0.075) {
  ry_chromosex <- 'consistent with XY but not XX'
}


# 12. Create a new variable to compare whether sex estimations from R_y and R_x match
if (!is.na(rx_chromosex) && !is.na(ry_chromosex) && rx_chromosex == ry_chromosex) {
  # if sex estimation from R_y and R_x match, assign R_x sex estimation (chosen arbitrarily) to new consensus sex estimation variable
  consensus_chromosex <- rx_chromosex
} else {
  consensus_chromosex <- 'no consensus'
}

# 13. Print data to console, only print header if the file does not yet exist
output_file <- "chromosex_summary.txt"

if (!file.exists(output_file)) {
  cat("File_Name\tTotal_Mapped_Reads\tReads_Mapped_to_Sex Chromosomes\tReads_Mapped_to_X\tReads_Mapped_to_Y\tR_y\tRy_Lower_95_CI_Bound\tRy_Upper_95_CI_Bound\tR_y_Chromosex\tR_x\tRx_Lower_95_CI_Bound\tRx_Upper_95_CI_Bound\tP_Value\tR_x_Chromosex\tConsensus_Chromosex\n", file = output_file)
}
  
  cat(base_file_name, total_sample_mapped, reads_mapping_sex_chromos, reads_mapping_X, reads_mapping_Y,
      round(ry, digits = 4), round(ry_lowerbound, digits = 4), round(ry_upperbound, digits = 4), ry_chromosex, 
      round(rx, digits = 4), round(rx_lowerbound, digits = 4), round(rx_upperbound, digits = 4), p_value, rx_chromosex, 
      consensus_chromosex, "\n", file = output_file, append = TRUE, sep = "\t")
