# **R_y and R_x Statistic Calculation Script**

Author: Laura Pott

## **Description**

This script calculates two statistics used to estimate chromosomal sex based on the number of reads that map to the X and Y chromosomes in a sample of sequencing reads:

1. **R_y Statistic**: Ratio of reads mapped to the Y chromosome divided by the sum of reads mapped to the X and Y chromosomes
   - Calculated following the method outlined in Skoglund et al. (2013)
   - Requires >100,000 total mapped reads
   - Reference: [DOI: 10.1016/j.jas.2013.07.004](https://doi.org/10.1016/j.jas.2013.07.004)
2. **R_x Statistic**: Ratio of normalized X chromosome reads to normalized autosomal reads
   - Calculated following the method described in Mittnik et al. (2016).
   - Requires >1,000 total mapped reads and a significant p-value from regression testing
   - Reference: [DOI: 10.1371/journal.pone.0163019](https://doi.org/10.1371/journal.pone.0163019)

The script outputs an assignment of **XX**, **XY**, **consistent with XX but not XY**, and **consistent with XY but not XX**.

**Note:** This script is not designed to detect intersex chromosomal combinations.



## **Requirements**

**Dependencies**

- `samtools`: Used for processing BAM files and generating `idxstats` files
- `R`: Required to run the R script for statistical calculations. All calculations are done in base R

**Environment**

- Unix/Linux terminal
- The script is designed to work on HPC systems where the `module` command is available. If this isn't the case, modify the shell script to call `samtools` and `R` differently. Ensure `samtools` and `R` are installed and in your `PATH`



## **Workflow**

### **Step 1: Prepare Input Data**

- Place all BAM files in a designated input folder

### **Step 2: Edit the Shell Script**

1. Open the shell script (`chromosex_shell.sh`) and specify the following variables:
   - `input_folder`: Path to the folder containing the BAM files
   - `output_summary_file`: Name of the summary output file
2. Ensure the shell script and R script (`chromosex.r`) is in the same directory as the input BAM files

### **Step 3: Run the Shell Script**

Run the shell script from the terminal

```
./chromosex_shell.sh
```

### **Step 4: Interpreting Output**

- The final output is a tab-delimited summary file with columns including:
  - File name
  - Total mapped reads
  - Mapped reads to X and Y chromosomes
  - Calculated R_y and R_x values
  - Chromosomal sex assignments from R_y, R_x, and consensus sex assignment
