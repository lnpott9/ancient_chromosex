# shell script to be paired with chromosex.r for processing BAM files

# 1. copy all deduplicated bam files to input folder

# 2. designate input folder
input_folder="<path/to/input/folder>"

# 3. designate output file name
export output_summary_file="pialq_phaseii_chromosex_summary.txt"

# 4. load R and samtools
module load R 
module load samtools

# 5. loop through all bam files and convert them to idxstats, then run the idxstat file through the R script
for bam_file in "$input_folder"/*.bam; do
	bai_file="${bam_file%.bam}.bai"
	if [ ! -f "$bai_file" ]; then
		samtools index "$bam_file"		
	fi
	output_file="${bam_file%.bam}.idxstats"
	samtools idxstats "$bam_file" > "$output_file"
	Rscript chromosex.r "$output_file"
done

# 6. remove intermediate idxstats files
rm $input_folder/*idxstats
