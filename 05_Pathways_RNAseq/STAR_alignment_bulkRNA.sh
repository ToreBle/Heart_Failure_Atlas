#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=64G
#SBATCH --job-name=STAR
#SBATCH --time=24:00:00
#SBATCH --output="STAR_output_%A_%a.out"
#SBATCH --array=1-6

export PATH=/home/uv525372/programmes/FastQC/:$PATH
export PATH=/home/uv525372/STAR-2.7.11b/source:$PATH
export PATH=/home/uv525372/samtools-1.16.1/bin:$PATH
export PATH=/home/uv525372/programmes/htslib:$PATH
export PATH=/home/uv525372/programmes/bcftools/:$PATH

# Sample information
input="./heart-atlas/bulkRNA/Cardiomyocytes/analysis/samples.txt"
IDs=`awk -v id="${SLURM_ARRAY_TASK_ID}" 'NR==id{print $1}' $input`
SAMPLE=$(echo "$IDs" | awk '{print $1}')

# Path information
PATH_FASTQ=$"/data/Kramann_lab/RawData/bulkRNA"
FASTQ_1="$PATH_FASTQ"/"$SAMPLE"/*_1.fq.gz
FASTQ_2="$PATH_FASTQ"/"$SAMPLE"/*_2.fq.gz

OUT=$"./heart-atlas/bulkRNA/Cardiomyocytes"
OUT_DIR="$OUT"/"$SAMPLE"/"$SAMPLE"
GENOME_DIR="./genomes/hg38/"
SJDB_GTF="./genomes/Homo_sapiens.GRCh38.98.chr.gtf"

# Print sample
echo "$SAMPLE"
echo "$OUT"

# STAR alignment
STAR --runThreadN $SLURM_CPUS_PER_TASK \
     --genomeDir "$GENOME_DIR" \
     --sjdbGTFfile "$SJDB_GTF" \
     --readFilesIn <(gunzip -c $FASTQ_1) <(gunzip -c $FASTQ_2) \
     --outFileNamePrefix "$OUT_DIR" \
     --outSAMattributes MD NH

# Convert into BAM
samtools view -b $OUT_DIR"Aligned.out.sam" > $OUT_DIR".bam"
samtools sort $OUT_DIR".bam" -o $OUT_DIR"_sorted.bam"
samtools index $OUT_DIR"_sorted.bam"
samtools flagstat $OUT_DIR"_sorted.bam" > $OUT_DIR"_flagstat_alignment.txt"

# Mark duplicates with Picard
picard MarkDuplicates INPUT=$OUT_DIR"_sorted.bam" OUTPUT=$OUT_DIR"_deduplicated.bam" METRICS_FILE=$OUT_DIR".txt" REMOVE_DUPLICATES=true
samtools sort $OUT_DIR"_deduplicated.bam" -o $OUT_DIR"_sorted.bam"
samtools index $OUT_DIR"_sorted.bam"
samtools flagstat $OUT_DIR"_sorted.bam" > $OUT_DIR"_flagstat_deduplicate.txt"

# Clean up intermediate files
rm $OUT_DIR"_deduplicated.bam" $OUT_DIR".bam"

# Clean up files
rm -r $OUT_DIR"_STARgenome"
rm $OUT_DIR"SJ.out.tab" $OUT_DIR"Aligned.out.sam"