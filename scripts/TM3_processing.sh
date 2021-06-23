###
# TM3 with umi for mouse organoid
# download raw data from vbcf links
~/scripts/ngs_tools/download_ngs_links.sh URL_download.txt

mkdir -p ngs_raw/BAMs
mv *.bam ngs_raw/BAMs/

# bam2fastq
~/scripts/ngs_tools/bam2fastq.sh
~/scripts/ngs_tools/run_fastqc.sh

# umitools_extract 
~/scripts/rnaseq/umitools_extract.sh

# change folder names for the next steps
cd ngs_raw
mv FASTQs FASTQs_raw
mv FASTQs_umi/ FASTQs
cd ..

# tirmming adaptor for quant-seq
~/scripts/rnaseq/trimming_fastq.sh

# index hg38 genome for Star
#sbatch ~/scripts/RNAseq/index_hg38_star.sh

# mapping with star
~/scripts/rnaseq/mapping_star.sh -g hg38 -D ngs_raw/FASTQs_trimmed

# remove duplication using umitools
~/scripts/rnaseq/umitools_dedup.sh

# counting reads and umis for genes using only protein-coding gene annotation
~/scripts/rnaseq/htseq_counts.sh -G /groups/cochella/jiwang/Genomes/Human/hg38/annotation/Genes/genes.gtf
~/scripts/rnaseq/htseq_counts.sh -D /groups/cochella/jiwang/Projects/Ariane/R10331_quantseq_hg/BAMs_umi -G /groups/cochella/jiwang/Genomes/Human/hg38/annotation/Genes/genes.gtf

# multiqc 
ml load multiqc/1.7-foss-2018b-python-2.7.15
multiqc .
