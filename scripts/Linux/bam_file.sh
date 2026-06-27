#!/bin/bash

### Step 2: Align reads and create sorted BAM

# HISAT2 index prefix (must match what you used in hisat2-build)
INDEX_PREFIX="path/to/store/genome"

# Input FASTQ files
READ1="path/to/sample_1_trimmed_new.fastq.gz"
READ2="path/to/sample_2_trimmed_new.fastq.gz"

# Output BAM file
SORTED_BAM= $outputdir

# Run HISAT2 and pipe directly to samtools
hisat2 -x $INDEX_PREFIX -1 $READ1 -2 $READ2 \
  | samtools view -b \
  | samtools sort -o $SORTED_BAM

# Index BAM file
samtools index $SORTED_BAM

echo "Alignment complete. Sorted BAM file: $SORTED_BAM"
