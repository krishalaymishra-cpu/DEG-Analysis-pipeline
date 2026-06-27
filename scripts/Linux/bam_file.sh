### Step 1: create genome index using HISAT2 indexing tool using reference genome. the reference genome i used was Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
hisat2-build reference.fa genome_index

### Step 2: align for bam file generation using parallel processing
#!/bin/bash
#path to hisat2 index
index=/path/to/hisat_index/*.ht2
fastq_dir= path/to/trimmed_files/
out_dir=/mnt/c/users/bionet3/onedrive/desktop/gse138109/hisat_res
mkdir -p $out_dir

#loop through all _1 files
for R1 in $fastq_dir/*_1.trimmed.fastq.gz
do
        sample=$(basename $R1 _1.trimmed.fastq.gz)
        R2=$fastq_dir/${sample}_2.trimmed.fastq.gz
        echo "processing $sample .."

        #align hisat2 pair reads
        hisat2 -p 16 \
               -x $index \
               -1 $R1 -2 $R2 \
               -S $out_dir/${sample}.sam \
               2> $out_dir/${sample}.log
        #convert sam->bam
        samtools view -@ 16 -bS $out_dir/${sample}.sam | samtools sort
        samtools index $out_dir/${sample}.sorted.bam

        #optional cleanup
        rm $out_dir/${sample}.sam

done
