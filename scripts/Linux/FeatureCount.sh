# after creating the bam file for piared ends, the feature count matrix was generated 
featureCounts -T 8 -p -a annotation.gtf -o counts.txt sample.sorted.bam

