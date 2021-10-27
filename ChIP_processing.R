

trim_galore --stringency 3 --fastqc --paired R1.fastq R2.fastq 

bowtie2 -p 64 -x /data1/genomes/Mus_musculus_UCSC_mm10/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 *.R1.fastq -2 *.R2.fastq -S output.sam > output.sam.out 2>&1&
  

  
# convert sam to bam
samtools view -@ 64 -bS SAM.sam -o SAM.bam

# sort bam file
samtools sort -@ 64 SAM.bam SAM.sorted

# filter duplicated reads
samtools rmdup $SAM.sorted.bam $SAM.sorted.filtered.bam


# index the sorted bam file
samtools index $SAM.secondsorted.filtered.bam

macs2 callpeak --nomodel -t ChIP.bam -c Control.bam --broad -g mm --broad-cutoff 0.1 -f BAMPE