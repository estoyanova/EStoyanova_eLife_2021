
# QC and trimming of FastQ Files 
trim_galore --stringency 3 --fastqc --paired $R1.fastq $R2.fastq 

# Alignment of raw reads
#### Using bowtie2 
bowtie2 -p 64 -x Sequence/Bowtie2Index/genome -X 2000 --no-mixed --no-discordant -1 *.R1.fastq -2 *.R2.fastq -S output.sam > output.sam.out 2>&1&

# Post-alignment processing
samtools view -Sb -@ 10 in.sam -o out.bam
samtools sort -@ 64 out.bam out.sorted 
samtools rmdup out.sorted.bam out.sorted.rmdup.bam

#### Selection of subnucleosomal reads 
samtools view out.sorted.rmdup.bam | awk '(sqrt($9*$9)<100)' > out.sam
samtools view -H *.sorted.rmdup.bam > header.txt
cat header.txt out.sam | samtools view -Sb -o $100nt.bam

# BigWig for browser visualization

bamCompare -b1 ATAC.sorted.dupsrem.bam  -o sample.bw --normalizeUsingRPKM --minMappingQuality 30 --ignoreForNormalization chrX chrM chrY --blackListFileName ~/mm10.blacklist.bed -p 32 &

# Peak calling using macs2

macs2 callpeak --nomodel -t input.sorted.rmdup.100nt.bam -f BAMPE -n output_name --call-summits -g mm -B -q 0.05 --outdir outdir_name/ 

#### Removing peaks that map to blacklisted regions

bedtools intersect -v -a my.bed -b mm10.blacklist.bed | grep -v chrY | grep -v chrM > my.filtered.bed

