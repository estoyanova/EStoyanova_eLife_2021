
# QC and trimming of FastQ Files 
trim_galore --stringency 3 --fastqc --paired --clip_R1 5 --clip_R2 $R1.fastq $R2.fastq 

# Alignment with Bismark
bismark --bowtie2 -p 4 --multicore 4 bismark_genome_build/ -1 R1.fq -2 R2.fq 
deduplicate_bismark -p --bam inputfile.bam


# Sequencing quality control with CEGX custom
docker	run	-v=`pwd`:/Data	-it	cegx_bsexpress_0.6	auto_bsExpress

# Downstream processing with Methpipe from Smith lab
to-mr -o bismark.deduplicated.bam.mr -m bismark deduplicated.bam 
LC_ALL=C sort -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o out.mr.sorted  bismark.deduplicated.bam.mr
methcounts -c genome.fa -o output.allC.meth *.sorted.mr &
mlml -u *BS*.meth -m *OX*.meth -o *.mlml.txt


#Merging symmetric CpGs, necessary for MethylSeekR
symmetric-cpgs -o symmetric*.meth -v -m input.meth &
