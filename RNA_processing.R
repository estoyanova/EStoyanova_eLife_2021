# QC
fastqc R1.fastq R2.fastq

# Alignment with STAR
STAR --genomeDir genom_sj --readFilesIn $fastq1 $fastq2 --runThreadN 32 --outFileNamePrefix PREFIX --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --genomeLoad LoadAndKeep --limitBAMsortRAM 20000000000 --outBAMsortingThreadN 6 --readFilesCommand zcat --outFilterMismatchNmax 999 --alignMatesGapMax 1000000 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 60 --outFilterMismatchNoverLmax 0.05

# QC with Picard
java -jar CollectRnaSeqMetrics.jar I=input.bam O=output_metrics.txt REF_FLAT=Picard_mm10_refFlat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=Picard_mm10_rrna_intervalList.txt &


# Transcript quantification with Salmon 

# generate transcriptome
library("BSgenome.Mmusculus.UCSC.mm10")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("org.Mm.eg.db")
allT <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene,columns=c("gene_id","tx_name"))
myT <- getSeq(BSgenome.Mmusculus.UCSC.mm10, allT)
names(myT) <-allT$tx_name
writeXStringSet(myT, "R_ucsc_transcriptome.fa")

# Salmon

for fn in rnaseq_fastq/*;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
~/Salmon-latest_linux_x86_64/bin/salmon quant -i R_ucsc_transcriptome_index -l A \
-1 ${fn}/${samp}_R1_001.fastq.gz \
-2 ${fn}/${samp}_R2_001.fastq.gz \
-p 32 -o quants_R_ucsc/${samp}_quant
done 


#### Reading in quant.sf files from Salmon output. 
directory = "~/quants_R_ucsc"
files = file.path(directory, list.files(directory), "quant.sf")
names(files) = paste0(list.files(directory))

#### Need to create a key between UCSC transcript ID and gene SYMBOL
k = keys(TxDb.Mmusculus.UCSC.mm10.knownGene, keytype = "TXNAME")
tx2gene = biomaRt::select(TxDb.Mmusculus.UCSC.mm10.knownGene, k, "GENEID", "TXNAME")
tx2gene = na.omit(tx2gene) # must remove NA otherwise gene level summarization fails


#### Use tximport to do transcript level summarization

library("tximport")
txi = tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")

#### Using DESeq2 for differential gene expression analysis. Rough scheme of experiment design. I've omitted particular details as the study is still unpublished.  

library(DESeq2)
sampleTable = data.frame(condition = c( "...", each = 2),
                         celltype = c( "...", each = 2))
rownames(sampleTable) = colnames(txi$counts)
dds = DESeqDataSetFromTximport(txi, sampleTable, design = ~condition + celltype + condition:celltype)
dds$group = factor(paste0(dds$celltype, dds$condition))
design(dds) = ~ group
dds = estimateSizeFactors(dds)
dds = DESeq(dds, betaPrior=TRUE) 
