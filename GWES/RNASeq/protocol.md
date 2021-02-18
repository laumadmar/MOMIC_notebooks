# RNASeq protocol

1. Download reference genome 
`sudo wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`

2. Download annotation file
`wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz`

3. Align reads using STAR

    3.1 Generate genome index - once per genome type
`sudo STAR --runMode genomeGenerate --genomeDir STARindex --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.94.gtf -- sjdbOverhang 100 --runThreadN 15`

    3.2 Alignment. This step has to be done for each individual FASTQ file. If the fasta files are paired-end, combine them in a single call.
`STAR --genomeDir $path/STARindex --readFilesIn $files --readFilesCommand zcat --outFileNamePrefix $sample --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts --runThreadN $nthreads-1`

4. Quality control of aligned reads
- look at Log.final.out file, it contains alignment statistics. An alignment of RNA-seq reads is usually considered to have succeeded if the mapping rate (uniquely mapped reads) is >70%
- more statistics and plots with scripts/QCplotsSTAR.R

5. Read Quantification. Gene-based read counting
Done in point 3.2 adding to STAR parameter `--quantMode GeneCounts`. Second column of result file (ReadsPerGene.out.tab) coincide with result from htseq-count.


Following http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf