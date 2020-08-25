snaptools align-paired-end  \
	--input-reference=/nas01/Genome/bwaindex/Macaca_mulatta_revised.fa  \
	--input-fastq1=ATAC-110225_demultiplexed_R1.fq.gz  \
	--input-fastq2=ATAC-110225_demultiplexed_R2.fq.gz  \
	--output-bam=ATAC-library-110225.bam  \
	--aligner=bwa  \
	--path-to-aligner=/home/devdata/Software/soft/anaconda3/bin/  \
	--read-fastq-command=zcat  \
	--min-cov=0  \
	--num-threads=5  \
	--if-sort=True  \
	--tmp-folder=./  \
	--overwrite=TRUE

# fetchChromSizes /Macaca_fascicularis_5.0.dna.toplevel_release91.fa > mf.chrom.sizes
# samtools faidx input.fa
# cut -f1,2 /nas01/Genome/bwaindex/Macaca_mulatta_revised.fa.fai > mm80.chrom.sizes

snaptools snap-pre  \
	--input-file=ATAC-library-110225.bam  \
	--output-snap=ATAC-library-110225.snap  \
	--genome-name=hg38  \
	--genome-size=mm80.chrom.sizes  \
	--min-mapq=10  \
	--min-flen=0  \
	--max-flen=1000  \
	--keep-chrm=TRUE  \
	--keep-single=FALSE  \
	--keep-secondary=FALSE  \
	--overwrite=True  \
	--min-cov=10  \
	--verbose=True

snaptools snap-add-bmat	\
	--snap-file=ATAC-library-110225.snap	\
	--bin-size-list 5000 10000 50000 100000	\
	--verbose=True
