#!/bin/bash
# HIV (env) plasmid sequencing analysis with optim_assembly.sh

# arguments
# 1. runfolder with gzipped fastq files
# 2. plasmid backbone to be excluded in fasta format
# 3. reference sequence in fasta format

list=$(find ${1}/Data/Intensities/BaseCalls | grep fastq.gz )

# variables
reads_limit=100000
exp_lgth=3000

# index the plasmid backbone and the reference
echo "*** indexing plasmid backbone and reference ***"
smalt index -k 7 -s 2 plasmid $2
smalt index -k 7 -s 2 reference $3


for i in $list; do
	
	echo "*** sample " $i " ***"
	# extract sample name
	# splits on "/" first and then on "_"
	arrIN=(${i//// })
	arrLEN=${#arrIN[@]}
	LAST=$((arrLEN - 1))
	FILENAME=${arrIN[$LAST]}
	ELS=(${FILENAME//_/ })
	sample=${ELS[0]}

	mkdir $sample
	cd $sample
		
	# sample fraction of reads
	seqtk sample $i $reads_limit > reads_limit.fastq

	# align against the plasmid (indexed already)
	echo "*** aligning reads to plasmid backbone ***"
	smalt map -n 28 -x -y 0.9 -c 0.9 -f samsoft -o reads_plasmid.sam ../plasmid reads_limit.fastq

	# convert to bam and return only unmapped reads with samtools view -f 4
	samtools view -Su reads_plasmid.sam | samtools sort - reads_plasmid_sorted
	samtools view -h -f 4 reads_plasmid_sorted.bam > reads_insert.sam

	# convert to fastq
	seqret reads_insert.sam fastq::reads_insert.fastq

	# optim assembly
	echo "*** optim assembly ***"
	python3.4 ../../../bin/optimassembly.py -f reads_insert.fastq -r ../$3 -l $exp_lgth > consensus.fasta
	
	# copy files and add sample to name of sequence in contigs.fasta and save in ../
	sed 's/NODE/'$sample'_NODE/' consensus.fasta > ../${sample}_consensus.fasta
	
	# remove temp files, comment out to debug
	#rm reads_limit.fastq
	rm reads_plasmid.sam
	rm reads_plasmid_sorted.bam
	#rm reads_insert.fastq
	
	cd ../

done
exit
