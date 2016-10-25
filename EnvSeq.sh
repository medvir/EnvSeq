#!/bin/bash
# HIV-1 env plasmid sequencing analysis with optim_assembly.sh
# Uses the backbone of the cloning plasmid to subtract reads

### arguments
# 1. MiSeq runfolder with gzipped fastq files
# 2. plasmid backbone to be excluded in fasta format
# 3. reference sequence in fasta format

run=$1
bb=$2
ref=$3

### settings
reads_limit=100000
expected_length=3500

### index the plasmid backbone and the reference
smalt index -k 7 -s 2 plasmid $bb
smalt index -k 7 -s 2 reference $ref

### loop over all fastq files
list=$(find $run/Data/Intensities/BaseCalls | grep fastq.gz)
for i in $list; do
	
	echo sample $i
	### extract sample name, splits on "/" first and then on "_"
	arrIN=(${i//// })
	arrLEN=${#arrIN[@]}
	LAST=$((arrLEN - 1))
	FILENAME=${arrIN[$LAST]}
	ELS=(${FILENAME//_/ })
	sample=${ELS[0]}

	mkdir $sample
	cd $sample
		
	### limit number of reads
	seqtk sample $i $reads_limit > reads_sample.fastq

	### align against the plasmid backbone
	echo aligning reads to plasmid backbone
	smalt map -n 28 -x -y 0.9 -c 0.9 -f samsoft -o reads_plasmid.sam ../plasmid reads_sample.fastq

	### convert to bam and return only unmapped reads with samtools view -f 4
	samtools view -Su reads_plasmid.sam | samtools sort - reads_plasmid_sorted
	samtools view -h -f 4 reads_plasmid_sorted.bam > reads_insert.sam

	### convert to fastq
	seqret reads_insert.sam fastq::reads_insert.fastq

	### optim assembly
	echo optim assembly of insert reads
	python3.4 ../../../bin/optimassembly.py -f reads_insert.fastq -r ../$ref -l $expected_length > consensus.fasta
	sed 's/NODE/'$sample'_optim/' consensus.fasta > ../${sample}_cons.fasta
		
	### de novo assembly velvet
	#echo velvet assembly of insert reads
	#velveth $sample 29 -short -fastq reads_insert.fastq
	#velvetg $sample -min_contig_lgth 2000
	#sed 's/NODE/'$sample'_velvet/' $sample/contigs.fa >> ../${sample}_velvet_cons.fasta
	
	### remove temp files
	rm reads_sample.fastq
	rm reads_plasmid.sam
	rm reads_plasmid_sorted.bam
	rm reads_insert.sam
	rm reads_insert.fastq
	
	cd ../
done
