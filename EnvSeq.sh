#!/bin/bash

### defaults
script_dir=$( dirname "$(readlink -f "$0")" )
backbone=${script_dir}/pcDNA3_bb.fasta
reference=${script_dir}/HXB2.fasta
reads_limit=100000
expected_length=3500

### arguments
if [ $# == 0 ]; then
	echo
	echo 'EnvSeq.sh [options] ...'
	echo
	echo '-r, --reference        refernce (default HXB2.fasta)'
	echo '-b, --backbone         plasmid backbone to be subtracted (default pcDNA3_bb.fasta'
	echo '-n, --reads_limit      limit number of reads (default 100000)'
	echo '-l, --expected_length  expected insert length (default 3500)'
	echo
	exit
fi

while [[ $# -gt 1 ]]; do
	key="$1"
	case $key in
		-r|--reference)
		reference="$2"
		shift # past argument
		;;
		-b|--backbone)
		backbone="$2"
		shift # past argument
		;;
		-n|--reads_limit)
		reads_limit="$2"
		shift # past argument
		;;
		-l|--expected_length)
		expected_length="$2"
		shift # past argument
		;;
		*)
		# unknown option
		;;
	esac
	shift # past argument or value
done

if [[ -n $1 ]]; then
    sample_dir=$1
fi

backbone=$( readlink -f $backbone )
reference=$( readlink -f $reference )

echo
echo -e 'sample_dir \t' $sample_dir
echo -e 'script_dir \t' $script_dir
echo -e 'backbone \t' $backbone
echo -e 'reference \t' $reference
echo -e 'reads_limit \t' $reads_limit
echo -e 'expected_length ' $expected_length
echo

### index the plasmid backbone and the reference
smalt index -k 7 -s 2 plasmid $backbone
smalt index -k 7 -s 2 reference $reference

### loop over all fastq files
list=$(find $sample_dir | grep fastq.gz)

for i in $list; do
	
	echo sample $i
	
	### extract sample name, splits on "/" first and then on "_"
	arrIN=(${i//// })
	arrLEN=${#arrIN[@]}
	LAST=$((arrLEN - 1))
	FILENAME=${arrIN[$LAST]}
	ELS=(${FILENAME//_/ })
	sample=${ELS[0]}
	
	### create and change to director for each sample
	mkdir $sample
	cd $sample
		
	### limit number of reads
	seqtk sample ../$i $reads_limit > reads_sample.fastq

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
	python3.4 ${script_dir}/optimassembly.py -f reads_insert.fastq -r $reference -l $expected_length > consensus.fasta
	sed 's/NODE/'$sample'_optim/' consensus.fasta > ../${sample}_cons.fasta
		
	### de novo assembly velvet NOT IN USE
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