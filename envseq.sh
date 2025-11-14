#!/bin/bash

### defaults
script_dir=$( dirname "$(readlink -f "$0")" )
backbone=${script_dir}/pcDNA3_bb.fasta
reference=${script_dir}/HXB2.fasta
sample_dir=/home/ubuntu/data/
reads_limit=ALL
expected_length=3500

### arguments
show_help() {
    echo
    echo "Usage: enveq.sh [options] [sample_dir]"
    echo
    echo "Options:"
    echo "  -r, --reference        Reference file (default: HXB2.fasta)"
    echo "  -b, --backbone         Plasmid backbone file to subtract (default: pcDNA3_bb.fasta)"
    echo "  -n, --reads_limit      Limit number of reads (default: ALL)"
    echo "  -l, --expected_length  Expected insert length (default: 3500)"
    echo "  -h, --help             Show this help message and exit"
    echo
    echo "If no arguments are given, defaults will be used."
	echo "If no sample directory is provided, /home/ubuntu/data/ will be used."
    echo
}

### handle help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
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
else
    echo "No sample directory provided — running with defaults."
fi

### convert relative to absolute path
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

	### create directory
	mkdir -p $sample

	### limit number of reads
	if [[ $reads_limit != "ALL" ]]; then
		seqtk sample $i $reads_limit > ${sample}/reads_sample.fastq
	else
		zcat "$i" > "${sample}/reads_sample.fastq"
	fi

	### change to directory
	cd $sample

	### align against the plasmid backbone
	echo aligning reads to plasmid backbone
	smalt map -n 28 -x -y 0.9 -c 0.9 -f samsoft -o reads_plasmid.sam ../plasmid reads_sample.fastq

	### convert to bam and return only unmapped reads with samtools view -f 4
	samtools view -Su reads_plasmid.sam | samtools sort -o reads_plasmid_sorted.bam -
	samtools view -h -f 4 reads_plasmid_sorted.bam > reads_insert.sam

	### convert to fastq
	seqret reads_insert.sam fastq::reads_insert.fastq

	### optim assembly
	echo optim assembly of insert reads
	python ${script_dir}/optimassembly.py -f reads_insert.fastq -r $reference -l $expected_length > consensus.fasta
	sed 's/consensus_contigs\|NODE/'$sample'_optim/' consensus.fasta > ../${sample}_cons.fasta

	### remove temp files
	rm reads_sample.fastq
	rm reads_plasmid.sam
	rm reads_plasmid_sorted.bam
	rm reads_insert.sam
	rm reads_insert.fastq

	cd ../
done
