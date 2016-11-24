# Envelope Sequencing
Package for HIV-1 env plasmid sequencing analysis with optim_assembly.

## EnvSeq.sh
`EnvSeq.sh` is the main script for HIV-1 env plasmid sequencing analysis. Uses the backbone of the cloning plasmid to subtract reads.
Uses pcDNA3 as a backbone vector and HXB2 as refernce strain. 
Input: Directory with (gzipped) fastq files

`optimassembly.py` is from [ozagordi](https://github.com/ozagordi)

`pcDNA3_bb.fasta` is the reference for the plasmid (backbone only!)

`HXB2.fasta` is the reference for HIV-1

### Usage
	usage: EnvSeq.sh [options] ...
	
	OPTIONS
	-r, --reference			reference (default HXB2.fasta)
	-b, --backbone			plasmid backbone to be subtracted (default pcDNA3_bb.fasta
	-n, --reads_limit		limit number of reads (default 100'000)
	-l, --expected_length	expected insert length (default 3500)
