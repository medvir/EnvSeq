# Envelope Sequencing
Package for HIV-1 env plasmid sequencing analysis with optimassembly.

## envseq.sh
`envseq.sh` is the main script for HIV-1 env plasmid sequencing analysis. First, it subtracts reads aligning to the backbone of the plasmid with smalt/samtools. Then it de novo assembles the insert using optimassembly. Default settings use pcDNA3 as a backbone vector and HXB2 as reference strain. Input is a directory containing gzipped fastq files. Output is a fast file SAMPLE_cons.fasta.

`optimassembly.py` is from [ozagordi](https://github.com/ozagordi)

`pcDNA3_bb.fasta` is the reference for the plasmid (backbone only!)

`HXB2.fasta` is the reference for HIV-1

### Usage
	usage: envseq.sh [options] ...
	
	OPTIONS
	-r, --reference		reference (default HXB2.fasta)
	-b, --backbone		plasmid backbone to be subtracted (default pcDNA3_bb.fasta)
	-n, --reads_limit	limit number of reads (default 100'000)
	-l, --expected_length	expected insert length (default 3500)
