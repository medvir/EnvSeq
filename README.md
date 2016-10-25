# Envelope Sequencing
Package (or how does one call this?) for HIV-1 env plasmid sequencing analysis with optim_assembly. Uses pcDNA3 as a backbone vector and HXB2 as refernce strain. 

## EnvSeq.sh
`EnvSeq.sh` is the main script for HIV-1 env plasmid sequencing analysis. Uses the backbone of the cloning plasmid to subtract reads.
Arguments:
1. MiSeq runfolder with gzipped fastq files
2. plasmid backbone to be excluded in fasta format
3. reference sequence in fasta format

## optimassembly.py
`optimassembly.py` is from ozagrodi

`pcDNA3_bb.fasta` is the reference for the plasmid (backbone only!)

`HXB2.fasta` is the reference for HIV-1
