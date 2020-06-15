# Envelope Sequencing
Script for HIV-1 env plasmid sequencing analysis.

## envseq.py
This script removes all sequences which map to the plasmid backbone given
(python variable plasmid_backbone), performs a denovo assembly using [VelvetOptimiser](https://github.com/tseemann/VelvetOptimiser)
and returns the sequence of the insert in the same orientation as the reference
sequence given (python variable plasmid_insert).

Caveats:  
- It's assumed that the sequence is clonal (no ambiguous positions)
- If multiple contigs are present, only the longest one is taken.
- If the identity of the pairwise sequence alignment to the given reference is below 50%, it's assumed that the contig is in reverse complement orientation, this could however also have other reasons.  

Regarding these points, it's therefore advised to look at the results (sequence lengths and maybe do a multiple sequence alignment) to spot obvious errors.

### Use conda environment from file
To ensure you have all dependencies needed for EnvSeq installed you can use the `environment.yml` file.  
First you need to have [Conda](https://conda.io/docs/install/quick.html) installed).  
With the command `conda env create -f <path>/environment.yml` you will create a copy of the EnvSeq environment.  
You enter the environment with the command `source activate EnvSeq` (and leave it with `source deactivate`).  
For more information visit following link to [Managing environments](https://conda.io/docs/using/envs.html).

### Usage
All files need to be in the same directory, this includes the reference sequences, sequencing reads and the envseq.py script.  
It can then be executed as follows:  
`python envseq.py`
