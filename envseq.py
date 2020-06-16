import os
import glob

# all files need to be in the same directory,
# this includes the reference sequences (in a references folder), sequencing reads and this script

# backbone is used to remove reads mapping to it
plasmid_backbone = "references/pcDNA3_1.fasta"
# insert sequence is used to find the right orientation of the contigs
plasmid_insert = "references/HXB2_insert.fasta"


os.system(f"smalt index -k 7 -s 2 references/plasmid {plasmid_backbone}")

# run the following script for all fastq files present in the current directory
for file in glob.glob("*.fastq*"):
    # extract the sample name from the full fastq filename
    sample = file.split("_")[0]

    os.system(f"mkdir {sample}")

    # only sample maximum 1'000'000 reads to increase speed
    os.system(f"seqtk sample {file} 1000000 > {sample}/reads_sample.fastq")

    # map sampled reads against the plasmid backbone and only extract the remaining
    os.system(f"smalt map -n 28 -x -y 0.9 -c 0.9 -f samsoft -o {sample}/reads_plasmid.sam references/plasmid {sample}/reads_sample.fastq")
    os.system(f"samtools view -Su {sample}/reads_plasmid.sam | samtools sort -o {sample}/reads_plasmid_sorted.bam -")
    os.system(f"samtools view -h -f 4 {sample}/reads_plasmid_sorted.bam > {sample}/reads_insert.sam")
    os.system(f"seqret {sample}/reads_insert.sam fastq::{sample}/reads_insert.fastq")

    # removing bam and sam files to save storage
    os.system(f"rm {sample}/*am")

    # run a denovo assembly with all those remaining reads
    os.system(f"VelvetOptimiser.pl -s 140 -e 150  -t 28 -d '{sample}/denovo_assembly' -f '-short -fastq {sample}/reads_insert.fastq'")

    # take longest contig and change the fasta header to the sample name
    # if there are multiple contigs present, some information might be lost
    # it's therefore advised to check the length of all final sequences
    os.system(f"seqkit sort -l -r {sample}/denovo_assembly/contigs.fa | seqkit head -n 1 -w 0 | sed 's/>.*/>{sample}/' > {sample}/{sample}.fasta")

    # align contig to the reference to check if orientation is the same or not
    os.system(f"needle {plasmid_insert} {sample}/{sample}.fasta -gapopen 10 -gapextend 0.5 -outfile {sample}/{sample}.aln")

    alignment = open(f"{sample}/{sample}.aln", "r")
    alignment_lines = alignment.readlines()
    alignment.close()

    # extract the percent identity and gaps to determine the current orientation
    identity_line = alignment_lines[23]
    identity = identity_line[identity_line.find("(")+1:identity_line.find("%")]
    identity = float(identity)

    gaps_line = alignment_lines[25]
    gaps = gaps_line[gaps_line.find("(")+1:gaps_line.find("%")]
    gaps = float(gaps)

    # if the identity is below 50% it's assumed the orientation is reverse complement
    # this could of course also have different reasons, it's therefore advised
    # to create a multiple sequence alignment of all final sequences to verify
    # the result
    if identity < 50:
        os.system(f"seqkit seq -v -p -r -t 'DNA' -w 0 {sample}/{sample}.fasta > {sample}/{sample}_rc.fasta")
        os.system(f"cp {sample}/{sample}_rc.fasta {sample}.fasta")
    else:
        os.system(f"cp {sample}/{sample}.fasta {sample}.fasta")
