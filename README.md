# CGM

Circular Genome Mapper is a tool to map sequencing reads from organisms with circular genomes to linear reference genomes.

As yet, it only supports single-end reads.  Support for paired-end reads is under development.

A Bowtie2 installation is required.

Useage:

python3 cgm.py -s reads.FASTQ -r referenceGenome.FASTA -o cgm
