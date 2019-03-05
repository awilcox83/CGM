# CGM

Circular Genome Mapper is a tool to map sequencing reads from organisms with circular genomes to linear reference genomes.  A Bowtie2 installation is required.

CGM extends the reference genome, allowing reads that span both ends to map correctly.  Reads are then split into two and their positions adjusted.

Useage:

```
cgm.py -s [reads in FASTQ] -r [reference genome in FASTA] -o [prefix for output files]
```

```
cgm.py -1 [Fwd reads FASTQ] -2 [Rev reads FASTQ] -r [reference genome in FASTA] -o [prefix for output files]
```
