# BLASTparse
Allow the identification of the species of many sequences and creates a phylogenetic tree with the closest species.

Parse the Single-file JSON BLAST results and filter the closest species for each sequences.
Creates a folder for each genus in which sequences of sample of that genus are saved.
Gathers the sequences from the same genus together and creates a phylogenetic tree.

## Installation

1. Install MEGA-CC software.
   -Check the installation of MEGACC with the command (Press Windows + R and type cmd): megacc --version.
2. Required libraries:
   -pandas
   -biopython
   -python-docx
   -openpyxl


## Instructions

1. From the NCBI website; BLAST all rDNA sequences in Fasta format at the same time.
   1. For 16S, select the 16S database from the rRNA/ITS Databases section.
   2. For 26S/28S, select the nucleotide collection from the Standard Database and check the 'Sequence from type material' button.
2. Download the Single-file JSON in the result section.
3. Execute BLASTparse.py
   1. Load in the fasta file and the JSON file.
   2. The script outputs BLAST_results.csv and a folder will all identified genera.
