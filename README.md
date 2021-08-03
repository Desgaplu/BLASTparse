# BLASTparse
Allow the identification of the species of many sequences and creates a phylogenetic tree with the closest species.

Parse the Multiple-file JSON BLAST results and filter the closest species for each sequences.
Creates a folder for each genus in which sequences of sample of that genus are saved.
Gathers the sequences from the same genus together and creates a phylogenetic tree.

1) From the NCBI website; BLAST all 16S rDNA sequences in Fasta format at the same time.
2) Download the Multiple-file JSON in the result section.
3) BLASTparse.py:
4) Load in the fasta file and the first JSON file.
5) The script outputs BLAST_results.csv and a folder will all identified genera.
