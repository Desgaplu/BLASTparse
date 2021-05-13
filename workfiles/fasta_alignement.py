# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 22:05:08 2020

@author: pldesgagne

Align Testing
"""

from Bio.Align.Applications import MuscleCommandline
# Muscle was installed in PATH
# from http://www.drive5.com/muscle/downloads.htm

#Align sequences with MUSCLE (using parameters to make the alignment
#process as fast as possible)
muscle_exe = r"muscle.exe"
in_file = r"Bacillus.fas"
out_file = r"Bacillus_aligned.fas"
muscle_cline = MuscleCommandline(muscle_exe,
                                 input=in_file, 
                                 out=out_file)
#muscle_cline()

#muscle_tree = MuscleCommandline(cmd='maketree', input=out_file, out='pedio_tree.nwk', cluster='neighborjoining' )
#muscle -maketree -in pediotest_aligned.fas -out pedio_tree.phy -cluster neighborjoining


# from Bio.Phylo.TreeConstruction import DistanceCalculator
# from Bio import AlignIO
# aln = AlignIO.read(out_file, 'fasta')
# print(aln)

# calculator = DistanceCalculator('identity')
# dm = calculator.get_distance(aln)

# from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
# constructor = DistanceTreeConstructor(calculator, 'nj')
# tree = constructor.build_tree(aln)

# from Bio.Phylo import write
# write(tree, "pediotree.nwk", 'newick')