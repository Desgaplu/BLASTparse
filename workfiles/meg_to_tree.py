# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 17:36:30 2020

@author: pldesgagne
"""

import os
cline_align = 'megacc -a muscle_align_nucleotide.mao -d Lactiplantibacillus.fas -o test.meg'
cline = 'megacc -a infer_NJ_nucleotide.mao -d test.meg -o testtree.nwk'

#os.system('megacc -h') #output 0 if work, 1 if doesnt

# stream = os.popen('megacc -h')
# output = stream.read()
# print(output)

os.system(cline_align)
os.system(cline)