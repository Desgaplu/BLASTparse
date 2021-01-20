#! python3
# -*- coding: utf-8

"""
BLASTparse.py
Created on Mon Dec 14 14:54:15 2020

Descrition:
Parse a multiple JSON BLAST output for the best match of each sequence and
uses the fasta sequences of the samples to classify them in folder names with
the closest genus identified with the BLAST search.

Output a TSV file containing the sample name, samples seq length,
closest species, %identity, species accession number and the sample fasta
sequence for each samples.
It then filter the TSV file to create a folder for each genus and add a fasta
file containing all sample sequences belonging to that genus.

@author: pldesgagne
"""

import tkinter as tk  # using the open filedialog
from tkinter import filedialog  # using the open filedialog
import ntpath  # get the filename from path string
import json
import os
import pandas as pd
from Bio.Align.Applications import MuscleCommandline

# Muscle.exe need to be in the same folder as this script
# from http://www.drive5.com/muscle/downloads.htm

# MEGACC (MEGA command console) need to be installed on computer
# The script use megacc command line from os.system() to create NJ trees
# The parameters of the tree creation are saved in file infer_NJ_nucleotide.mao
# This file is created with MEGAX in prototype mode

"""
TODO:
    -Ajust parameters for the selection of closest match
        -Current:   -Unique species name
                    -At least 5 hits
                    -Stop when %identity < 95%
                    -Max of 20 hits

    -Ajust length (current:300bp) for what is considered a short sequence

    -***Add an option to change parameters
    
    -***Add a function that create a word document with genus figure names

    -***Create a csv file with the species count
    
    -Create better instructions with a pause before the files query
    
DONE:
    -Add output options:
        -Summary only
        -Fasta sequences separated by genus
        -Fasta sequences separated by genus aligned with closest TYPE strain +
            newick NJ tree output.
        -Add the species rename function similar to MEGAparse.py script
        
    -Add percentage of identity on the 5 closest matches
    
    -Convert the script to execute all functions at same level 
        for better readability
"""


# Initialise tkinter to enable the uses of filedialog
root = tk.Tk()
root.withdraw()  # Prevent a empty window to be opened


def get_fasta_dict(fasta_data):
    """ Return a dictionnay of samples names and their fasta sequence """
    fasta_list = fasta_data.strip('\n').split('\n')
    fasta_seqs = {}
    for i in range(1, len(fasta_list), 2):
        fasta_seqs[fasta_list[i-1].strip('>')] = fasta_list[i-1] + '\n' + fasta_list[i]

    return fasta_seqs


def parse_all_json(data_list, dirname, fasta_seqs):
    """
    Fetch the BLAST results information for each samples (1 per JSON file)

    Parameters
    ----------
    data_list : JSON
        contains a list of JSON file names.
    dirname : STR
        Path of the JSON files.

    Returns
    -------
    results : list of lists
        Each list contains the BLAST info of each samples.
        [name, length, species, identity, accession, fasta_seq, matches]

    """
    
    def find_best_matches(data):
        """
        Finds the best matches according to the following parameters
        - unique species name
        - at least 5 results no matter percentage
        - stop when below 94%
        - stop at 20 results
    
        Return the list of best matches fasta seq with a list of corresponding % identity
        """
        species_list = [] # all species seen
        counter = 0
        matches = {}
        matches_percentage = {}
        for hit in data['BlastOutput2']['report']['results']['search']['hits']:
            # Ignore if species has already been added
            species = hit['description'][0]['sciname']
            if species in species_list:
                continue
    
            # query_to = hit['hsps'][0]['query_to']
            identity = hit['hsps'][0]['identity']
            align_len = hit['hsps'][0]['align_len']
            per_identity = identity/align_len
            # cover = align_len/query_to
            # if cover > 1:
            #     cover = 1
    
            # Conditions check
            if (counter >= 5 and per_identity < 0.95) or counter >= 20:
                break
    
            # Add species fasta seq in matches
            #title = hit['description'][0]['title']
            acc_num = hit['description'][0]['accession']
            # Original fasta format
            #matches[acc_num] = '>' + acc_num + ' ' + title + '\n' + hit['hsps'][0]['hseq']
            # Clean format
            matches[acc_num] = '>' + species + ' ' + acc_num + '\n' + hit['hsps'][0]['hseq']
            matches_percentage[acc_num] = f'{species} {per_identity:.2%}'
            species_list.append(species)
            counter += 1
        
        matches_percentage_ajusted = '\n'.join(matches_percentage.values())
    
        return matches, matches_percentage_ajusted
    
    
    # Loops all file and extract best match from each
    results = []
    for file in data_list['BlastJSON']:
        # load file
        #print(f"{dirname}/{file['File']}")
        with open(f"{dirname}/{file['File']}") as json_file:
            data = json.load(json_file)
        # extract information
        name = data['BlastOutput2']['report']['results']['search']['query_title']
        length = data['BlastOutput2']['report']['results']['search']['query_len']
        identity = data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'][0]['identity']
        align_len = data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'][0]['align_len']
        per_identity =  identity/align_len
        accession = data['BlastOutput2']['report']['results']['search']['hits'][0]['description'][0]['accession']
        species = data['BlastOutput2']['report']['results']['search']['hits'][0]['description'][0]['sciname']
        query_to = data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'][0]['query_to']
        cover = align_len/query_to
        if cover > 1:
            cover = 1
        matches, matches_per = find_best_matches(data)
        # save information
        results.append([name, length, species, f'{per_identity:.2%}', f'{cover:.2%}', accession, fasta_seqs[name], matches, matches_per])

    return results


def create_dir_files(dirname, df, do_alignment):
    """
    Create 1 folder per unique genus which contain a fasta file with all
    sample sequences belonging to that genus
    
    If do_alignment, the closest matches are added to the fasta file and
    an aligned fasta file is saved separatly. Also, a newick tree of the
    alignment is saved.
    """
    
    
    def muscle_align(in_file):
        """
        Align a fasta file (in_file) with MUSCLE, output a aligned fasta file
        Then it uses the aligned fasta file to create a NJ tree
        The tree is saved in newick tree file
        """
        
        
        def newick_tree(out_file):
            """
            Uses the aligned fasta file to create a NJ tree
            The tree is saved in newick tree file
            """
            
            print(f'NJ tree of {out_file[:-4]}_Tree.nwk')
            nj_cline = f'megacc -a infer_NJ_nucleotide.mao -d "{out_file}" -o "{out_file[:-4]}_Tree.nwk" -n' # Add -n to remove summary
            os.system(nj_cline)
            
            # stream = os.popen(f'megacc -a infer_NJ_nucleotide.mao -d "{out_file}" -o "{out_file[:-4]}.nwk"')
            # output = stream.read()
            # print(output)
            
            
        #Align sequences with MUSCLE
        muscle_exe = "muscle.exe"
        out_file = f'{in_file[:-4]}_aligned.fas'
        # Creating the command line
        muscle_cline = MuscleCommandline(muscle_exe,
                                         input=in_file,
                                         out=out_file)
        # Lauch the muscle command line
        print(f'MUSCLE sequences alignment of {out_file}')
        muscle_cline()
    
        # Create a NJ newick  tree
        newick_tree(out_file)
    
    
    try:
        # Create a genus field
        df['Genus'] = df['Species'].str.split().apply(lambda x: x[0])
        # Create the main folder
        os.mkdir(f'{dirname}/Genus/')
        # For each genus
        for genus in df['Genus'].unique():
            # Make a folder for that genus
            os.mkdir(f'{dirname}/Genus/{genus}')
            print(f'{dirname}/Genus/{genus} folder created.')
            # Filter all sequences of that genus
            dflong = df.loc[(df['Genus'] == genus) & (df['Length'] > 300)]
            dfshort = df.loc[(df['Genus'] == genus) & (df['Length'] <= 300)]
            seq_genus = dflong['Sequence']
            seq_genus_short = dfshort['Sequence']
            # Add matches for all samples of same genus in a dictionnary
            # Duplicates will overwrite each other since acc_num is used as unique key
            type_strains_long = {k:v for d in dflong['Matches'] for k, v in d.items()}
            type_strains_short = {k:v for d in dfshort['Matches'] for k, v in d.items()}

            # Save sequences in a txt file with the genus name
            if len(seq_genus) > 0:
                # Save the sequences in fasta format for long sequence
                with open(f'{dirname}/Genus/{genus}/{genus}.fas', 'w') as file:
                    file.write('\n'.join(seq_genus.to_list()))
                    if do_alignment:
                        file.write('\n')
                        file.write('\n'.join(type_strains_long.values()))
                    print(f'{dirname}/Genus/{genus}/{genus}.fas file created.')
                # # Save the TYPE strain acc num in a txt file for long sequences
                # with open(f'{dirname}/Genus/{genus}/matches_acc_num.txt', 'w') as file:
                #     file.write('\n'.join(type_strains_long.keys()))
                #     print(f'{dirname}/Genus/{genus}/matches_acc_num.txt file created.')
                # MUSCLE alignment of the sequences, saving in a new fasta file
                if do_alignment:
                    muscle_align(f'{dirname}/Genus/{genus}/{genus}.fas')

            if len(seq_genus_short) > 0:
                # Save the sequences in fasta format for short sequence
                with open(f'{dirname}/Genus/{genus}/{genus}_short.fas', 'w') as file:
                    file.write('\n'.join(seq_genus_short.to_list()))
                    if do_alignment:
                        file.write('\n')
                        file.write('\n'.join(type_strains_short.values()))
                    print(f'{dirname}/Genus/{genus}/{genus}_short.fas file created.')
                # # Save the TYPE strain acc num in a txt file for short sequences
                # with open(f'{dirname}/Genus/{genus}/matches_short_acc_num.txt', 'w') as file:
                #     file.write('\n'.join(type_strains_short.keys()))
                #     print(f'{dirname}/Genus/{genus}/matches_short_acc_num.txt file created.')
                # MUSCLE alignment of the sequences, saving in a new fasta file
                if do_alignment:
                    muscle_align(f'{dirname}/Genus/{genus}/{genus}_short.fas')

    except OSError:
        print ("Creation of the directory or file failed. The script does not overwrite existing folders.")
    else:
        print('------------------------------------------------')
        print ("Successfully created all directories and files")


# ----------------------------------Main---------------------------------------

print('*Do not execute this script on server since it creates folders.*')
print('Locally unzip the multiple-files JSON BLAST result.')
print("Select the main JSON file; It's the file with no number at the end.\n")

# Ask for the main JSON file path
json_file_path = filedialog.askopenfilename(
    title="Select an JSON main File",
    filetypes=(("JSON files", "*.json"), ("All files", "*.*"))
    )


print('Select the FASTA file used for the BLAST search:')
# Ask for the fasta file path
fasta_file_path = filedialog.askopenfilename(
    title="Select an FASTA main File",
    filetypes=(("FASTA files", "*.fas"), ("All files", "*.*"))
    )

try:
    # Loading main JSON file
    with open(json_file_path) as json_file:
        data_list = json.load(json_file)

    # Loading FASTA sequences
    with open(fasta_file_path) as fasta_file:
        fasta_data = fasta_file.read()

    dirname = ntpath.dirname(json_file_path) # Name of the current folder
    
    # Dict for easy search by index name
    fasta_seqs = get_fasta_dict(fasta_data)
    
    # Parsing the BLAST multiple-file JSON results
    results = parse_all_json(data_list, dirname, fasta_seqs)

    # Make into dataframe for easier filtering options
    df = pd.DataFrame(results, columns=['Name', 'Length', 'Species',
                                        '%identity', '%cover', 'Acc Num',
                                        'Sequence', 'Matches', 'Matches Percentage'])
    # Save BLAST results dataframe into a csv file
    columns = ['Name', 'Species', '%identity', '%cover', 'Length', 'Matches Percentage']
    df[columns].to_csv(f'{dirname}/BLASTresults.csv', index=False)
    print(f'\n\tResults summary saved in {dirname}/BLASTresults.csv')


    print()
    print('Do you wish to classify the sequences in genus specific folders?')
    print('1) Fasta sequences in genus specific folders')
    print('2) Fasta sequences in genus specific folders aligned with closest TYPE strains and NJ tree')
    print('3) Change parameters (TODO)')
    print('x) Skip and Exit')
    query = input('Choice: ')
    if query.lower() == '1':
        # Create folders and files with fasta sequences, alignments and TYPE strain acc nums
        create_dir_files(dirname, df, do_alignment=False)
    elif query.lower() == '2':
        create_dir_files(dirname, df, do_alignment=True)


except FileNotFoundError:
    print()
    print('ERROR: File selection was cancelled.')

except PermissionError:
    print()
    print('ERROR: BLASTresults.csv file is already open. Please close file.')

input("Press Enter to exit.")