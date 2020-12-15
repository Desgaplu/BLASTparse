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
# Muscle was installed in PATH
# from http://www.drive5.com/muscle/downloads.htm

"""
TODO:
    -Ajust parameter for the selection of closest match
    -Ajust parameter for what is considered a short sequence
"""


# Initialise tkinter to enable the uses of filedialog
root = tk.Tk()
root.withdraw()  # Prevent a empty window to be opened


def find_best_matches(data):
    """
    Finds the best matches according to the following parameters
    - unique species name
    - at least 5 results no matter percentage
    - stop when below 94%
    - stop at 20 results
    
    list of fasta sequences
    """
    species_list = [] # all species seen
    counter = 0
    matches = {}
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
        species_list.append(species)
        counter += 1
    
    return matches


def muscle_align(in_file):
    #Align sequences with MUSCLE
    muscle_exe = "muscle.exe"
    out_file = f"{in_file[:-4]}_aligned.fas"
    # Creating the command line
    muscle_cline = MuscleCommandline(muscle_exe,
                                     input=in_file, 
                                     out=out_file)
    # Lauch the muscle command line
    print(f'MUSCLE sequences alignment of {out_file}')
    muscle_cline()


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
        [name, length, species, identity, accession, fasta_seq]

    """
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
        matches = find_best_matches(data)
        # save information
        results.append([name, length, species, per_identity, cover, accession, fasta_seqs[name], matches])
        #results.append(name + '\t' + str(length) + '\t' + species + '\t' + f"{identity:.4f}" + '\t' + accession + '\t' + fasta_seqs[name])
        
    return results


def get_fasta_dict(fasta_data):
    # Return a dictionnay of samples names and their fasta sequence
    fasta_list = fasta_data.strip('\n').split('\n')
    fasta_seqs = {}
    for i in range(1, len(fasta_list), 2):
        fasta_seqs[fasta_list[i-1].strip('>')] = fasta_list[i-1] + '\n' + fasta_list[i]
    
    return fasta_seqs


def create_dir_files(dirname, df):
    # Create a folder per unique genus which contain a txt file with all
    # samples sequences belonging to that genus
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
            type_strains_long = {k:v for d in dflong['Matches'] for k, v in d.items()}
            type_strains_short = {k:v for d in dfshort['Matches'] for k, v in d.items()}
            
            # Save sequences in a txt file with the genus name
            if len(seq_genus) > 0:
                # Save the sequences in fasta format for long sequence
                with open(f'{dirname}/Genus/{genus}/{genus}.fas', 'w') as file:
                    file.write('\n'.join(seq_genus.to_list()))
                    file.write('\n')
                    file.write('\n'.join(type_strains_long.values()))
                    print(f'{dirname}/Genus/{genus}/{genus}.fas file created.')
                # Save the TYPE strain acc num in a txt file for long sequences
                with open(f'{dirname}/Genus/{genus}/matches_acc_num.txt', 'w') as file:
                    file.write('\n'.join(type_strains_long.keys()))
                    print(f'{dirname}/Genus/{genus}/matches_acc_num.txt file created.')
                # MUSCLE alignment of the sequences, saving in a new fasta file
                muscle_align(f'{dirname}/Genus/{genus}/{genus}.fas')
                    
            if len(seq_genus_short) > 0:
                # Save the sequences in fasta format for short sequence
                with open(f'{dirname}/Genus/{genus}/{genus}_short.fas', 'w') as file:
                    file.write('\n'.join(seq_genus_short.to_list()))
                    file.write('\n')
                    file.write('\n'.join(type_strains_short.values()))
                    print(f'{dirname}/Genus/{genus}/{genus}_short.fas file created.')
                # Save the TYPE strain acc num in a txt file for short sequences
                with open(f'{dirname}/Genus/{genus}/matches_short_acc_num.txt', 'w') as file:
                    file.write('\n'.join(type_strains_short.keys()))
                    print(f'{dirname}/Genus/{genus}/matches_short_acc_num.txt file created.')
                # MUSCLE alignment of the sequences, saving in a new fasta file
                muscle_align(f'{dirname}/Genus/{genus}/{genus}_short.fas')
                
    except OSError:
        print ("Creation of the directory or file failed. The script does not overwrite existing folders.")
    else:
        print('------------------------------------------------')
        print ("Successfully created all directories and files")


# --------------------------------------------------------------------------

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
                                        'Sequence', 'Matches'])
    # Save BLAST results dataframe into a csv file
    df.to_csv(f'{dirname}/OUTPUT.csv', index=False)
    print(f'\n\tResults saved in {dirname}/OUTPUT.csv')
    
    
    print()
    # Create folders and files with fasta sequences, alignments and TYPE strain acc nums
    query = input('Do you wish to classify the sequences in genus specific folders? (y/n)')
    if query.lower() == 'y':
        create_dir_files(dirname, df)
    
except FileNotFoundError:
    print()
    print('ERROR: File selection was cancelled.')
    
except PermissionError:
    print()
    print('ERROR: OUTPUT file is already open. Please close file.')
    
input("Press Enter to exit.")