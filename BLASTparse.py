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

# Initialise tkinter to enable the uses of filedialog
root = tk.Tk()
root.withdraw()  # Prevent a empty window to be opened


def parse_all_json(data_list, json_path, fasta_seqs):
    """
    Fetch the BLAST results information for each samples (1 per JSON file)
    
    Parameters
    ----------
    data_list : JSON
        contains a list of JSON file names.
    json_path : STR
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
        #print(f"{json_path}/{file['File']}")
        with open(f"{json_path}/{file['File']}") as json_file:
            data = json.load(json_file)
        # extract information
        name = data['BlastOutput2']['report']['results']['search']['query_title']
        length = data['BlastOutput2']['report']['results']['search']['query_len']
        identity = data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'][0]['identity'] / data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'][0]['align_len']
        accession = data['BlastOutput2']['report']['results']['search']['hits'][0]['description'][0]['accession']
        species = data['BlastOutput2']['report']['results']['search']['hits'][0]['description'][0]['sciname']
        # save information
        results.append([name, length, species, identity, accession, fasta_seqs[name]])
        #results.append(name + '\t' + str(length) + '\t' + species + '\t' + f"{identity:.4f}" + '\t' + accession + '\t' + fasta_seqs[name])
        
    return results


def get_fasta_dict(fasta_data):
    # Return a dictionnay of samples names and their fasta sequence
    fasta_list = fasta_data.strip('\n').split('\n')
    fasta_seqs = {}
    for i in range(1, len(fasta_list), 2):
        fasta_seqs[fasta_list[i-1].strip('>')] = fasta_list[i-1] + '&' + fasta_list[i]
    
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
            seq_genus = df.loc[(df['Genus'] == genus) & (df['Length'] > 300)]['Sequence'].str.replace('&','\n')
            seq_genus_short = df.loc[(df['Genus'] == genus) & (df['Length'] <= 300)]['Sequence'].str.replace('&','\n')
            # Save sequences in a txt file with the genus name
            if len(seq_genus) > 0:
                with open(f'{dirname}/Genus/{genus}/{genus}.txt', 'w') as file:
                    file.write('\n'.join(seq_genus.to_list()))
                    print(f'{dirname}/Genus/{genus}/{genus}.txt file created.')
            if len(seq_genus_short) > 0:
                with open(f'{dirname}/Genus/{genus}/{genus}_short.txt', 'w') as file:
                    file.write('\n'.join(seq_genus_short.to_list()))
                    print(f'{dirname}/Genus/{genus}/{genus}_short.txt file created.')
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
    
    dirname = ntpath.dirname(json_file_path)
    fasta_seqs = get_fasta_dict(fasta_data) # Dict for easy search by index name
    results = parse_all_json(data_list, dirname, fasta_seqs)
    
    # Make into dataframe for easier filtering options
    df = pd.DataFrame(results, columns=['Name', 'Length', 'Species', 
                                        '%identity', 'Acc Num',
                                        'Sequence'])
    # Save dataframe into csv file
    df.to_csv(f'{dirname}/OUTPUT.csv', index=False)
    print(f'\n\tResults saved in {dirname}/OUTPUT.csv')
    
    
    print()
    query = input('Do you wish to classify the sequences in genus specific folders? (y/n)')
    if query.lower() == 'y':
        # Create folders and files with fasta sequences
        create_dir_files(dirname, df)
    
except FileNotFoundError:
    print()
    print('ERROR: File selection was cancelled.')
    
except PermissionError:
    print()
    print('ERROR: OUTPUT file is already open. Please close file.')
    
input("Press Enter to exit.")