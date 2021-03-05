#!/usr/bin/env python
# coding: utf-8

# ### Cleaning the result table to make it report ready:
# Adding the percentages next to each species name and surrounding it with parenthesis.
# Creating a table with a genus count.



import pandas as pd


df = pd.read_csv("TotalTableCSV.csv")

# Adding percentages to species line that does have any
def addpercentage(x):
    if '%' not in x[0]:
        return x[0] + ' ' + x[1]
    else:
        return x[0]
    

spec1 = df[['Species', '%identity']].apply(addpercentage, axis=1)

# Splitting before adding parenthesis
spec2 = spec1.str.split('\n').apply(lambda x: [i.split() for i in x])


# Adding parenthesis around percentages
def addparenthesis(species):
    return [['('+item+')' if '%' in item else item for item in line] for line in species]
spec3 = spec2.apply(addparenthesis)

# Rejoining in a single string
spec_final = spec3.apply(lambda x: '\n'.join([' '.join(i) for i in x]))

# List of genus
genus = spec3.apply(lambda x: x[0][0])

df['Genus'] = genus
df['Species'] = spec_final

# Save to CSV file
df[['Name', 'Kingdom', 'Genus', 'Species']].to_csv("TotalTableCSV_MOD.csv")

# Save the genus count in a CSV file
df['Genus'].value_counts().to_csv('Genus count.csv')


# ### Genus list table (ordered the same way as the genus count)
# Also could've done it in excel... by removing the count column from 
#the 'genus count.csv' and saving it in a new file -_-
genus_count = genus.value_counts()
sorted_genus_count = genus_count.sort_index().sort_values(kind='mergesort', ascending=False)
sorted_genus_count.to_csv('Genus_list.csv')

