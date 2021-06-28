import csv
import os
import pandas as pd

def main():
    # needs to be in the format of 'run4June2021'
    file_date = 'run4June2021'
    file_directory = '/home/lilia/beebiome/beebiome-scripts/xml_parse_output_' + file_date + '/'
    synonyms = pd.read_csv(file_directory + 'apoidea-synonyms.csv')
    synonyms['count'] = [int(x) for x in synonyms['count'].fillna(0)]
    # I'm not sure if the taxid field can ever be a float but heres a line just in case which would 
    # handle the floats into int operation. Uncomment it if you want.
    # synonyms['taxid'] = [int(x) for x in synonyms['taxid']]
    for f_name in os.listdir(file_directory):
        hosts = pd.read_csv(file_directory + f_name)['host']
        for synonym_index, synonym in enumerate(synonyms['host']):
            synonyms['count'][synonym_index] = int(synonyms['count'][synonym_index]) + list(hosts).count(synonym)
        synonyms.to_csv(file_directory + 'apoidea-synonyms.csv', index=False)

if __name__ == "__main__":
    main()