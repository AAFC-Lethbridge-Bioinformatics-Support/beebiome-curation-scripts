import csv
import xmltodict
import os
import requests
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np

# rearrange columns and have sample accession first
# make a beta column

def selectResults(biosamples_path):
    biosamples = pd.read_csv(biosamples_path)
    results = []

    for i, row in biosamples.iterrows():
        result = {}
        result['sample_name'] = biosamples.at[i,'accession']
        result['project_name'] = biosamples.at[i, 'has_proj']
        result['biosample_accession'] = biosamples.at[i,'accession']
        result['sra_accession'] = biosamples.at[i, 'has_sra']
        results.append(result)
    return results

def savetoCSV(data, filename):
    fields = data[0].keys()
    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = fields)
        writer.writeheader()  
        writer.writerows(data)

def main(): 
    for f_name in os.listdir('/home/lilia/beebiome/beebiome-scripts/xml_parse_output_run1June2021/'):
        if (f_name.startswith('Apoidea_biosample_')):
            print(f_name + '.csv is processing...')
            data = selectResults('/home/lilia/beebiome/beebiome-scripts/xml_parse_output_run1June2021/' + f_name)
            #print(data)
            savetoCSV(data, '/home/lilia/beebiome/beebiome-scripts/xml-parse-output/' + f_name + '_upload.csv') 
        
if __name__ == "__main__":
  
    # calling main function
    main()
