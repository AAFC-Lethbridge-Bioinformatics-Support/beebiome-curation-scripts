import csv
import xmltodict
import os
import requests
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np


def parseXML(xmlfile, beta_data):
    tree = ET.parse(xmlfile)

    root = tree.getroot()

    #print(root) ----> BioSampleSet at some memory location
    data = []
    number_of_samples = len(root)
    for document_index, document_summary in enumerate(root.findall('./')):
        bioproject = {}
        bioproject['uid'] = document_summary.get("uid")
                
        for attribute in document_summary:
            if attribute.tag == 'Project':
                bioproject['project_xml'] = xmltodict.parse(ET.tostring(attribute).decode('utf-8'))
                for id in attribute.findall('./ProjectID/'):
                    if id.tag == 'ArchiveID':
                        bioproject['project_accession'] = id.get("accession")
                        if bioproject['project_accession'] in beta_data.loc[:, "BioProject acc"].values:
                            bioproject['in_beta'] = 'Yes'
                        else:
                            bioproject['in_beta'] = 'No'
                        break
                for target in attribute.findall('./ProjectType/ProjectTypeSubmission/Target/'):
                    if target.tag == 'Organism':
                        bioproject['tax_id'] = target.get("taxID")
                        break
                else:
                    bioproject['tax_id'] = 'Not found'
                for descr in attribute.findall('./ProjectDescr/'):
                    if descr.tag == 'Publication':
                        bioproject['pmid'] = descr.get("id")
                        break
                else: 
                    bioproject['pmid'] = 'Not found'
            else:
                bioproject['last_update'] = attribute.get("last_update")
                bioproject['submission_id'] = attribute.get("submission_id")
                bioproject['submitted'] = attribute.get("submitted")
        data.append(bioproject)
        #percent = (sample_index/number_of_samples * 100) 
        #print('Done ' + str(sample_index + 1) + '/' + str(number_of_samples) + ' samples or ' + str(percent) + '% of the set.')
    return data

def savetoCSV(data, filename):
    fields = ['project_accession', 'tax_id', 'pmid', 'in_beta', 'uid', 'last_update', 'submission_id', 'submitted', 'project_xml']
    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = fields)
        writer.writeheader()  
        writer.writerows(data)

def main(): 
    beta_data = pd.read_csv(r'/home/lilia/beebiome-taxonomy-scripts/old-site-data.csv')
    for f_name in os.listdir('/home/lilia/beebiome/beebiome-update/data/Apoidea(copy)/'):
        if f_name.startswith('Apoidea_bioproject.') and (not (f_name.startswith('Apoidea_bioproject.12') 
        or f_name.startswith('Apoidea_bioproject.2'))):
            print(f_name + '.csv is processing...')
            data = parseXML('/home/lilia/beebiome/beebiome-update/data/Apoidea(copy)/' + f_name, beta_data)
            savetoCSV(data, '/home/lilia/beebiome/beebiome-scripts/results (march)/' + f_name + '.csv')    
            df = pd.read_csv(r'/home/lilia/beebiome/beebiome-scripts/results (march)/'+ f_name + '.csv')
            new_df = df.reindex(columns=['project_accession', 'tax_id', 'pmid', 'in_beta', 'uid', 'last_update', 'submission_id', 'submitted', 'project_xml'])
            new_df.to_csv(r'/home/lilia/beebiome/beebiome-scripts/results (march)/'+ f_name + '.csv', index=False)        

if __name__ == "__main__":
  
    # calling main function
    main()

