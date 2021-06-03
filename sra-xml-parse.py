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

def parseXML(xmlfile, beta_data):
    tree = ET.parse(xmlfile)

    root = tree.getroot()

    #print(root) ----> EXPERIMENT_PACKAGE_SET at some memory location
    # you can make a list to like collect stuff
    data = []
    number_of_samples = len(root)
    for experiment_package_index, experiment_package in enumerate(root.findall('./')):
        experiment = {}
        # bioproject['uid'] = document_summary.get("uid")
        # print(experiment_package) ----> EXPERIMENT_PACKAGE at some memory location     
        for attribute in experiment_package:
            # print(attribute) -----> EXPERIMENT, SUBMISSION, Organization, STUDY, SAMPLE, Pool, RUN_SET
            tag = attribute.tag 
            
            if tag == 'EXPERIMENT':
                experiment['experiment_accession'] = attribute.get("accession")
                experiment['experiment_alias'] = attribute.get("alias")
                experiment['experiment_xml'] = xmltodict.parse(ET.tostring(attribute).decode('utf-8'))
            elif tag ==  'SUBMISSION':
                experiment['submission_lab'] = attribute.get("lab_name")
                experiment['submission_center'] = attribute.get("center_name")
                experiment['submission_accession'] = attribute.get("accession")
                experiment['submission_alias'] = attribute.get("alias")
                experiment['submission_xml'] = xmltodict.parse(ET.tostring(attribute).decode('utf-8'))
            elif tag == 'Organization':
                experiment['organization_type'] = attribute.get("type")
                experiment['organization_xml'] = xmltodict.parse(ET.tostring(attribute).decode('utf-8'))
            elif tag == 'STUDY':
                experiment['study_center'] = attribute.get("center_name")
                experiment['study_alias'] = attribute.get("alias")
                experiment['study_accession'] = attribute.get("accession")
                experiment['study_xml'] = xmltodict.parse(ET.tostring(attribute).decode('utf-8'))
            elif tag == 'SAMPLE':
                experiment['sample_alias'] = attribute.get("alias")
                experiment['sample_accession'] = attribute.get("accession")
                if experiment['sample_accession'] in beta_data.loc[:, "BioSample acc"].values:
                    experiment['in_beta'] = 'Yes'
                else:
                    experiment['in_beta'] = 'No'
                experiment['sample_xml'] = xmltodict.parse(ET.tostring(attribute).decode('utf-8'))
            elif tag == 'Pool':
                experiment['pool_xml'] = xmltodict.parse(ET.tostring(attribute).decode('utf-8'))
            elif tag == 'RUN_SET':
                experiment['run_set_xml'] = xmltodict.parse(ET.tostring(attribute).decode('utf-8'))
        data.append(experiment)
        #percent = (sample_index/number_of_samples * 100) 
        #print('Done ' + str(sample_index + 1) + '/' + str(number_of_samples) + ' samples or ' + str(percent) + '% of the set.')
    return data

def savetoCSV(data, filename):
    fields = data[0].keys()
    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = fields)
        writer.writeheader()  
        writer.writerows(data)

def main(): 
    beta_data = pd.read_csv(r'/home/lilia/beebiome-taxonomy-scripts/old-site-data.csv')
    for f_name in os.listdir('/home/lilia/beebiome/beebiome-update/data/Apoidea/'):
        if (f_name.startswith('Apoidea_sra.')):
            print(f_name + '.csv is processing...')
            data = parseXML('/home/lilia/beebiome/beebiome-update/data/Apoidea/' + f_name, beta_data)
            f_name = f_name.replace('.', '_')
            savetoCSV(data, '/home/lilia/beebiome/beebiome-scripts/xml_parse_output_run31May2021/' + f_name + '.csv') 
            df = pd.read_csv(r'/home/lilia/beebiome/beebiome-scripts/xml_parse_output_run31May2021/'+ f_name + '.csv')
            new_df = df.reindex(columns=['sample_accession', 'in_beta', 'experiment_accession', 'experiment_alias', 'experiment_xml', 'submission_lab', 'submission_center',
                                    'submission_accession', 'submission_alias', 'submission_xml', 'organization_type', 'organization_xml',
                                    'study_center', 'study_alias', 'study_accession', 'study_xml', 'study_xml', 'sample_alias', 
                                    'sample_xml', 'pool_xml', 'run_set_xml'])
            new_df.to_csv(r'/home/lilia/beebiome/beebiome-scripts/xml_parse_output_run31May2021/'+ f_name + '.csv', index=False)           

if __name__ == "__main__":
  
    # calling main function
    main()
