import csv
import xmltodict
import os
import requests
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np


def parseXML(xmlfile):
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
    for f_name in os.listdir('/home/lilia/beebiome/beebiome-update/data/Apoidea/'):
        if (f_name.startswith('Apoidea_sra.') and (not f_name.startswith('Apoidea_sra.2.'))):
            print(f_name + '.csv is processing...')
            data = parseXML('/home/lilia/beebiome/beebiome-update/data/Apoidea/' + f_name)
            savetoCSV(data, '/home/lilia/beebiome/beebiome-scripts/xml-parse-output/' + f_name + '.csv')            

if __name__ == "__main__":
  
    # calling main function
    main()
