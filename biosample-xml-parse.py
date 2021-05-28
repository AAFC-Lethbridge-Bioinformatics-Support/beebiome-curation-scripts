import csv
import xmltodict
import os
import requests
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np

def makeTaxonomyDict(root):
    tax_dict = {}
    for taxon in root.findall('./'):
        temp_taxid = 'empty'
        for field in taxon:
            if field.tag == 'TaxId':
                temp_taxid = field.text
            elif field.tag == 'ScientificName':
                tax_dict[field.text] = temp_taxid
            else:
                continue
    return tax_dict

# Function to print sum
def checkApoidea(dict, key):
    if key in dict:
        return dict[key]
    else:
        return "Not in Apoidea"

def parseXML(xmlfile, taxonomy_dict, beta_data):
    #parses the read-in xmlfile and populates a list of dicts
    tree = ET.parse(xmlfile)

    root = tree.getroot()

    # print(root) ----> BioSampleSet at some memory location
    # you can make a list to like collect stuff
    data = []
    number_of_samples = len(root)
    for sample_index, item in enumerate(root.findall('./')):
        biosample = {}
        # print(item) ---> BioSample at some memory location

        # you can collate this stuff into a dictionary/array and then return it
        # you can also use dict = item.attrib to get a dictionary of all the xml 'attributes'
        #print(item.get("access"))
        biosample['access'] = item.get("access")
        biosample['pub_date'] = item.get("publication_date")
        biosample['last_update'] = item.get("last_update")
        biosample['sub_date'] = item.get("submission_date")
        biosample['id'] = item.get("id")
        biosample['accession'] = item.get("accession")
        if biosample['accession'] in beta_data.loc[:, "BioSample acc"].values:
            biosample['in_beta'] = 'Yes'
        else:
            biosample['in_beta'] = 'No'

        for child in item:
            # print(child) ---> for ex Element 'Attributes'
            if child.tag == 'Attributes':
                hosted = False
                for attribute in child:
                    # attributes need to be treated differently, that is we need to check for 
                    # the attribute tag with attribute name host
                    if attribute.get("attribute_name") == 'host':
                        # new column host with the host name
                        biosample['host'] = attribute.text
                        biosample['taxid'] = checkApoidea(taxonomy_dict, attribute.text)
                        hosted = True
                if hosted == False:
                    # new column host with N/A as the host
                    biosample['host'] = "Not avaliable"
                    biosample['taxid'] = "Not avaliable"
                # then all attributes are turned into xml text in utf-8
                #biosample['attributes'] = ET.tostring(child).decode('utf-8')
                biosample['attributes'] = xmltodict.parse(ET.tostring(child).decode('utf-8'))
            else:
                # then all other tags are turned into xml test in utf-8
                #biosample[child.tag.lower()] = ET.tostring(child).decode('utf-8')
                biosample[child.tag.lower()] = xmltodict.parse(ET.tostring(child).decode('utf-8'))
        data.append(biosample)
        #percent = (sample_index/number_of_samples * 100) 
        #print('Done ' + str(sample_index + 1) + '/' + str(number_of_samples) + ' samples or ' + str(percent) + '% of the set.')
    return data

def savetoCSV(data, filename):
  
    # specifying the fields for csv file
    fields = ['taxid', 'access', 'pub_date', 'last_update', 'sub_date', 'id', 'accession', 'in_beta', 'ids', 'description', 'owner', 'models', 'package', 'host', 'attributes', 'links', 'status']
  
    # writing to csv file
    with open(filename, 'w') as csvfile:
  
        # creating a csv dict writer object
        writer = csv.DictWriter(csvfile, fieldnames = fields)
  
        # writing headers (field names)
        writer.writeheader()
  
        # writing data rows
        writer.writerows(data)

def main(): 
    apoidea_taxonomy_tree = ET.parse('/home/lilia/beebiome/beebiome-update/data/Apoidea(copy)/Apoidea_taxonomy.xml')
    apoidea_taxonomy_root = apoidea_taxonomy_tree.getroot()
    apoidea_taxonomy_dict = makeTaxonomyDict(apoidea_taxonomy_root)
    #print(apoidea_taxonomy_dict)
    #print(checkApoidea(apoidea_taxonomy_dict, "Sudila"))
    beta_data = pd.read_csv(r'/home/lilia/beebiome-taxonomy-scripts/old-site-data.csv')
    for f_name in os.listdir('/home/lilia/beebiome/beebiome-update/data/Apoidea(copy)/'):
        if (f_name.startswith('Apoidea_biosample.') and (not f_name.startswith('Apoidea_biosample.2.'))):
            print(f_name + '.csv is being processed...')
            data = parseXML('/home/lilia/beebiome/beebiome-update/data/Apoidea(copy)/' + f_name, apoidea_taxonomy_dict, beta_data)
            savetoCSV(data, '/home/lilia/beebiome/beebiome-scripts/results (march)/' + f_name + '.csv')
            df = pd.read_csv(r'/home/lilia/beebiome/beebiome-scripts/results (march)/'+ f_name + '.csv')
            new_df = df.reindex(columns=['host', 'taxid', 'in_beta', 'access', 'pub_date', 'last_update', 'sub_date', 'id', 'accession', 
            'ids', 'description', 'owner', 'models', 'package', 'attributes', 'links', 'status'])
            new_df.to_csv(r'/home/lilia/beebiome/beebiome-scripts/results (march)/'+ f_name + '.csv', index=False)
            #compiled_data.extend(parseXML('/home/lilia/beebiome/beebiome-update/data/Apoidea/' + f_name))
    #print(compiled_data)
    #savetoCSV(compiled_data, 'test.csv')

if __name__ == "__main__":
  
    # calling main function
    main()
