import csv
import xmltodict
import os
import requests
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np


def parseXML(xmlfile):
    #parses the read-in xmlfile and populates a list of dicts
    tree = ET.parse(xmlfile)
    tax_tree = ET.parse('/home/lilia/beebiome/beebiome-update/data/Apoidea/Apoidea_taxonomy.xml')

    root = tree.getroot()
    tax_root = tax_tree.getroot()
    # print(root) ----> BioSampleSet at some memory location
    # you can make a list to like collect stuff
    data = []
    number_of_samples = len(tax_root)
    for sample_index, item in enumerate(root.findall('./')):
        biosample = {}
        # print(item) ---> BioSample at some memory location

        # you can collate this stuff into a dictionary/array and then return it
        # you can also use dict = item.attrib to get a dictionary of all the xml 'attributes'
        #print(item.get("access"))
        biosample['access'] = item.get("access")
        #print(item.get("publication_date"))
        biosample['pub_date'] = item.get("publication_date")
        #print(item.get("last_update"))
        biosample['last_update'] = item.get("last_update")
        #print(item.get("submission_date"))
        biosample['sub_date'] = item.get("submission_date")
        #print(item.get("id"))
        biosample['id'] = item.get("id")
        #print(item.get("accession"))
        biosample['accession'] = item.get("accession")

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
                        for taxon in tax_root.findall('./'):
                            #print(taxon.tag)
                            temp_taxid = 'temp meaningless value'
                            for field in taxon:
                                #print(field.tag)
                                if field.tag == 'TaxId':
                                    temp_taxid = field.text
                                    #print(temp_taxid)
                                if field.tag == 'ScientificName':
                                    #print(attribute.text.lower())
                                    if field.text.lower() == attribute.text.lower():
                                        #print(temp_taxid + ' ' + field.text)
                                        biosample['taxid'] = temp_taxid
                                        #print(temp_taxid)
                                        break
                            else:
                                biosample['taxid'] = "Not in Apoidea"
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
    fields = ['taxid', 'access', 'pub_date', 'last_update', 'sub_date', 'id', 'accession', 'ids', 'description', 'owner', 'models', 'package', 'host', 'attributes', 'links', 'status']
  
    # writing to csv file
    with open(filename, 'w') as csvfile:
  
        # creating a csv dict writer object
        writer = csv.DictWriter(csvfile, fieldnames = fields)
  
        # writing headers (field names)
        writer.writeheader()
  
        # writing data rows
        writer.writerows(data)

def main(): 
    #compiled_data = [] 
    for f_name in os.listdir('/home/lilia/beebiome/beebiome-update/data/Apoidea/'):
        if f_name.startswith('Apoidea_biosample.'):
            data = parseXML('/home/lilia/beebiome/beebiome-update/data/Apoidea/' + f_name)
            savetoCSV(data, '/home/lilia/beebiome/beebiome-scripts/xml-parse-output/' + f_name + '.csv')
            df = pd.read_csv(r'/home/lilia/beebiome/beebiome-scripts/xml-parse-output/'+ f_name + '.csv')
            new_df = df.reindex(columns=['host', 'taxid', 'access', 'pub_date', 'last_update', 'sub_date', 'id', 'accession', 
            'ids', 'description', 'owner', 'models', 'package', 'attributes', 'links', 'status'])
            new_df.to_csv(r'/home/lilia/beebiome/beebiome-scripts/xml-parse-output/'+ f_name + '.csv', index=False)
            print(f_name + '.csv is done being processed.')
            #compiled_data.extend(parseXML('/home/lilia/beebiome/beebiome-update/data/Apoidea/' + f_name))
    #print(compiled_data)
    #savetoCSV(compiled_data, 'test.csv')

if __name__ == "__main__":
  
    # calling main function
    main()
