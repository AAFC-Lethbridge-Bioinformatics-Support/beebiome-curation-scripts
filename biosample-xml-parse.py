import csv
from typing import TYPE_CHECKING
import xmltodict
import os
import re
import requests
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np

matches_to_exclude = ['sp.', 'aff.', 'nr.', 'n.', 'subsp.', 'group sp.', 'x', 'S.']
colors_list = ['blue', 'orange', 'red', 'black', 'white', 'green', 'purple', 'yellow']

def calcStringDistance(taxonomy_xml_key, sliced_key):
    taxonomy_xml_sliced_key = taxonomy_xml_key.replace('(', '').replace(')', '').split()
    taxonomy_xml_sliced_key = [s for s in taxonomy_xml_sliced_key if s not in matches_to_exclude]
    # It might be incorrect to assume order doesn't matter
    length_of_match = 0
    for key_fragment in sliced_key:
        if key_fragment in taxonomy_xml_sliced_key:
            length_of_match = length_of_match + len(key_fragment)
    taxonomy_xml_sliced_key_length = len(''.join(taxonomy_xml_sliced_key))
    return length_of_match / taxonomy_xml_sliced_key_length

def makeTaxonomyDict(root, genera_list):
    tax_dict = {}
    for taxon_index, taxon in enumerate(root.findall('./')):
        # What's the difference between the TaxID attribute vs. the ParentTaxID
        temp_taxid = None
        for field in taxon:
            if field.tag == 'TaxId':
                temp_taxid = field.text
            elif field.tag == 'ScientificName':
                tax_dict[field.text] = temp_taxid
                sliced_scientific_name = field.text.split()
                if sliced_scientific_name[0] == "unclassified":
                    if sliced_scientific_name[1] not in genera_list:
                        genera_list.append(sliced_scientific_name[1])
                elif sliced_scientific_name[0] not in genera_list:
                    genera_list.append(sliced_scientific_name[0])
            elif field.tag == 'OtherNames':
                # This code can be uncommented when I am sure there are only two types of other names
                #for disp_name in field.findall('Name/DispName'):
                #    tax_dict[disp_name.text] = temp_taxid
                #for synonym in field.findall('Synonym'):
                #    tax_dict[synonym.text] = temp_taxid
                for other_name in field:
                    if other_name.tag == 'Name':
                        name_fields = other_name.findall('./DispName')
                        if len(name_fields) > 1:
                            print("Multiple display names?")
                            exit()
                        tax_dict[name_fields[0].text] = temp_taxid
                    # Maybe add extra cases to like handle GenbankCommonName and CommonName
                    # cuz for example, maybe i should add functionality for where if the first
                    # like modifier of a phrase doesn't match such as "red mason bee" vs "mason bee"
                    # it becomes the species. 
                    elif other_name.tag in ['Synonym', 'Includes', 'GenbankCommonName', 'CommonName',
                        'EquivalentName', 'GenbankSynonym', 'BlastName']:
                        tax_dict[other_name.text] = temp_taxid
                   # else:
                    #    print("THERE IS ANOTHER TYPE OF NAME!!!!")
                     #   print("SEE:" + other_name.tag + " OF " + str(temp_taxid) + " TAXID AND #" + str(taxon_index) + " TAXON")
                      #  exit()
            else:
                continue
    return tax_dict

def doesGenusMatch(entry_sliced_key, sliced_key, genera_list):
    for key_fragment in entry_sliced_key:
        if key_fragment in genera_list and key_fragment in sliced_key:
            return True
    return False

def checkApoidea(taxonomy_xml_dict, key, host_subs, genera_list):
    # Assuming doing this is allowed:
    # key_wo_brackets_and_sp = my_string.replace(')', '').replace('(', '').replace(' sp.', '')
    # replace line below with:
    # if key_wo_brackets_and_sp in dict:
    if key in taxonomy_xml_dict:
        return taxonomy_xml_dict[key]
    elif key in host_subs:
        return "Found inexact match"
    else:
        sliced_key = key.split()
        # Perhaps I should consider treating the text inside of brackets differently and not splitting it by spaces
        # instead maybe I should keep the stuff inside brackets as one key_fragment and then look for exact matches 
        # instead
        sliced_key = [s.translate({ord('('): None, ord(')'): None}) for s in sliced_key]
        # Ummm this might not be the best solution because sp. is not the only abbreviated metadata tag there's
        # nr. n. aff. subsp.
        # Host "Licopersicon de solanum" was matched to "Astata graeca de Beaumont, 1965" 2494830
        # ^ Maybe it's worth just ignoring results if the length of the similarity is too small
        sliced_key = [s for s in sliced_key if s not in matches_to_exclude]
        #print(key_split)
        key_contains_genus = False
        for key_fragment_index, key_fragment in enumerate(sliced_key):
            if key_fragment in genera_list:
                key_contains_genus = True

        max_similarity_key = None
        max_similarities = 0

        max_similarity_sliced_key = None

        max_similarity_key_contains_genus = False

        for taxonomy_xml_key in taxonomy_xml_dict:
            taxonomy_xml_sliced_key = taxonomy_xml_key.replace('(', '').replace(')', '').split()
            taxonomy_xml_sliced_key = [s for s in taxonomy_xml_sliced_key if s not in matches_to_exclude]

            taxonomy_xml_key_contains_genus = False
            for taxonomy_xml_key_fragment_index, taxonomy_xml_key_fragment in enumerate(taxonomy_xml_sliced_key):
                if taxonomy_xml_key_fragment in genera_list:
                    taxonomy_xml_key_contains_genus = True
            
            similarities = 0

            for key_fragment in sliced_key:
                # this might be useful to turn into case insensitive
                if key_fragment in taxonomy_xml_sliced_key:
                    similarities = similarities + 1
            # Handle case where no max similarity yet
            if max_similarity_key == None:
                max_similarities = similarities
                max_similarity_key = taxonomy_xml_key
                max_similarity_sliced_key = taxonomy_xml_sliced_key
                max_similarity_key_contains_genus = taxonomy_xml_key_contains_genus
            elif similarities > max_similarities: 
                if taxonomy_xml_key_contains_genus and max_similarity_key_contains_genus:
                    key_genus_match_taxonomy_xml = doesGenusMatch(taxonomy_xml_sliced_key, sliced_key, genera_list)
                    key_genus_match_max_similarity = doesGenusMatch(max_similarity_sliced_key,sliced_key, genera_list)
                    if key_genus_match_taxonomy_xml:
                        max_similarities = similarities
                        max_similarity_key = taxonomy_xml_key
                        max_similarity_sliced_key = taxonomy_xml_sliced_key
                        max_similarity_key_contains_genus = True
                    elif (not key_genus_match_taxonomy_xml and key_genus_match_max_similarity):
                        continue
                    else:
                        max_similarities = similarities
                        max_similarity_key = taxonomy_xml_key
                        max_similarity_sliced_key = taxonomy_xml_sliced_key
                        max_similarity_key_contains_genus = True
                elif not taxonomy_xml_key_contains_genus and max_similarity_key_contains_genus:
                    key_genus_match_taxonomy_xml = doesGenusMatch(taxonomy_xml_sliced_key, sliced_key, genera_list)
                    key_genus_match_max_similarity = doesGenusMatch(max_similarity_sliced_key,sliced_key, genera_list)
                    if key_genus_match_max_similarity:
                        continue
                    else:
                        max_similarities = similarities
                        max_similarity_key = taxonomy_xml_key
                        max_similarity_sliced_key = taxonomy_xml_sliced_key
                        max_similarity_key_contains_genus = True
                else:
                    max_similarities = similarities
                    max_similarity_key = taxonomy_xml_key
                    max_similarity_sliced_key = taxonomy_xml_sliced_key
                    max_similarity_key_contains_genus = True
            elif similarities == max_similarities and similarities != 0:
                #I'm pretty sure this will always take species (just "Megachile") as a result instead of like
                #  Megachile sp. (leafcutter bee) -> Megachile lerma,,,,, if there are not other matches
                # cuz if you think about it part/tiniest whole is bigger than part/larges whole
                key_genus_match_taxonomy_xml = doesGenusMatch(taxonomy_xml_sliced_key, sliced_key, genera_list)
                key_genus_match_max_similarity = doesGenusMatch(max_similarity_sliced_key,sliced_key, genera_list)
                max_worse = (calcStringDistance(max_similarity_key, sliced_key) < 
                calcStringDistance(taxonomy_xml_key, sliced_key))
                if key_genus_match_taxonomy_xml == key_genus_match_max_similarity:
                    if key_genus_match_taxonomy_xml:
                        if max_worse:
                            max_similarities = similarities
                            max_similarity_key = taxonomy_xml_key
                            max_similarity_sliced_key = taxonomy_xml_sliced_key
                            max_similarity_key_contains_genus = True
                        else:
                            continue
                    else:
                        if max_worse:
                            max_similarities = similarities
                            max_similarity_key = taxonomy_xml_key
                            max_similarity_sliced_key = taxonomy_xml_sliced_key
                            max_similarity_key_contains_genus = True
                        else:
                            continue
                else:
                    # if both have different bool values of match upper
                    if key_genus_match_taxonomy_xml:
                        max_similarities = similarities
                        max_similarity_key = taxonomy_xml_key
                        max_similarity_sliced_key = taxonomy_xml_sliced_key
                        max_similarity_key_contains_genus = True
                    else:
                        continue
            else:
                # if max similarities are bigger than new similarities
                key_genus_match_taxonomy_xml = doesGenusMatch(taxonomy_xml_sliced_key, sliced_key, genera_list)
                key_genus_match_max_similarity = doesGenusMatch(max_similarity_sliced_key,sliced_key, genera_list)   
                if key_genus_match_taxonomy_xml and not key_genus_match_max_similarity:
                    max_similarities = similarities
                    max_similarity_key = taxonomy_xml_key
                    max_similarity_sliced_key = taxonomy_xml_sliced_key
                    max_similarity_key_contains_genus = True
                else:
                    continue
        # Then if there are actually no partial matches in the taxonomy file, return "Not in Apoidea"
        if max_similarities == 0:
            return "Not in Apoidea"
        else:
            if max_similarities == 1:
                for term in max_similarity_sliced_key:
                    if term in sliced_key:
                        if term == 'de' or term in colors_list:
                            return "Not in Apoidea"
                        else:
                            break

            string_to_return = str('Host "' + key + '" was matched to "' + max_similarity_key + '" ' + taxonomy_xml_dict[max_similarity_key])
            print(string_to_return)
            host_subs[key] = string_to_return
            return "Found inexact match"

def parseXML(xmlfile, taxonomy_dict, beta_data, auto_no, host_subs, genera_list):
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
            #load if in xml return autoNo
            biosample['in_beta'] = 'Yes'
            biosample['load'] = 'inBeta'
        else:
            biosample['in_beta'] = 'No'
            biosample['load'] = ''

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
                        if str(attribute.text) == 'NA':
                            biosample['host'] = 'not applicable'
                        if re.match(auto_no, biosample['host'].casefold()):
                            biosample['load'] = 'autoNo'
                            biosample['taxid'] = 'Not in Apoidea'
                        else:
                            #biosample['load'] stays equal to ''
                            biosample['taxid'] = checkApoidea(taxonomy_dict, attribute.text, host_subs, genera_list)

                            #LM 03Dec2021, set biosample to autoNo if host is Not in Apoidea
                            if biosample['taxid'] == 'Not in Apoidea':
                                biosample['load'] = 'autoNo'
                        hosted = True
                if hosted == False:
                    # new column host with N/A as the host
                    biosample['host'] = "Not available"
                    biosample['taxid'] = "Not available"
                    biosample['load'] = 'autoNo' #LM 03Dec2021, set biosample to autoNo
                    #biosample['load'] stays equal to ''
                # then all attributes are turned into xml text in utf-8
                #biosample['attributes'] = ET.tostring(child).decode('utf-8')
                biosample['attributes'] = xmltodict.parse(ET.tostring(child).decode('utf-8'))
            elif child.tag == 'Ids':
                for id in child:
                    if id.get("db") == 'SRA':
                        biosample['has_sra'] = id.text
                        break
                else:
                    biosample['has_sra'] = 'no sra'
                    if biosample['load'] == '':
                        biosample['load'] = 'autoNo'  #LM: might be yes on NCBI, check for a fix latter
                biosample[child.tag.lower()] = xmltodict.parse(ET.tostring(child).decode('utf-8'))

            elif child.tag == 'Links':
                for link in child:
                    if link.get("target") == 'bioproject':
                        biosample['has_proj'] = link.get("label") 
                        break
                else:
                    biosample['has_proj'] = 'no bioproject'
                biosample[child.tag.lower()] = xmltodict.parse(ET.tostring(child).decode('utf-8'))
            else:
                # then all other tags are turned into xml test in utf-8
                #biosample[child.tag.lower()] = ET.tostring(child).decode('utf-8')
                biosample[child.tag.lower()] = xmltodict.parse(ET.tostring(child).decode('utf-8'))
        data.append(biosample)
        #print(sample_index)
        #percent = (sample_index/number_of_samples * 100) 
        print('Done ' + str(sample_index + 1) + '/' + str(number_of_samples) + ' samples') # or ' + str(percent) + '% of the set.')
    return data

def savetoCSV(data, filename):
    # specifying the fields for csv file
    fields = ['taxid', 'access', 'pub_date', 'last_update', 'sub_date', 'id', 'accession', 'in_beta', 'load', 'has_sra', 'has_proj', 'ids', 'description', 'owner', 'models', 'package', 'host', 'attributes', 'links', 'status']
    # writing to csv file
    with open(filename, 'w') as csvfile:
        # creating a csv dict writer object
        writer = csv.DictWriter(csvfile, fieldnames = fields)
        # writing headers (field names)
        writer.writeheader()
        # writing data rows
        writer.writerows(data)

def main():
    import sys, os
    
    # LM 03Dec2021, for security set the space to work with the Apoidea folder set as argument
    # Get the Apoidea folder to process
    arguments = len(sys.argv) - 1
    if arguments == 0:
        print("Error, no argumet provided. Run this program with the full path to the Apoidea folder as argument.")
        exit()
    dnIn = sys.argv[1]
    dnIn = os.path.normpath(dnIn)
    print("Selected folder to process: " + dnIn)
    
    # Set the working directory
    dnWork = os.getcwd() #os.path.basename(os.getcwd())  
    dn, dName  = os.path.split(dnWork)

    # Set the output folder
    dnOut = dnIn + '\\xml_parse_output_CSVs\\'
    if not os.path.exists(dnOut):
        os.makedirs(dnOut)

    apoidea_taxonomy_tree = ET.parse(dnIn + '\\Apoidea_taxonomy.xml')
    apoidea_taxonomy_root = apoidea_taxonomy_tree.getroot()
    apoidea_genera = []
    apoidea_taxonomy_dict = makeTaxonomyDict(apoidea_taxonomy_root, apoidea_genera)
    host_subs = pd.read_csv(dnWork + '\\host_subs.csv', header=None, index_col=0, squeeze=True).to_dict()
    #print(apoidea_taxonomy_dict) 
    #print(checkApoidea(apoidea_taxonomy_dict, "Sudila"))
    beta_data = pd.read_csv(dnWork + '\\in-beta.csv')
    auto_no = pd.read_csv(dnWork + '\\auto-no.csv')
    auto_no = [str(x) for x in auto_no['autoNo if host in'].tolist()]
    auto_no = [x.casefold() for x in auto_no]
    auto_no = "(" + ")|(".join(auto_no) + ")"

    #print(auto_no)
    for f_name in os.listdir(dnIn):
        #and (not f_name.startswith('Apoidea_biosample.1.'))
        if (f_name.startswith('Apoidea_biosample.')):
            #print(f_name)
            print(f_name + '.csv is being processed...')
            data = parseXML(dnIn + '\\' + f_name, apoidea_taxonomy_dict, beta_data, auto_no, host_subs, apoidea_genera)
            f_name = f_name.replace('.', '_')
            savetoCSV(data, dnOut + f_name + '.csv')
            df = pd.read_csv(dnOut + f_name + '.csv')
            new_df = df.reindex(columns=['host', 'taxid', 'in_beta', 'load', 'has_sra', 'has_proj', 'access', 'pub_date', 'last_update', 'sub_date', 'id', 'accession', 
            'ids', 'description', 'owner', 'models', 'package', 'attributes', 'links', 'status'])
            new_df.to_csv(dnOut + f_name + '.csv', index=False)
            pd.DataFrame.from_dict(data=host_subs, orient='index').to_csv(dnWork + 'host_subs.csv', header=False)

if __name__ == "__main__":
    # calling main function
    main()
