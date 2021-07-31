# beebiome-taxonomy-scripts
Some scripts which are being used to parse the BioSampleSet .xmls.

## Script methodology

1. Read in all the XML files in which the BioSample data is stored.
2. Create a dataframe where each column is a field in the BioSample.
3. Read in the ‘host name’ field of each BioSample -- labelled as ‘host’ in the ‘Attributes’ child element of the BioSample XML.
4. If the ‘host name’ is on a word list of ‘auto-not-in-Apoidea-supergroup’, remove this BioSample from the dataframe.
5. Check if the ‘host name’ provided in the field has any exact matches with any entries in the NCBI Taxonomy Browser. If it does, fill in the taxid with the corresponding one from the Taxonomy database.
6. If there are no exact matches, start looking for inexact matches.
  * Process the ‘host name’ to not have any illegal strings (e.g. ‘sp.’, ‘!’, ‘aff.’, any colors).
  * Find the entry in the NCBI Taxonomy Browser which has the most and longest (in length) matches to the processed ‘taxid’. 
  * Highlight these entries with inexact matches for review (tag them with the suggested taxid).
7. For each ‘host name’ with no inexact matches, remove this BioSample from the dataframe.
8. The highlighted inexact matches are removed as well from the dataframe -- but kept in a special file for future review.
9. Output the remaining entries in the dataframe to a .csv which can be funneled into the website backend database.




