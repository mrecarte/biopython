import urllib.request
import sys; print(sys.version)
import platform; print(platform.python_implementation()); print(platform.platform())
import Bio; print(Bio.__version__)
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import MissingExternalDependencyError
from Bio import BiopythonWarning

import pandas as pd

import urllib
import urllib.parse
from urllib.request import urlopen

from xmltodict import parse

# import the fasta file
my_query = SeqIO.read("test.fasta", format="fasta")

print("Query:")

print("\n")

print(my_query)

print("\n")

print("Starting BLAST query...")

result_handle = NCBIWWW.qblast("blastp", "swissprot", my_query.seq, expect=200000, 
                               hitlist_size=500, gapcosts='9 1', matrix_name="PAM30", 
                               filter='F',genetic_code=1, word_size=2, threshold=11, 
                               short_query=True, entrez_query='human [ORGN]', 
                               service="phi", phi_pattern="L-[RKN]-G-G") 
print("\n")
print("BLAST query complete. Saving results.")

# save the BLAST results in an xml file
blast_result = open("my_blast50_swiss.xml", "w")
blast_result.write(result_handle.read())
blast_result.close()
result_handle.close()

# open the results from above

result_handle = open("my_blast50_swiss.xml")

blast_records = NCBIXML.parse(result_handle)

blast_records = list(blast_records)

# this is some code I wrote to parse the xml file and get relevant values
# it's not very nice code! it should probably be a function
print("\n")
print("Organizing BLAST results and preparing for Uniprot search.")

print("\n")

name_list = []
accession_list = []
align_list = []
bit_list = []
e_list = []
frame_list = []
gaps_list = []
identities_list = []
match_list = []
num_align_list = []
pos_list = []
score_list = []
q_start_list = []
q_end_list = []
s_start_list = []
s_end_list = []
GO_list = []
function_list = []

for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            name_list.append(alignment.title)
            accession_list.append(alignment.accession)
            align_list.append(hsp.align_length)
            bit_list.append(hsp.bits)
            e_list.append(hsp.expect)
            frame_list.append(hsp.frame)
            gaps_list.append(hsp.gaps)
            identities_list.append(hsp.identities)
            match_list.append(hsp.match)
            num_align_list.append(hsp.num_alignments)
            pos_list.append(hsp.positives)
            score_list.append(hsp.score)
            q_start_list.append(hsp.query_start)
            q_end_list.append(hsp.query_end)
            s_start_list.append(hsp.sbjct_start)
            s_end_list.append(hsp.sbjct_end)

url = 'https://www.uniprot.org/uploadlists/'

full_protein_name_list = []
function_list = []

print("Starting UniProt query.")

print("\n")

print("Protein names and functions from UniProt:")

print("\n")

for i in range(len(accession_list)):
    
    params = {
        'from':'ACC+ID',
        'to':'ACC',
        'format':'xml', 
        'query': accession_list[i] #iterate through all accession numbers in accession_list
    }
    data = urllib.parse.urlencode(params).encode("utf-8")
    request = urllib.request.Request(url)
    request.add_header('User-Agent', 'Python %s')
    response = urlopen(request, data)
    page = response.read()
    page = page.decode()
    # parse the xml and convert to dict
    uniprot_dict = parse(page)
    # get the full protein name from the dict
    try: 
        if isinstance(uniprot_dict['uniprot']['entry']['protein']['recommendedName']['fullName'], str):
            print(uniprot_dict['uniprot']['entry']['protein']['recommendedName']['fullName'])
            full_protein_name_list.append(uniprot_dict['uniprot']['entry']['protein']['recommendedName']['fullName'])
        else:
            print(uniprot_dict['uniprot']['entry']['protein']['recommendedName']['fullName']['#text'])
            full_protein_name_list.append(uniprot_dict['uniprot']['entry']['protein']['recommendedName']['fullName']['#text'])
    except KeyError:
        print("key error from name")
        full_protein_name_list.append("null")
    # get the protein function from the dict
    try:
        print('\n')
        if isinstance(uniprot_dict['uniprot']['entry']['comment'][0]['text'], str):
            print(uniprot_dict['uniprot']['entry']['comment'][0]['text'])
            function_list.append(uniprot_dict['uniprot']['entry']['comment'][0]['text'])
        else:
            print(uniprot_dict['uniprot']['entry']['comment'][0]['text']['#text'])
            function_list.append(uniprot_dict['uniprot']['entry']['comment'][0]['text']['#text'])
    except KeyError:
        print("key error from comment")
        function_list.append("null")

    print('________')

print("UniProt search complete.")

print("\n")

print("Getting GO terms:")

print("\n")


# get the GO terms for each hit

GO_list = []
for i in range(len(accession_list)):
        
    params = {
        'from':'ACC+ID',
        'to':'ACC',
        'format':'xml', #this is where you can indicate file type; xml and txt work for me
        'query': accession_list[i] #this is the first accession number; once it works I can send all accession numbers
    }
    data = urllib.parse.urlencode(params).encode("utf-8")
    request = urllib.request.Request(url)
    request.add_header('User-Agent', 'Python %s')
    response = urlopen(request, data)
    page = response.read()
    page = page.decode()
    # use above code to parse page
    uniprot_dict = parse(page)

    # this code works for getting GO annotations!
    GO_list_each = []

    for ref in range(len(uniprot_dict['uniprot']['entry']['dbReference'])):
        #print(uniprot_dict['uniprot']['entry']['dbReference'][i]['@type'])
        if uniprot_dict['uniprot']['entry']['dbReference'][ref]['@type'] == 'GO':
            GO_list_each.append(uniprot_dict['uniprot']['entry']['dbReference'][ref]['@id'])
            print(uniprot_dict['uniprot']['entry']['dbReference'][ref]['@id'])
    GO_list.append(GO_list_each)

df = pd.DataFrame(
    {'title': name_list,
    'accession': accession_list,
    'align_length': align_list,
    'bits': bit_list,
    'expect': e_list,
    'frame': frame_list,
    'gaps': gaps_list,
    'identities': identities_list,
    'match': match_list,
    'num_alignments': num_align_list,
    'positives': pos_list,
    'score': score_list,
     'query_start': q_start_list,
     'query_end': q_end_list,
     'sbjct_start': s_start_list,
     'sbjct_end': s_end_list,
     'protein_name': full_protein_name_list,
     'function': function_list,
     'gene_ontology' : GO_list
     
     
    })

print(df.head())

print('\n')

print("removing entries with gaps")

print('\n')

df_no_gaps = df[df['gaps'] == 0]



print("sorting according to bit score")

print('\n')

df_no_gaps.sort_values(by=['bits'])

print("Saving results to csv file.")

df_no_gaps.to_csv("results_uniprot.csv")

import pandas as pd
import numpy as np

import bioservices
from bioservices import QuickGO
g = QuickGO(verbose=False)
go = QuickGO()

df = pd.read_csv("results_uniprot.csv")

# intitialize empty lists
GO_list = []
protein_list = []
res_list = []

# iterate through columns for uniprot results csv

for i,m in zip(df['gene_ontology'], df['protein_name']):

    for j in eval(i):
        GO_list.append(j)
        protein_list.append(m)

# get full GO entries from QuickGO

for annotation in GO_list:
    res = go.get_go_terms(annotation)
    print("\n\n")
    print(annotation)
    print(res)
    res_list.append(res)

print(len(res_list))
print(len(GO_list))
print(len(protein_list))

# create dataframe from lists

output_df = pd.DataFrame(
    {"protein_name": protein_list,
    "GO_number": GO_list,
    "GO_entry": res_list
    })

# save output to csv

output_df.to_csv("GO_entries_0315.csv")

'''
Biological Process describes biological goals accomplished by one or more ordered assemblies of 
molecular functions. A biological process is not equivalent to a pathway. Specifically it does not 
represent any of the dynamics or dependencies that would be required to describe a pathway.
'''

for item in res_list:
    for content in item:
        for values in content.values():
            if values == "biological_process":
                print(content['id'],':',content['definition']['text'])

'''
Molecular Function describes activities, such as catalytic, binding or transporter 
activities, at the molecular level. GO molecular function terms describe activities rather 
than the entities (complexes, gene products or molecules) that perform the actions. 
Typically direct assays such as enzyme kinetics measurements or binding studies can be 
used to infer molecular function annotations. In addition sequence comparison methods are often used to predict 
the molecular function of a gene product because functions are often associated with conserved protein domains.
'''

for item in res_list:
    for content in item:
        for values in content.values():
            if values == "molecular_function":
                print(content['id'],':',content['definition']['text'])

'''
Cellular Component describes locations, at the levels of subcellular structures and macromolecular 
complexes. Experiments informing Cellular Component annotations include fluorescence microscopy and 
co-fractionation of complex members.
'''

for item in res_list:
    for content in item:
        for values in content.values():
            if values == "cellular_component":
                print(content['id'],':',content['definition']['text'])
