#!/usr/bin/python
# Code for retrieving PDFs based on a URL

# Import packages and modules
import argparse
import os
import urllib
import requests
import csv
import sys
import re


# Parsing user input
args = argparse.ArgumentParser(description='This program retrieves pdfs based on a url',
epilog="All PDFs should be in the output folder")

args.add_argument('-i',dest="infile",type=str,
                   help='Input file containing protein IDs and URLs',required=True) 

#args.print_help()
argp = args.parse_args()

# Declare output directory based on file name
dpath = argp.infile.split("_")[0].split("/")[1]+"_"+argp.infile.split("_")[3]
print dpath
directory = "/Users/manasa/Documents/Work/TTT/13_Biopep_DMinde/"+dpath
print directory

if not os.path.exists(directory):
    os.makedirs(directory)


# We start by reading in the file with Uniprot IDs and URLs for IDR regions by protein
f = csv.reader(open(argp.infile,'rU'),delimiter="\t")
col_names = f.next()
#print(col_names)

# Extracting uniprot id and url column index and storing information 
idcol = col_names.index('uniprot.id')
genecol = col_names.index('gene.names')
urlcol = col_names.index('protter.url')


# Save Uniprot IDs as well as metadata 
for x in f:
    url = x[urlcol]
    outfile = directory+"/"+dpath+"_"+x[idcol]+"_"+x[genecol]+".pdf"
    print(outfile)
    urllib.urlretrieve(url, outfile)


'''
# Takes 3:30 mins for ~3000 proteins
'''