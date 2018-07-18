#!/usr/bin/python
# Code for annotating Nucleolar RNABP with Intrinsically disordered regions

# Import packages and modules
import argparse
import csv
import d2p2
import protinfo
import numpy as np
from time import gmtime, strftime
import sys
import re

import collections
import imp
imp.reload(protinfo)
imp.reload(d2p2)


# Parsing user input
args = argparse.ArgumentParser(description='This program adds IDRs from certain callers to a uniprot protein list',
epilog="All IDRs should be in the output file")

args.add_argument('-i',dest="infile",type=str,
                   help='Input file containing protein IDs',required=True)
args.add_argument('-o', dest="outfile", type=str,
                   help='Name for output file')
args.add_argument('-c', dest = "consensus",type=int,default=3,
                   help='The minimum caller consensus you require for an IDR to qualify. Default = 3')
args.add_argument('-l',dest = "minlen",type=int,default=20,
                   help='The minimum length of IDR called by all consensus callers for an IDR to qualify. Default = 20')
args.add_argument('-w',dest = "whitelist",nargs='+',default=None,
                   help='The callers you want to include for IDR calling. Default is None which means all callers are included')
args.add_argument('-b',dest = "blacklist",nargs='+',default=None,
                   help='The callers you want to include for IDR calling. Default is None')  

#args.print_help()
argp = args.parse_args()


# We start by reading in the file with Uniprot IDs for Mark and David's protein list
f = csv.reader(open(argp.infile,'rU'),delimiter="\t")
col_names = f.next()
#print(col_names)

# Extracting unipto id column index and storing information with uniprot.id as key 
idcol = col_names.index('uniprot.id')
uniprot_ids = []
metadata = {}

# Save Uniprot IDs as well as metadata 
for x in f:
    uniprot_ids.append(x[idcol])
    metadata[x[idcol]] = x

# Are there any tools in D2P2 that you want or do not want reported
# eg : VL-XT, VSL2b, PrDOS, PV2, Espritz and IUPred

tools_whitelist = argp.whitelist
tools_blacklist = argp.blacklist
consensus_req = argp.consensus # Minimum number of predictors that need to agree on a base being called disordered
min_block_size = argp.minlen # Minimum length of disordered region all of which should be called by at least 'consensus_req' number of disorder predictors
print tools_whitelist, tools_blacklist, consensus_req, min_block_size

# Declare outfile
temp = argp.infile.split("_")[3]

if tools_whitelist:
    outfile = temp+"_c"+str(consensus_req)+"_len"+str(min_block_size)+"_"+'.'.join(tools_whitelist)+"_outfile.txt"
else:
    outfile = temp+"_c"+str(consensus_req)+"_len"+str(min_block_size)+"_all-callers_outfile.txt"
    
outfile = "Python-output/"+outfile
print "Infile : "+argp.infile
print "Outfile : "+outfile

# If user inputs a specific name, use that
if argp.outfile:
    outfile = argp.outfile

outf = open(outfile,'w')
outf.write('\t'.join(col_names)+"\t"+"\t".join(map(str, ("idrs", "consensus","overall.consensus","protein_length","total_idr_length", "fraction_idr", "number.of.idrs","protter.url"))) + "\n")

# Creating IDR data containers to hold the information for proteins that do and do not have IDR annotations
idr_data_available = set()
idr_data_unavailable = set()
sequence_length_idr = {}

blocks = collections.defaultdict(list) # Blocks of consensus calls on disorder
med_consensus = collections.defaultdict(list) # Mean consensus for the blocks of consensus calls
n = 0 # number of ids that we have iterated through

# Looping through the set of uniprot ids provided above and adding disorder annotations
# d2p2_entry is an item of class d2p2 that has the characteristics 

for d2p2_entry in d2p2.iterator(uniprot_ids):
    
    # Build the predictor list that will be used to call consensus on IDR sequence
    d2p2_entry.rebuildConsensus(
        tools_whitelist=tools_whitelist, tools_blacklist=tools_blacklist)
    
    # Set minimum peptide length and consensus support for IDRs to be called
    d2p2_entry.setIDRs(consensus_req, min_block_size)
    
    # Add the IDRs to a list
    #print(d2p2_entry.idrs)
    blocks[d2p2_entry.name] = d2p2_entry.idrs
    
    ku_con = np.array(d2p2_entry.consensus)
    g = d2p2_entry.idrs
    med = []
    
    # Calculate median consensus for each IDR in block and save it in med_consensus
    #for rang in g:
    #    calc = np.median(ku_con[int(rang[0]):int(rang[1])+1])
    #    med.append(float("%.2f" % round(calc,2)))
    #        
    #med_consensus[d2p2_entry.name] = med 
    
    for rang in g:
        #if tools_whitelist is not None:
        calc = np.median(ku_con[int(rang[0]):int(rang[1])+1])
        #print calc
        med.append(float("%.2f" % round(calc,2)))    
        #else:
        #    med.append(1)
        med_consensus[d2p2_entry.name] = med
    
    # Add successful entries to data available
    idr_data_available.add(d2p2_entry.name)
    sequence_length_idr[d2p2_entry.name] = d2p2_entry.length
    
    # Checking how long it takes to process requests
    n += 1
    if n%500 == 0:
        print('proteins done: %i %s' % (
                n, strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    

# Add all uniprot ids without IDRs to idr_data_unavailable
idr_data_unavailable = set(uniprot_ids).difference(idr_data_available)

# Check output 
print "IDR data available for "+str(len(idr_data_available))+" proteins." # ID available
print "IDR data missing for "+str(len(idr_data_unavailable))+" proteins." # ID unavailable


# Prints results to outfile
for protein in idr_data_available:
    total_idr = 0
    protein_length=0
    fraction_idr = 0
    for block in blocks[protein]:
        total_idr += (block[1] - block[0])
        protein_length = sequence_length_idr[protein]
        fraction_idr = total_idr/float(protein_length)
    idr_ranges = ','.join(map(str, ('-'.join(str(e)for e in x) for x in blocks[protein])))
    medcon = ",".join(str(e) for e in med_consensus[protein])
    
    # Trying to reconstruct the url for Protter using idr, biotin and PTMs
    tw = tools_whitelist
    if tw is None:
    	tw = ['VLXT','VSL2b','PrDOS','PV2','IUPred-S','IUPred-L','Espritz-N','Espritz-X','Espritz-D']
    gene_name = metadata[protein][col_names.index('uniprot.id')]+" : "+metadata[protein][col_names.index('gene.names')]+" : "+metadata[protein][col_names.index('protein.names')]+" : " + str(",".join(tw))
    #print(gene_name+"\n\n")
    protter_url = "http://wlab.ethz.ch/protter/create?up="+protein+"&tm=auto&mc=lightsalmon&lc=blue&tml=numcount&numbers&legend&title="+gene_name+"&n:IDRs,fc:red,bc:peachpuff="+idr_ranges+"&n:Phosphorylation,s:box,cc:aqua,fc:midnightblue,bc:cornflowerblue="+metadata[protein][col_names.index('Phosphorylation')]+"&n:Ubiqutination,s:diamond,cc:black,fc:olive,bc:greenyellow="+metadata[protein][col_names.index('Ubiqutination')]+"&n:Acetylation,s:box,cc:lightpink,fc:lightpink,bc:deeppink="+metadata[protein][col_names.index('Acetylation')]+"&n:Biotins,s:circ,cc:white,fc:forestgreen,bc:forestgreen="+metadata[protein][col_names.index('loc.hits')]+"&format=pdf"
    protter_url = re.sub("=NA&n:","=&n:", protter_url)   
    #print(protter_url)
    
    
    # Write results to outfile
    outf.write(("\t".join(metadata[protein])+"\t"+"\t".join(map(str, (idr_ranges,medcon,np.median(med_consensus[protein]),protein_length, total_idr, fraction_idr, len(blocks[protein]),protter_url)))) + "\n")
    
outf.close()

             
# Run a test protein using d2p2 database
protein = idr_data_available.pop()
bl = blocks[protein]
print("IDRs for " + protein+ " are : "+','.join(map(str, ('-'.join(str(e)for e in x) for x in bl))))
print("Length of " +protein+" : "+str(sequence_length_idr[protein]))

# Takes 3:30 mins for ~3000 proteins