#!/usr/bin/env python3
import os
import pickle
import sys
import errno
import argparse
import csv
from math import log, exp
from Bio import SeqIO
import operator
import datetime
import time
from parse_blast_xml import parser


#1. Input files and output directory

def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="BLAST-based viral completeness verification")
    parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-f', required = True, help='Input fasta file')
    required_args.add_argument('-o', required = True, help='Output directory')
    optional_args = parser.add_argument_group('optional arguments')
    optional_args.add_argument('-t', help='Number of threads')   
    optional_args.add_argument('-thr', help='Completeness threshold (0.0-1.0), default = 0.9')   
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


args = parse_args(sys.argv[1:])
base = os.path.basename(args.f)
name_file = os.path.splitext(base)[0]
dirname = os.path.dirname(os.path.abspath(__file__))
outdir = args.o

try:
    os.makedirs(outdir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise
    
name = os.path.join(outdir, name_file)


if args.t:
        threads = str(args.t)
else:
        threads = str(20)

if args.thr:
        threshold = float(args.thr)
else:
        threshold = 0.9


real_len = {} # real length of input contigs
for record in SeqIO.parse(open(args.f),"fasta"):
    real_len[record.id] = len(record)



with open(dirname + "/data/viral_genomes_train.pkl", 'rb') as f:
    genomes = pickle.load(f)

with open(dirname + "/data/viral_genomes_len_train.pkl", 'rb') as f:
    genomes_len = pickle.load(f)

with open(dirname + "/data/proteins_to_genomes_train.pkl", 'rb') as f:
    proteins_to_genomes = pickle.load(f)

with open(dirname + "/data/viral_genomes_train.pkl", 'rb') as f:
    train = pickle.load(f)



#2. Gene prediction and protein extraction using Prodigal

print ("Gene prediction...")
res = os.system ("prodigal -p meta -i " + args.f + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
if res != 0:
   print ("Prodigal run failed")
   exit(1)   


#3. Blast predicted proteins against RefSeq viral proteins from the train dataset
print ("Running BLAST...")
blastdb = dirname + "/blast_db/train_proteins_refseq"
os.system ("blastp  -query " + name+"_proteins.fa" + " -db " + blastdb + " -evalue 1e-06 -outfmt 5 -out "+name+".xml -num_threads "+threads)



#4. Parse blast output

parser(name+".xml", outdir)

#5. For each contig:
#a) We take each protein and get list of viruses corresponding to protein-to-protein hits (from "proteins_to_genomes.pkl" - need to fix).
#b) Then we combine these viral hits into the single set corresponding to each protein, and assign these sets to contig.


def get_blast_hits(blast_xml):
  #Input - blast xml, output - dict of viruses and corresponding viral hits by proteins
  result_handle = open(blast_xml)
  from Bio.Blast import NCBIXML
  blast_records = NCBIXML.parse(result_handle)
  prot_to_viruses = {} 
  for record in blast_records:
      query_id = record.query.split()[0]
      if query_id not in prot_to_viruses:
        prot_to_viruses[query_id] = []
        
      for alignment in record.alignments:
        prot_to_viruses[query_id] += [proteins_to_genomes[alignment.title.split()[1]]]



  virus_to_viruses={}
  for i in prot_to_viruses:
        virus_name = i.split()[0].rsplit("_", 1)[0]
        if virus_name not in virus_to_viruses:
            virus_to_viruses[virus_name] = {}        
        if len(prot_to_viruses[i]) > 0:
          virus_to_viruses[virus_name][i]=set(prot_to_viruses[i])

  return(virus_to_viruses)


queries =  get_blast_hits(name+".xml")
    

#6. Build corresponding probability distribution

eps=0.00001
queries_log_prob = {}


# We take set of RefSeq viruses from training dataset ("train.pkl")
# a) Create dict of probabilities - for each virus from training dataset assign eps or add 1/n where n is num of viral hits for each protein
# b) Turn probabilities into logs - single value for each train virus.

aposteriori = {} 
for k in train:
 aposteriori[k] = {}
 
 for contig in queries:
  aposteriori[k][contig] = []
  for protein in queries[contig]:
   n = len(queries[contig][protein])

   if k in queries[contig][protein]:
    aposteriori[k][contig] += [1/n - eps]
   else:
    aposteriori[k][contig] += [eps*n/(len(train)-n)] 


for contig in queries:    
 p_vir_genes={}
 for virus in aposteriori:
  sum_log = 0
  for prob in aposteriori[virus][contig]:
      sum_log += log(prob)
        
  p_vir_genes[virus] = sum_log
 
 queries_log_prob[contig] = p_vir_genes

# Get argmax, print output
final_table = []
complete_list = []
partial_list = []

for i in queries_log_prob:
  print ("Query: ", i, real_len[i])
  maxValue = max(queries_log_prob[i].items(), key=operator.itemgetter(1))
  compl =  real_len[i]/genomes_len[maxValue[0]][0]
  if compl >= threshold:
     result = "Full-length"
     complete_list.append(i)
  else:
     result = "Partial"
     partial_list.append(i)


  print (result, "{:.1%}".format(compl),  maxValue[0], genomes_len[maxValue[0]][0], genomes_len[maxValue[0]][1])
  final_table.append([i, real_len[i], "{:.1%}".format(compl), result, maxValue[0], genomes_len[maxValue[0]][0], genomes_len[maxValue[0]][1]])  

if not os.path.exists(outdir + "/Prediction_results_fasta/"):
    os.mkdir(outdir + "/Prediction_results_fasta/")
complete_fasta = []
partial_fasta = []
for record in SeqIO.parse(open(args.f),"fasta"):
    if record.id in complete_list:
        complete_fasta.append(record)
    if record.id in partial_list:
        partial_fasta.append(record)

SeqIO.write(complete_fasta, outdir + "/Prediction_results_fasta/complete_viruses.fasta", 'fasta')
SeqIO.write(partial_fasta, outdir + "/Prediction_results_fasta/partial_viruses.fasta", 'fasta') 


result_file = name + "_result_table.csv"
with open(result_file, 'w') as output:
         writer = csv.writer(output, lineterminator='\n')
         writer.writerows(final_table)

