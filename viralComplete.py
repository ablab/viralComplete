import os
import pickle
import sys
import csv
from math import log, exp
from Bio import SeqIO
import operator
import datetime
import time
from parse_blast_xml import parser
blastdb = "/Bmo/ncbi_nt_database/nt"

### to add: argparse

#1. Input files and output directory

base = os.path.basename(sys.argv[1])
name_file = os.path.splitext(base)[0]

print(os.getcwd()+"/"+sys.argv[2])
if not os.path.exists(os.getcwd()+"/"+sys.argv[2]):
    os.mkdir(os.getcwd()+"/"+sys.argv[2])
name = os.path.join(os.getcwd()+"/"+sys.argv[2],  name_file)


#if not os.path.exists(os.getcwd()+"/"+name_file+"_viralComplete_NB_out"):
#    os.mkdir(os.getcwd()+"/"+name_file+"_viralComplete_NB_out")
#name = os.path.join(os.getcwd()+"/"+name_file+"_viralComplete_NB_out",  name_file)


#name = os.path.splitext(os.path.basename(sys.argv[1]))[0]
threads = str(20)


real_len = {} # real length of input contigs
for record in SeqIO.parse(open(sys.argv[1]),"fasta"):
    real_len[record.id] = len(record)



with open(os.path.dirname(os.path.abspath(__file__)) + "/viral_genomes_train.pkl", 'rb') as f:
    genomes = pickle.load(f)

with open(os.path.dirname(os.path.abspath(__file__)) + "/viral_genomes_len_train.pkl", 'rb') as f:
    genomes_len = pickle.load(f)

with open(os.path.dirname(os.path.abspath(__file__)) + "/proteins_to_genomes_train.pkl", 'rb') as f:
    proteins_to_genomes = pickle.load(f)
    

with open(os.path.dirname(os.path.abspath(__file__)) + "/viral_genomes_train.pkl", 'rb') as f:
    train = pickle.load(f)



#2. Gene prediction and protein extraction using Prodigal

print ("Gene prediction...")
res = os.system ("prodigal -p meta -i " + sys.argv[1] + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
if res != 0:
   print ("Prodigal run failed")
   exit(1)   


#3. Blast predicted proteins against RefSeq viral proteins from the train dataset
print ("Running BLAST...")
blastdb = "/Nancy/mrayko/viruses/viralComplete/train_proteins_refseq"
os.system ("blastp  -query " + name+"_proteins.fa" + " -db " + blastdb + " -evalue 1e-06 -outfmt 5 -out "+name+".xml -num_threads "+threads)

  # + " -num_alignments 5" )


#4. Parse blast output
#parser(name+".xml", os.getcwd()+"/"+name+"_viralComplete_NB_out")
#print(name)
#os.mkdir(os.getcwd()+"_viralComplete_NB_out")
#parser(name+".xml", os.getcwd()+"_viralComplete_NB_out")
parser(name+".xml", sys.argv[2])


#5. For each contig:
#a) We take each protein and get list of viruses corresponding to protein-to-protein hits (from "proteins_to_genomes.pkl" - need to fix).
#b) Then we combine these viral hits into the single set corresponding to each protein, and assign these sets to contig.

### To fix - take into account number of hits and e-value/bitscore?

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
  print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
  print("first part done")


  virus_to_viruses={}
  for i in prot_to_viruses:
        virus_name = i.split()[0].rsplit("_", 1)[0]
        if virus_name not in virus_to_viruses:
            virus_to_viruses[virus_name] = {}        
        if len(prot_to_viruses[i]) > 0:
          virus_to_viruses[virus_name][i]=set(prot_to_viruses[i])
  print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
  print("second part done")

  return(virus_to_viruses)


print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
queries =  get_blast_hits(name+".xml")
    

#6. Build corresponding probability distribution

eps=0.00001
queries_log_prob = {}


# We take set of RefSeq viruses from training dataset ("train.pkl")
# a) Create dict of probabilities - for each virus from training dataset assign eps or add 1/n where n is num of viral hits for each protein
# b) Turn probabilities into logs - single value for each train virus.

print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
print("first part done")




#for contig in queries:
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




#   for virus in queries[contig][protein]:
 #   aposteriori[virus][contig] += [1/n - eps] 

# aposteriori = {}
# for k in train:
 #   aposteriori[k] = []
 #for protein in queries[contig]:
 #  n = len(queries[contig][protein])
#   print (contig +  "No hits:" + str(eps*n/(len(train)-n)))
#   print (contig +  "Hit:" + str(1/n - eps))

  # for k in train:
  #      if k in queries[contig][protein]:
   #         aposteriori[k] += [1/n - eps]
            #print(aposteriori[k])

    #    else:
     #       aposteriori[k] += [eps*n/(len(train)-n)] 

print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
print("third part done")

# print(aposteriori)

for contig in queries:    
 p_vir_genes={}
 for virus in aposteriori:
  sum_log = 0
  for prob in aposteriori[virus][contig]:
  #  sum_log = 0
   # for prob in aposteriori[i]:
      sum_log += log(prob)
        
  p_vir_genes[virus] = sum_log
 
 queries_log_prob[contig] = p_vir_genes

print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
print("fourth part done")

# Get argmax, print output
final_table = []

for i in queries_log_prob:
  print ("Query: ", i, real_len[i])
  #maxValue = max(queries_log_prob[i].items(), key=operator.itemgetter(1))[1]
  #print ([key for key in queries_log_prob[i].keys() if queries_log_prob[i][key]==maxValue])
  #print ([genomes_len[key][0] for key in queries_log_prob[i].keys() if queries_log_prob[i][key]==maxValue])
  maxValue = max(queries_log_prob[i].items(), key=operator.itemgetter(1))
  if genomes_len[maxValue[0]][0]<= real_len[i]/0.9 and genomes_len[maxValue[0]][0]>= 0.9*real_len[i]:  


    print ("Full-length", maxValue[0], genomes_len[maxValue[0]][0], genomes_len[maxValue[0]][1])
    final_table.append([i, real_len[i], "Full-length", maxValue[0], genomes_len[maxValue[0]][0], genomes_len[maxValue[0]][1]])
  else:
    print ("Partial", maxValue[0], genomes_len[maxValue[0]])
    final_table.append([i, real_len[i], "Partial", maxValue[0], genomes_len[maxValue[0]][0], genomes_len[maxValue[0]][1]])
  


  #print (maxValue[0], maxValue[1], genomes_len[maxValue[0]])

#          final_table.append([k[0], genomes_len[k[0]][0], genomes_len[k[0]][1], k[1]])
  

result_file = name + "_result_table.csv"
with open(result_file, 'w') as output:
         writer = csv.writer(output, lineterminator='\n')
         writer.writerows(final_table)

