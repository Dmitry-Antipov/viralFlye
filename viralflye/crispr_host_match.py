#!/usr/bin/python3
import sys
import os
from Bio import SeqIO


import argparse
parser = argparse.ArgumentParser(description=f'''Script for host predictions based on CRISPR spacers.
It takes metaFlye result as an input, extracts circular and linear isolated contigs,
predicts viruses and CRISPR spacers and matches them using BLAST.
''',
formatter_class=argparse.RawTextHelpFormatter)
parser._action_groups.pop()
required_args = parser.add_argument_group('required arguments')
required_args.add_argument('-f', required = True, help='metaFlye output directory')
required_args.add_argument('-o', required = True, help='Output directory')
required_args.add_argument('--hmm', help='Path to Pfam-A HMM database for viralVerify script')    
            
    
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args(sys.argv[1:])

flye = args.f
outdir = args.o
hmm = args.hmm

os.system(f"mkdir {outdir}")



print ("Keep contigs over 1 kb...")
os.system (f"seqtk seq -L 1000 {flye}/assembly.fasta  > {outdir}/assembly_1k.fasta")

print ("Running minced...")
os.system (f"minced  -spacers  {outdir}/assembly_1k.fasta {outdir}/minced_out")

print ("Selecting cirular contigs...")

os.system(f"{os.path.dirname(os.path.abspath(__file__))}/extract_circulars_metaflye.py {flye}/assembly_info.txt > {outdir}/circular_names.txt") 

contigs = SeqIO.parse(open(outdir + "/assembly_1k.fasta"), 'fasta')
circular_contigs = []
circular_names = set()
for line in open(f"{outdir}/circular_names.txt"):
  circular_names.add(line.strip())
for contig in contigs:
  seq_id = str(contig.id)
  if seq_id in circular_names:
    circular_contigs.append(contig)
SeqIO.write(circular_contigs, open(outdir+ "/circular.fasta", mode='w'), 'fasta')


print ("Selecting isolated components")
os.system(f"python3 {os.path.dirname(os.path.abspath(__file__))}/extract_components_from_gfa.py {flye}/assembly_graph.gfa > {outdir}/isolated.fasta")
os.system(f"cat {outdir}/circular.fasta {outdir}/isolated.fasta > {outdir}/merged.fasta")


print ("Predicting viruses...")

os.system(f"viralverify.py  -f {outdir}/merged.fasta -p -o {outdir}/viralverify --hmm {hmm} -t 15")

print ("Blasting...")

os.system (f"makeblastdb -dbtype nucl -in {outdir}/viralverify/Prediction_results_fasta/merged_virus.fasta -out {outdir}/viruses_db")
os.system (f"blastn -query {outdir}/minced_out_spacers.fa -db {outdir}/viruses_db -evalue 1E-5 -outfmt 6  -out {outdir}/blast.out -max_target_seqs 5 -num_threads 10  -max_hsps 5 -perc_identity 0.9  -task blastn-short")

print(f"Verification results can be found in {outdir}/blast.out")
