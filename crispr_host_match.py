#!/usr/bin/python3
import sys
import os
from Bio import SeqIO


import argparse
parser = argparse.ArgumentParser(description=f'''Script for host predictions based on CRISPR spacers.  
It takes viralFlye result as an input, predicts CRISPR spacers and matches them with viruses using BLAST.
''',
formatter_class=argparse.RawTextHelpFormatter)
parser._action_groups.pop()
required_args = parser.add_argument_group('required arguments')
required_args.add_argument('-f', required = True, help='viralFlye output directory')
required_args.add_argument('-o', required = True, help='Output directory')
            
    
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args(sys.argv[1:])

flye = args.f
outdir = args.o

os.system(f"mkdir {outdir}")

print ("Keep contigs over 1 kb...")
os.system (f"seqtk seq -L 1000 {flye}/assembly.fasta  > {outdir}/assembly_1k.fasta")

print ("Running minced...")
os.system (f"minced  -spacers  {outdir}/assembly_1k.fasta {outdir}/minced_out")


print ("Collecting viruses...")
os.system (f"cat {flye}/linears_viralFlye.fasta  {flye}/components_viralFlye.fasta {flye}/circulars_viralFlye.fasta > {outdir}/viruses.fa")


print ("Blasting...")

os.system (f"makeblastdb -dbtype nucl -in {outdir}/viruses.fa -out {outdir}/viruses_db")
os.system (f"blastn -query {outdir}/minced_out_spacers.fa -db {outdir}/viruses_db -evalue 1E-5 -outfmt 6  -out {outdir}/blast.out -max_target_seqs 5 -num_threads 10  -max_hsps 5 -perc_identity 0.9  -task blastn-short")

print(f"Verification results can be found in {outdir}/blast.out")
