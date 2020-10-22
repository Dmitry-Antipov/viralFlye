#!/usr/bin/python3
import sys
import os

if len(sys.argv) != 3:
    print (f"Usage: {sys.argv[0]} <assembly.stats file from metaflye assembly> <flye/spades> <outdir>")
    exit()


genome = sys.argv[1]
#viruses = sys.argv[2]
assembly_type = sys.argv[2]
outdir = sys.argv[3]


print ("Keep contigs over 1 kb...")
    res = os.system (f"/Nancy/mrayko/Libs/seqtk/seqtk seq -L 1000 {genome}  > {outidr}/assembly_1k.fasta")
    if res != 0:
        print ("seqtk run failed")
        exit(1)   


print ("Running minced...")
    res = os.system (f" /Nancy/mrayko/Libs/minced/minced  -spacers  {outdir}/assembly_1k.fasta minced_out")
    if res != 0:
        print ("Minced run failed")
        exit(1)   


print ("Selecting cirular contigs...")


print ("Predicting viruses...")


print ("Blasting...")

os.system (f"makeblastdb -dbtype nucl -in {outdir}/viruses.fa -out {outdir}/viruses_db")
os.system (f"blastn -query {outdir}/minced_out_spacers..fa -db {outdir}/viruses_db -evalue 1E-5 -outfmt 6  -out {outdir}/blast.out -max_target_seqs 3 -num_threads 10  -max_hsps 5 -perc_identity 0.9  -task blastn-short")



