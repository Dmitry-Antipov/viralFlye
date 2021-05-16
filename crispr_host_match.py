#!/usr/bin/python3
import sys
import os
from Bio import SeqIO

if len(sys.argv) != 4:
    print (f"Usage: {sys.argv[0]} <flye_assembly_folder> <spades_assembly_folder> <outdir>")
    exit()


flye = sys.argv[1]
spades = sys.argv[2]
outdir = sys.argv[3]

os.system(f"mkdir {outdir}")

def is_circular_knp(seq, min_match=50, max_match=200):
    def prefix(s):
        v = [0] * len(s)
        for i in range(1, len(s)):
            k = v[i - 1]
            while k > 0 and s[k] != s[i]:
                k = v[k - 1]
            if s[k] == s[i]:
                k += 1
            v[i] = k
        return v
    s = seq[:max_match] + "$" + seq[-max_match:]
    ps = prefix(s)
    return ps[-1] >= min_match


print ("Keep contigs over 1 kb...")
os.system (f"seqtk seq -L 1000 {flye}/assembly.fasta  > {outdir}/flye_assembly_1k.fasta")
os.system (f"seqtk seq -L 1000 {spades}/scaffolds.fasta > {outdir}/spades_assembly_1k.fasta")


print ("Running minced...")
os.system (f"minced  -spacers  {outdir}/flye_assembly_1k.fasta {outdir}/flye_minced_out")
os.system (f"minced  -spacers  {outdir}/spades_assembly_1k.fasta {outdir}/spades_minced_out")


print ("Selecting cirular contigs...")

os.system(f"{os.path.dirname(sys.argv[0])}/extract_circulars_metaflye.py {flye}/assembly_info.txt > {outdir}/flye_circular.names") 

# Flye
contigs = SeqIO.parse(open(outdir + "/flye_assembly_1k.fasta"), 'fasta')
circular_contigs = []
circular_names = set()
for line in open(f"{outdir}/flye_circular.names"):
  circular_names.add(line.strip())
for contig in contigs:
  seq_id = str(contig.id)
  if seq_id in circular_names:
    circular_contigs.append(contig)
SeqIO.write(circular_contigs, open(outdir+ "/flye_circular.fasta", mode='w'), 'fasta')

# Illlumina
contigs = SeqIO.parse(open(outdir + "/spades_assembly_1k.fasta"), 'fasta')
circular_contigs = []
for contig in contigs:
    seq = contig.seq
    if len(seq) >= 500 and is_circular_knp(seq):
        circular_contigs.append(contig)
SeqIO.write(circular_contigs, open(outdir+ "/spades_circular.fasta", mode='w'), 'fasta')


print ("Selecting isolated components")
os.system(f"python3 {os.path.dirname(sys.argv[0])}/extract_components_from_gfa.py {flye}/assembly_graph.gfa  > {outdir}/flye_isolated.fasta")
os.system(f"python3 {os.path.dirname(sys.argv[0])}/extract_components_from_gfa.py {spades}/assembly_graph_with_scaffolds.gfa  > {outdir}/spades_isolated.fasta")
os.system(f"cat {outdir}/flye_circular.fasta {outdir}/flye_isolated.fasta > {outdir}/flye_merged.fasta")
os.system(f"cat {outdir}/spades_circular.fasta {outdir}/spades_isolated.fasta > {outdir}/spades_merged.fasta")

# Intersection

os.system(f"makeblastdb -dbtype nucl -in {outdir}/spades_merged.fasta -out {outdir}/spades_circ_db")
os.system(f"blastn -query {outdir}/flye_merged.fasta -db {outdir}/spades_circ_db -outfmt \"6 qseqid sseqid qlen slen length pident evalue qcovs qcovhsp\" -out {outdir}/blast_circ.out")

intersect_contigs = []
intersect_names = set()
for line in open(f"{outdir}/blast_circ.out"):
    print(line.split())
    if int(line.split()[7]) > 90:
        intersect_names.add(line.split()[1]) 
print(intersect_names)
for contig in circular_contigs:
  print(contig.id)
  seq_id = str(contig.id)
  if seq_id in intersect_names:
    intersect_contigs.append(contig)
SeqIO.write(intersect_contigs, open(f"{outdir}/intersect.fasta", mode='w'), 'fasta')


print ("Predicting viruses...")

os.system(f"viralverify.py -f {outdir}/intersect.fasta -p -o {outdir}/viralverify --hmm /Nancy/mrayko/db/pfam/Pfam-A.hmm -t 15")
os.system(f"viralverify.py -f {outdir}/flye_merged.fasta -p -o {outdir}/viralverify_flye --hmm /Nancy/mrayko/db/pfam/Pfam-A.hmm -t 15")
os.system(f"viralverify.py -f {outdir}/spades_merged.fasta -p -o {outdir}/viralverify_spades --hmm /Nancy/mrayko/db/pfam/Pfam-A.hmm -t 15")

print ("Blasting...")

os.system (f"makeblastdb -dbtype nucl -in {outdir}/viralverify/Prediction_results_fasta/intersect_virus.fasta -out {outdir}/viruses_db")
os.system (f"blastn -query {outdir}/flye_minced_out_spacers.fa -db {outdir}/viruses_db -evalue 1E-5 -outfmt 6  -out {outdir}/flye_blast.out -max_target_seqs 3 -num_threads 10  -max_hsps 5 -perc_identity 0.9  -task blastn-short")
os.system (f"blastn -query {outdir}/spades_minced_out_spacers.fa -db {outdir}/viruses_db -evalue 1E-5 -outfmt 6  -out {outdir}/spades_blast.out -max_target_seqs 3 -num_threads 10  -max_hsps 5 -perc_identity 0.9  -task blastn-short")

