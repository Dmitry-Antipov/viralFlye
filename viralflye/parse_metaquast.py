#!/usr/bin/python
import sys
import os
import random
import subprocess
import shutil
import re
from genericpath import isdir, exists
from os.path import join
from joblib import Parallel, delayed
import viralflye.check_mash



class reference_stats:
    def __init__(self, name):
        self.name = name
        self.largest_alignment = 0    
        self.reference_length = 0
        self.genome_fraction = 0
        self.largest_contig = 0
        self.different_contigs = 0
    def __str__(self):
        return(self.name +" " + str(self.reference_length) + " " + str(self.largest_alignment) + " " + str(self.different_contigs)) 

def parse(work_dir) :
    res = {}
    dir = join(work_dir, "summary", "TSV")
    gf = join(dir, "Genome_fraction.tsv")
    for line in open(gf, "r"):
        arr = line.split()
        if arr[0] != "Assemblies":
            if arr[1] == "-":
                arr[1] = 0
            res[arr[0]] = reference_stats(arr[0])
            res[arr[0]].genome_fraction = float(arr[1])/ 100
    lf = join (dir, "Largest_contig.tsv")
    for line in open(lf, "r"):
        arr = line.split()        
        if arr[1] == "-":
           arr[1] = 0

        if arr[0] != "Assemblies" and arr[0] != "not_aligned":   
            res[arr[0]].largest_contig = int(arr[1])
    af = join (dir, "Largest_alignment.tsv")
#    print(af)
    for line in open(af, "r"):
        arr = line.split()        
        if arr[1] == "-":
            arr[1] = 0
        if arr[0] != "Assemblies" and arr[0] != "not_aligned":   
 #           print(line)
            res[arr[0]].largest_alignment = int(arr[1])
#runs_per_reference/GCA_012935815.1_ASM1293581v1_genomic/report.tsv     
    rpr = join(work_dir, "runs_per_reference")
    for dir in os.listdir(rpr):
        rep = join(rpr, dir, "report.tsv")
        for line in open (rep, "r"):
            arr = line.split("\t")
            if len(arr) == 2 and arr[0] == "Reference length":
                res[dir].reference_length = int(arr[1])


#1 20983 | 46409 25457 | 20983 20953 | 99.69 | GCA_012949105.1_ASM1294910v1_genomic_MN693581.1 contig_8767 | cs:Z::158-a:27-g:278-gg:4*ga:44-a:519*ca:11*gt:214-a:90-t:451-a:108*cg:269-a:410-t:938*ga:27-t:1551-a:311-t:778*ga:16-t:61-a:48*ag:780*ga:691-t:53-t:543-a:335-t:1027-g:244-t:429+a:358-a:57-a:117-a:66-a:1119-a:358-t:597*ag:692-a:395-t:119-a:481+g:113+c:5+t:204-a:21-a:139-t:137+a:369+a:12+t:348-a:252-a:507-t:161-a:8-a:374-t:443-t:1660-a:734*ga:155+tg:55*cg:46+g:121+c:143*ga:113*ta*at:34
    for genome in os.listdir(join(work_dir, "runs_per_reference")):
        coords = join(work_dir, "runs_per_reference", genome, "contigs_reports/minimap_output/assembly.coords.filtered")
        if os.path.isfile(coords):
            contigs = set()
            for line in open (coords, "r"):
                contig = line.split("|")[4].strip().split()[1]
                contigs.add(contig)
            res[genome].different_contigs = len(contigs)
            
#runs_per_reference/GCA_012953035.1_ASM1295303v1_genomic/contigs_reports/minimap_output/assembly.coords.filtered
    return res
def get_stats(work_dir) :
    stats = parse(work_dir)
    genome20 = set()
    genome50 = set()
    genome80 = set()
    alignment80 = set()
    large80 = set()
    for s in stats:
        ref = stats[s]
        if ref.genome_fraction < 0.20:
            genome20.add(s)
        if ref.genome_fraction < 0.50:
            genome50.add(s)
        if ref.genome_fraction < 0.80:
            genome80.add(s)
        else:
            if ref.genome_fraction * ref.reference_length * 0.8 > ref.largest_alignment:
       
                alignment80.add(s)
            if ref.largest_contig *0.8 > ref.reference_length:
                large80.add(s)
    print ("Format: reference_id ref_length largest_alignment different_contigs_in_quast_alignments")
    print ("Total genomes: " + str(len(stats)))
    print("Genomes with less than 0.2 genome fraction: "+ str(len(genome20)))
    for g in genome20:
        print(stats[g])
    print("Genomes with less than 0.5 genome fraction: "+ str(len(genome50)))
    for g in genome50:
        print(stats[g])
    print("Genomes with less than 0.8 genome fraction: " + str(len(genome80)))
    for g in genome80:
        print(stats[g])
    print ("Other problems are counted only for genomes with genome fraction more than 0.8")
    print("Fragmented genomes (largest alignment is less than 0.8 of aligned genome length): " + str(len(alignment80)))
    for g in alignment80:
        print(stats[g])
    print("Extended genomes (largest contig is more than 1/0.8 of reference length): " + str(len(large80)))
    for g in large80:
        print(stats[g])

def main():
    if len(sys.argv) != 2:
        print ("Usage: " + sys.argv[0] + "path to metaquast results for assembly")
        print("--fast --min-alignment 500 --min-identity 80 in metaquast please")
        exit()
    work_dir = sys.argv[1]
    get_stats(work_dir)


if __name__ == '__main__':
     main()

